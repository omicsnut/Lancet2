#include "lancet/cbdg/dot_renderer.h"

#include "lancet/base/rev_comp.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/dot_layers.h"
#include "lancet/cbdg/dot_overlay_factories.h"
#include "lancet/cbdg/dot_plan.h"
#include "lancet/cbdg/dot_walk_layers.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/probe_tracker.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

namespace {

using namespace std::string_view_literals;

// ── DOT Graph Preamble ──────────────────────────────────────────────────
// Visual style preserved verbatim from the legacy serialize_dot.cpp so the
// new renderer is byte-compatible with existing rendered assets at the
// `subgraph` boundary; only the per-node and per-edge attribute strings
// change to reflect layered overlays.
constexpr auto DOT_PREAMBLE = R"raw(strict digraph G {
graph [layout=neato,bgcolor=black,size="120,180",ratio=compress,rankdir=LR,overlap=vpsc,overlap_shrink=true,start=self];
node [style=filled,fontsize=2,width=2,height=2,fixedsize=false];
edge [color=gray,fontsize=8,fontcolor=floralwhite,len=3,fixedsize=false,headclip=true,tailclip=true];
)raw"sv;

constexpr auto DOT_FOOTER = "}\n}\n"sv;

// ── Role-based default fill colors (stack at z=0; layers override) ─────
constexpr auto COLOR_SHARED = "steelblue"sv;
constexpr auto COLOR_CASE_ONLY = "indianred"sv;
constexpr auto COLOR_CTRL_ONLY = "mediumseagreen"sv;
constexpr auto COLOR_DEFAULT = "lightblue"sv;

// ── Layer convention ───────────────────────────────────────────────────
constexpr auto ANCHOR_BORDER_COLOR = "goldenrod"sv;
constexpr auto PROBE_STRIPE_ACCENT = "orchid"sv;

// ── REF walk accent (walk index 0). ALT walks consume the LAB palette ──
constexpr auto REF_WALK_COLOR = "#FFFFFF"sv;

// ── Maximally-Distinct ALT Walk Color Palette ──────────────────────────
//
// 64 colors pre-computed by scripts/gen_walk_palette.py via k-means
// clustering in CIE L*a*b* (perceptually uniform) color space, then
// reordered by farthest-first traversal so any prefix has maximal pairwise
// separation. Quality (Delta-E = Euclidean L*a*b* distance, >10 = easily
// distinct): K=2: dE=191.7, K=8: 47.3, K=16: 29.8, K=64: 12.6.
// Walk indices beyond 64 cycle through the palette via modulo.
constexpr std::array<std::string_view, 64> ALT_WALK_COLORS = {
    // clang-format off
    "#32C124"sv, "#7D28D9"sv, "#EB3B15"sv, "#2CB1E2"sv, "#70641A"sv, "#DE5E99"sv, "#379B71"sv, "#3C56A8"sv,
    "#E59222"sv, "#EA2CC7"sv, "#8CAB27"sv, "#DC736C"sv, "#BA82E0"sv, "#EB2363"sv, "#A42B19"sv, "#2A6F1A"sv,
    "#953690"sv, "#444EDB"sv, "#DD995B"sv, "#792FAC"sv, "#8C491E"sv, "#8AAF5F"sv, "#21BCC9"sv, "#D16633"sv,
    "#A12D3D"sv, "#D462DF"sv, "#9C3762"sv, "#37BE55"sv, "#4487E6"sv, "#A6943C"sv, "#E3259E"sv, "#938FDB"sv,
    "#BB78B5"sv, "#AC721D"sv, "#724B97"sv, "#409F51"sv, "#3761D5"sv, "#DD2DE6"sv, "#CDA022"sv, "#657F25"sv,
    "#379620"sv, "#3371B2"sv, "#DF297D"sv, "#AC227A"sv, "#A05AE5"sv, "#70B926"sv, "#B322AE"sv, "#B1AA24"sv,
    "#E22743"sv, "#E07B99"sv, "#469FE4"sv, "#E27ED1"sv, "#3ABF81"sv, "#31B9A2"sv, "#D94C3F"sv, "#836FDD"sv,
    "#E66E1C"sv, "#B97449"sv, "#2B703A"sv, "#403AE3"sv, "#DF54B9"sv, "#AA2BD7"sv, "#D23D60"sv, "#D92120"sv,
    // clang-format on
};

auto SignChar(Kmer::Sign sign) -> char {
  return sign == Kmer::Sign::PLUS ? '+' : '-';
}

// ── Role-derived default node style. Layers can override any field ─────
auto RoleStyle(Node const& node) -> NodeStyle {
  if (node.IsShared()) return NodeStyle{.mFillColor = COLOR_SHARED};
  if (node.IsCaseOnly()) return NodeStyle{.mFillColor = COLOR_CASE_ONLY};
  if (node.IsCtrlOnly()) return NodeStyle{.mFillColor = COLOR_CTRL_ONLY};
  return NodeStyle{.mFillColor = COLOR_DEFAULT};
}

// ── Per-edge accumulated overlay (built once before rendering) ─────────
// Mirrors the colorList composition rule: every layer's color (in z-order)
// is appended; the final colorList is rendered as `color="A:B:C"`.
struct OverlayEdge {
  // ── 8B Align ──────────────────────────────────────────────────────────
  std::vector<std::string_view> mColors;

  // ── 4B Align ──────────────────────────────────────────────────────────
  u32 mPenWidth = 0;

  // ── 1B Align ──────────────────────────────────────────────────────────
  bool mDashed = false;
};

// Sort `layers` indirectly into ascending z_order, returning pointer indices.
template <class Layer>
auto SortLayersByZOrder(absl::Span<Layer const> layers) -> std::vector<Layer const*> {
  std::vector<Layer const*> ptrs;
  ptrs.reserve(layers.size());
  for (auto const& layer : layers) ptrs.push_back(&layer);
  std::ranges::sort(ptrs, {}, [](Layer const* layer_ptr) { return layer_ptr->mZOrder; });
  return ptrs;
}

// Compose all node layers into a per-NodeID resolved style. Each non-empty
// field from a higher-z layer overrides earlier ones; the role layer is
// added as a default in `RenderNodeDot`, not here.
auto BuildNodeOverlays(absl::Span<NodeLayer const> layers)
    -> absl::flat_hash_map<NodeID, NodeStyle> {
  absl::flat_hash_map<NodeID, NodeStyle> result;
  for (auto const* layer : SortLayersByZOrder(layers)) {
    for (NodeID const node_id : layer->mIds) {
      auto& slot = result[node_id];
      if (!layer->mStyle.mFillColor.empty()) slot.mFillColor = layer->mStyle.mFillColor;
      if (!layer->mStyle.mBorderColor.empty()) slot.mBorderColor = layer->mStyle.mBorderColor;
      if (!layer->mStyle.mStripeAccent.empty()) slot.mStripeAccent = layer->mStyle.mStripeAccent;
      if (layer->mStyle.mPenWidth != 0) slot.mPenWidth = layer->mStyle.mPenWidth;
      if (layer->mStyle.mPeripheries != 0) slot.mPeripheries = layer->mStyle.mPeripheries;
    }
  }
  return result;
}

// Compose all edge layers into per-LogicalEdge color lists. The colors
// vector preserves z-order so the rendered colorList stripes are stable
// across runs. PenWidth from the highest-z layer that sets it wins.
auto BuildEdgeOverlays(absl::Span<EdgeLayer const> layers)
    -> absl::flat_hash_map<LogicalEdge, OverlayEdge> {
  absl::flat_hash_map<LogicalEdge, OverlayEdge> result;
  for (auto const* layer : SortLayersByZOrder(layers)) {
    for (LogicalEdge const& key : layer->mEdges) {
      auto& slot = result[key];
      if (!layer->mStyle.mColor.empty()) slot.mColors.push_back(layer->mStyle.mColor);
      if (layer->mStyle.mPenWidth != 0) slot.mPenWidth = layer->mStyle.mPenWidth;
      if (layer->mStyle.mDashed) slot.mDashed = true;
    }
  }
  return result;
}

// Compose role default + layer overrides into a single resolved NodeStyle
// without mutating the input. Pulled out of RenderNodeDot so the per-node
// hot path reads as a sequence of named, well-typed steps.
auto ResolveNodeStyle(Node const& node, NodeID node_id,
                      absl::flat_hash_map<NodeID, NodeStyle> const& node_overlays) -> NodeStyle {
  NodeStyle resolved_style = RoleStyle(node);
  auto const overlay_iter = node_overlays.find(node_id);
  if (overlay_iter == node_overlays.end()) return resolved_style;

  NodeStyle const& overlay_style = overlay_iter->second;
  if (!overlay_style.mFillColor.empty()) resolved_style.mFillColor = overlay_style.mFillColor;
  if (!overlay_style.mBorderColor.empty()) resolved_style.mBorderColor = overlay_style.mBorderColor;
  if (!overlay_style.mStripeAccent.empty()) {
    resolved_style.mStripeAccent = overlay_style.mStripeAccent;
  }
  if (overlay_style.mPenWidth != 0) resolved_style.mPenWidth = overlay_style.mPenWidth;
  if (overlay_style.mPeripheries != 0) resolved_style.mPeripheries = overlay_style.mPeripheries;
  return resolved_style;
}

// Fill the per-thread complement buffer with the base-by-base complement of
// `canonical_seq`. The DOT label displays this as the second strand under
// the canonical sequence. Mathematically `reverse(RevComp(seq)) ==
// Complement(seq)` so we walk forward and complement each base directly,
// avoiding the two-string allocation dance the legacy code used.
void FillComplementBuffer(std::string_view canonical_seq, std::string& complement_buffer) {
  complement_buffer.clear();
  complement_buffer.reserve(canonical_seq.size());
  for (char const dna_base : canonical_seq) {
    complement_buffer.push_back(lancet::base::RevComp(dna_base));
  }
}

void RenderNodeDot(fmt::memory_buffer& dot_buffer, Node const& node, NodeID node_id,
                   absl::flat_hash_map<NodeID, NodeStyle> const& node_overlays) {
  auto const resolved_style = ResolveNodeStyle(node, node_id, node_overlays);

  // Per-thread reusable buffer for the complement-strand display string.
  // Single amortised allocation per worker thread for the entire run; no
  // per-node heap allocation on the hot path.
  static thread_local std::string complement_strand_buffer;
  auto const canonical_seq_view = node.SeqView();
  FillComplementBuffer(canonical_seq_view, complement_strand_buffer);

  auto const default_orientation_sign = SignChar(node.SignFor(Kmer::Ordering::DEFAULT));
  auto const node_seq_length = node.Length();
  auto const node_total_read_support = node.TotalReadSupport();

  // Required attributes branch on stripe accent: with-accent emits the
  // colon-paired fillcolor and the "filled,striped" style; without-accent
  // emits the plain fillcolor and "filled" style. Two paths, no
  // intermediate `std::string fill` allocation either way.
  if (resolved_style.mStripeAccent.empty()) {
    fmt::format_to(std::back_inserter(dot_buffer),
                   R"({} [shape=circle fillcolor="{}" style="filled")"
                   R"( label="{}\n{}\n {}:{}\nlength={}\ncoverage={}")",
                   node_id, resolved_style.mFillColor, canonical_seq_view, complement_strand_buffer,
                   node_id, default_orientation_sign, node_seq_length, node_total_read_support);
  } else {
    fmt::format_to(std::back_inserter(dot_buffer),
                   R"({} [shape=circle fillcolor="{}:{}" style="filled,striped")"
                   R"( label="{}\n{}\n {}:{}\nlength={}\ncoverage={}")",
                   node_id, resolved_style.mFillColor, resolved_style.mStripeAccent,
                   canonical_seq_view, complement_strand_buffer, node_id, default_orientation_sign,
                   node_seq_length, node_total_read_support);
  }

  // Optional attributes: border color, penwidth, peripheries.
  if (!resolved_style.mBorderColor.empty()) {
    fmt::format_to(std::back_inserter(dot_buffer), R"( color="{}")", resolved_style.mBorderColor);
  }
  if (resolved_style.mPenWidth != 0) {
    fmt::format_to(std::back_inserter(dot_buffer), " penwidth={}", resolved_style.mPenWidth);
  }
  if (resolved_style.mPeripheries != 0) {
    fmt::format_to(std::back_inserter(dot_buffer), " peripheries={}", resolved_style.mPeripheries);
  }
  fmt::format_to(std::back_inserter(dot_buffer), "]\n");
}

void RenderBaseEdgeDot(fmt::memory_buffer& dot_buffer, Edge const& conn) {
  // Base edges carry sign-pair labels and inherit `color=gray` from the
  // graph-level edge default. Overlay statements emitted in the second
  // pass merge into these via DOT strict-mode semantics, supplying color
  // and penwidth overrides without touching the labels.
  auto const src_sign_char = SignChar(conn.SrcSign());
  auto const dst_sign_char = SignChar(conn.DstSign());
  fmt::format_to(std::back_inserter(dot_buffer),
                 R"({} -> {} [taillabel="{}" headlabel="{}"])"
                 "\n",
                 conn.SrcId(), conn.DstId(), src_sign_char, dst_sign_char);
}

void RenderOverlayEdgeStatement(fmt::memory_buffer& dot_buffer, NodeID src_node_id,
                                NodeID dst_node_id, OverlayEdge const& overlay) {
  if (overlay.mColors.empty()) return;

  auto const color_list_string = absl::StrJoin(overlay.mColors, ":");
  fmt::format_to(std::back_inserter(dot_buffer), R"({} -> {} [color="{}")", src_node_id,
                 dst_node_id, color_list_string);
  if (overlay.mPenWidth != 0) {
    fmt::format_to(std::back_inserter(dot_buffer), " penwidth={}", overlay.mPenWidth);
  }
  if (overlay.mDashed) {
    fmt::format_to(std::back_inserter(dot_buffer), R"( style="dashed")");
  }
  fmt::format_to(std::back_inserter(dot_buffer), "]\n");
}

void RenderToBuffer(fmt::memory_buffer& dot_buffer,
                    absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                    DotPlan const& plan) {
  fmt::format_to(std::back_inserter(dot_buffer), "{}", DOT_PREAMBLE);
  fmt::format_to(std::back_inserter(dot_buffer), "subgraph {} {{\n", plan.mSubgraphName);

  auto const node_overlays = BuildNodeOverlays(plan.mNodeLayers);
  auto const edge_overlays = BuildEdgeOverlays(plan.mEdgeLayers);

  // First pass: nodes + base edges.
  for (auto const& [node_id, node_ptr] : graph) {
    if (node_ptr->GetComponentId() != plan.mCompId) continue;
    RenderNodeDot(dot_buffer, *node_ptr, node_id, node_overlays);
    for (Edge const& conn : *node_ptr) {
      RenderBaseEdgeDot(dot_buffer, conn);
    }
  }

  // Second pass: overlay edge statements, expanded to BOTH DOT directions
  // so anti-parallel mirror splines under `neato` carry the same style.
  for (auto const& [logical_edge_key, overlay] : edge_overlays) {
    RenderOverlayEdgeStatement(dot_buffer, logical_edge_key.mLo, logical_edge_key.mHi, overlay);
    if (logical_edge_key.mLo != logical_edge_key.mHi) {
      RenderOverlayEdgeStatement(dot_buffer, logical_edge_key.mHi, logical_edge_key.mLo, overlay);
    }
  }

  fmt::format_to(std::back_inserter(dot_buffer), "{}", DOT_FOOTER);
}

}  // namespace

auto SerializeToDotString(GraphNodeTable const& graph, DotPlan const& plan) -> std::string {
  fmt::memory_buffer dot_buffer;
  RenderToBuffer(dot_buffer, graph, plan);
  return fmt::to_string(dot_buffer);
}

// ============================================================================
// Walk-layer factory and palette accessor (declared in dot_walk_layers.h)
// ============================================================================

auto WalkColor(usize const walk_index) -> std::string_view {
  if (walk_index == 0) return REF_WALK_COLOR;
  return ALT_WALK_COLORS[(walk_index - 1) % ALT_WALK_COLORS.size()];
}

auto MakeWalkLayers(absl::Span<std::vector<Edge> const> walks) -> std::vector<EdgeLayer> {
  std::vector<EdgeLayer> layers;
  layers.reserve(walks.size());

  // z-order assigns later walks a higher z so their colors append later in
  // the colorList — all walks remain visible thanks to graphviz's parallel-
  // stripe rendering, but the ordering is stable run-to-run.
  for (usize idx = 0; idx < walks.size(); ++idx) {
    EdgeLayer layer;
    layer.mName = idx == 0 ? "walk-ref"sv : "walk-alt"sv;
    layer.mZOrder = 30 + static_cast<i32>(idx);
    layer.mStyle.mColor = WalkColor(idx);
    layer.mStyle.mPenWidth = 2;
    for (Edge const& conn : walks[idx]) {
      layer.mEdges.insert(LogicalEdge::FromEdge(conn));
    }
    layers.emplace_back(std::move(layer));
  }
  return layers;
}

// ============================================================================
// Overlay-layer factories (declared in dot_overlay_factories.h)
// ============================================================================

auto MakeAnchorLayer(NodeIDPair const& source_and_sink) -> NodeLayer {
  NodeLayer layer;
  layer.mName = "anchors"sv;
  layer.mZOrder = 10;
  layer.mIds.insert(source_and_sink[0]);
  layer.mIds.insert(source_and_sink[1]);
  layer.mStyle.mBorderColor = ANCHOR_BORDER_COLOR;
  layer.mStyle.mPenWidth = 3;
  layer.mStyle.mPeripheries = 2;
  return layer;
}

auto MakeProbeLayer(ProbeTracker const& tracker, usize const comp_id,
                    ProbeTracker::NodeTable const& nodes) -> NodeLayer {
  NodeLayer layer;
  layer.mName = "probe"sv;
  layer.mZOrder = 20;
  layer.mIds = tracker.GetHighlightNodeIds(nodes, comp_id);
  layer.mStyle.mStripeAccent = PROBE_STRIPE_ACCENT;
  return layer;
}

// ============================================================================
// ReconstructRefWalk (declared in dot_walk_layers.h)
//
// Walks `ref_node_ids` left-to-right keeping only canonical k-mer IDs that
// survive in this component, deduplicates consecutive duplicates, then for
// each consecutive pair looks up the connecting `Edge` (carrying EdgeKind
// + sign pair the renderer needs). Returns an empty walk when any pair is
// missing — a fragmented backbone yields no overlay rather than a partial
// one to avoid misleading the reader.
// ============================================================================
auto ReconstructRefWalk(absl::Span<NodeID const> ref_node_ids,
                        absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& nodes,
                        usize const comp_id) -> std::vector<Edge> {
  // First pass: filter ref_node_ids down to the IDs that survived pruning
  // and still belong to the requested component. Collapse consecutive
  // duplicates so a single backbone position spanning a self-loop or
  // hairpin doesn't force two entries.
  std::vector<NodeID> surviving_ref_ids;
  surviving_ref_ids.reserve(ref_node_ids.size());
  for (NodeID const ref_id : ref_node_ids) {
    auto const node_iter = nodes.find(ref_id);
    if (node_iter == nodes.end()) continue;
    if (node_iter->second->GetComponentId() != comp_id) continue;
    auto const is_consecutive_duplicate =
        !surviving_ref_ids.empty() && surviving_ref_ids.back() == ref_id;
    if (is_consecutive_duplicate) continue;
    surviving_ref_ids.push_back(ref_id);
  }

  std::vector<Edge> reconstructed_walk;
  if (surviving_ref_ids.size() < 2) return reconstructed_walk;
  reconstructed_walk.reserve(surviving_ref_ids.size() - 1);

  // Second pass: for each consecutive (curr_id, next_id) pair, locate the outgoing edge
  // from curr_id whose destination is next_id. If no such edge exists, the backbone
  // has fragmented and we return an empty walk — a partial walk would mislead the
  // renderer about which edges are backbone-confirmed.
  for (usize ref_idx = 0; ref_idx + 1 < surviving_ref_ids.size(); ++ref_idx) {
    auto const curr_node_id = surviving_ref_ids[ref_idx];
    auto const next_node_id = surviving_ref_ids[ref_idx + 1];

    auto const curr_node_iter = nodes.find(curr_node_id);
    if (curr_node_iter == nodes.end()) return {};
    auto const& curr_node_outgoing_edges = *curr_node_iter->second;

    auto const matches_next_node = [next_node_id](Edge const& candidate_edge) {
      return candidate_edge.DstId() == next_node_id;
    };
    auto const* const edge_to_next_node =
        std::ranges::find_if(curr_node_outgoing_edges, matches_next_node);

    if (edge_to_next_node == curr_node_outgoing_edges.end()) {
      return {};  // backbone fragmented; abandon
    }
    reconstructed_walk.push_back(*edge_to_next_node);
  }

  return reconstructed_walk;
}

}  // namespace lancet::cbdg
