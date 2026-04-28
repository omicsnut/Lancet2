#include "lancet/cbdg/dot_renderer.h"

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
#include "spdlog/fmt/bundled/ostream.h"

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <ranges>
#include <sstream>
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
    "#32C124"sv, "#7D28D9"sv, "#EB3B15"sv, "#2CB1E2"sv, "#70641A"sv, "#DE5E99"sv, "#379B71"sv,
    "#3C56A8"sv, "#E59222"sv, "#EA2CC7"sv, "#8CAB27"sv, "#DC736C"sv, "#BA82E0"sv, "#EB2363"sv,
    "#A42B19"sv, "#2A6F1A"sv, "#953690"sv, "#444EDB"sv, "#DD995B"sv, "#792FAC"sv, "#8C491E"sv,
    "#8AAF5F"sv, "#21BCC9"sv, "#D16633"sv, "#A12D3D"sv, "#D462DF"sv, "#9C3762"sv, "#37BE55"sv,
    "#4487E6"sv, "#A6943C"sv, "#E3259E"sv, "#938FDB"sv, "#BB78B5"sv, "#AC721D"sv, "#724B97"sv,
    "#409F51"sv, "#3761D5"sv, "#DD2DE6"sv, "#CDA022"sv, "#657F25"sv, "#379620"sv, "#3371B2"sv,
    "#DF297D"sv, "#AC227A"sv, "#A05AE5"sv, "#70B926"sv, "#B322AE"sv, "#B1AA24"sv, "#E22743"sv,
    "#E07B99"sv, "#469FE4"sv, "#E27ED1"sv, "#3ABF81"sv, "#31B9A2"sv, "#D94C3F"sv, "#836FDD"sv,
    "#E66E1C"sv, "#B97449"sv, "#2B703A"sv, "#403AE3"sv, "#DF54B9"sv, "#AA2BD7"sv, "#D23D60"sv,
    "#D92120"sv,
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

void RenderNodeDot(std::ostream& out, Node const& node, NodeID node_id,
                   absl::flat_hash_map<NodeID, NodeStyle> const& overlays) {
  // Compose: role default + layer overrides.
  NodeStyle style = RoleStyle(node);
  if (auto const itr = overlays.find(node_id); itr != overlays.end()) {
    NodeStyle const& over = itr->second;
    if (!over.mFillColor.empty()) style.mFillColor = over.mFillColor;
    if (!over.mBorderColor.empty()) style.mBorderColor = over.mBorderColor;
    if (!over.mStripeAccent.empty()) style.mStripeAccent = over.mStripeAccent;
    if (over.mPenWidth != 0) style.mPenWidth = over.mPenWidth;
    if (over.mPeripheries != 0) style.mPeripheries = over.mPeripheries;
  }

  // Build attribute string. fillcolor uses a colon-pair when a stripe
  // accent is set; style toggles between "filled" and "filled,striped".
  std::string fill;
  if (style.mStripeAccent.empty()) {
    fill = std::string(style.mFillColor);
  } else {
    fill = fmt::format("{}:{}", style.mFillColor, style.mStripeAccent);
  }
  std::string_view const style_attr = style.mStripeAccent.empty() ? "filled"sv : "filled,striped"sv;

  auto const dflt_seq = node.SequenceFor(Kmer::Ordering::DEFAULT);
  auto const oppo_seq = node.SequenceFor(Kmer::Ordering::OPPOSITE);
  auto const rev_oppo_seq = std::string(oppo_seq.crbegin(), oppo_seq.crend());

  // Required attributes: shape, fillcolor, style, label.
  fmt::print(out,
             R"({} [shape=circle fillcolor="{}" style="{}")"
             R"( label="{}\n{}\n {}:{}\nlength={}\ncoverage={}")",
             node_id, fill, style_attr, dflt_seq, rev_oppo_seq, node_id,
             SignChar(node.SignFor(Kmer::Ordering::DEFAULT)), node.Length(),
             node.TotalReadSupport());

  // Optional attributes: border color, penwidth, peripheries.
  if (!style.mBorderColor.empty()) fmt::print(out, R"( color="{}")", style.mBorderColor);
  if (style.mPenWidth != 0) fmt::print(out, " penwidth={}", style.mPenWidth);
  if (style.mPeripheries != 0) fmt::print(out, " peripheries={}", style.mPeripheries);
  out << "]\n";
}

void RenderBaseEdgeDot(std::ostream& out, Edge const& conn) {
  // Base edges carry sign-pair labels and inherit `color=gray` from the
  // graph-level edge default. Overlay statements emitted in the second
  // pass merge into these via DOT strict-mode semantics, supplying color
  // and penwidth overrides without touching the labels.
  fmt::print(out,
             R"({} -> {} [taillabel="{}" headlabel="{}"])"
             "\n",
             conn.SrcId(), conn.DstId(), SignChar(conn.SrcSign()), SignChar(conn.DstSign()));
}

void RenderOverlayEdgeStatement(std::ostream& out, NodeID src, NodeID dst,
                                OverlayEdge const& overlay) {
  if (overlay.mColors.empty()) return;
  auto const color_list = absl::StrJoin(overlay.mColors, ":");
  fmt::print(out, R"({} -> {} [color="{}")", src, dst, color_list);
  if (overlay.mPenWidth != 0) fmt::print(out, " penwidth={}", overlay.mPenWidth);
  if (overlay.mDashed) out << R"( style="dashed")";
  out << "]\n";
}

void RenderToOstream(std::ostream& out,
                     absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                     DotPlan const& plan) {
  out << DOT_PREAMBLE;
  fmt::print(out, "subgraph {} {{\n", plan.mSubgraphName);

  auto const node_overlays = BuildNodeOverlays(plan.mNodeLayers);
  auto const edge_overlays = BuildEdgeOverlays(plan.mEdgeLayers);

  // First pass: nodes + base edges.
  for (auto const& [node_id, node_ptr] : graph) {
    if (node_ptr->GetComponentId() != plan.mCompId) continue;
    RenderNodeDot(out, *node_ptr, node_id, node_overlays);
    for (Edge const& conn : *node_ptr) {
      RenderBaseEdgeDot(out, conn);
    }
  }

  // Second pass: overlay edge statements, expanded to BOTH DOT directions
  // so anti-parallel mirror splines under `neato` carry the same style.
  for (auto const& [key, overlay] : edge_overlays) {
    RenderOverlayEdgeStatement(out, key.mLo, key.mHi, overlay);
    if (key.mLo != key.mHi) {
      RenderOverlayEdgeStatement(out, key.mHi, key.mLo, overlay);
    }
  }

  out << DOT_FOOTER;
}

}  // namespace

auto SerializeToDotString(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                          DotPlan const& plan) -> std::string {
  std::ostringstream stream;
  RenderToOstream(stream, graph, plan);
  return std::move(stream).str();
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

}  // namespace lancet::cbdg
