#include "lancet/cbdg/serialize_dot.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "spdlog/fmt/bundled/ostream.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

namespace lancet::cbdg {

#ifdef LANCET_DEVELOP_MODE
auto ToString(GraphState const state) -> std::string {
  switch (state) {
    case GraphState::FIRST_LOW_COV_REMOVAL:
      return "low_cov_removal1";
    case GraphState::FOUND_REF_ANCHORS:
      return "found_ref_anchors";
    case GraphState::FIRST_COMPRESSION:
      return "compression1";
    case GraphState::SECOND_LOW_COV_REMOVAL:
      return "low_cov_removal2";
    case GraphState::SECOND_COMPRESSION:
      return "compression2";
    case GraphState::SHORT_TIP_REMOVAL:
      return "short_tip_removal";
    default:
      break;
  }

  return "fully_pruned";
}
#endif

namespace {

using namespace std::string_view_literals;

// ── DOT Graph Preamble ──────────────────────────────────────────────────
constexpr auto DOT_PREAMBLE = R"raw(strict digraph G {
graph [layout=neato,bgcolor=black,size="120,180",ratio=compress,rankdir=LR,overlap=vpsc,overlap_shrink=true,start=self];
node [style=filled,fontsize=2,width=2,height=2,fixedsize=false];
edge [color=gray,fontsize=8,fontcolor=floralwhite,len=3,fixedsize=false,headclip=true,tailclip=true];
)raw"sv;

constexpr auto DOT_FOOTER = "}\n}\n"sv;

// ============================================================================
// Node Colors
// ============================================================================
constexpr auto COLOR_BACKGROUND = "darkgray"sv;
constexpr auto COLOR_HIGHLIGHT = "orchid"sv;
constexpr auto COLOR_SHARED = "steelblue"sv;
constexpr auto COLOR_CASE_ONLY = "indianred"sv;
constexpr auto COLOR_CTRL_ONLY = "mediumseagreen"sv;
constexpr auto COLOR_DEFAULT = "lightblue"sv;

// ============================================================================
// Edge Styles
// ============================================================================
constexpr auto STYLE_DOTTED = "dotted"sv;
constexpr auto STYLE_SOLID = "solid"sv;
constexpr auto ATTR_HIGHLIGHT_COLOR = R"raw( color="goldenrod")raw"sv;
constexpr auto ATTR_NO_COLOR = ""sv;

// ── Format Templates ────────────────────────────────────────────────────
constexpr auto EDGE_FMT = R"raw({} -> {} [taillabel="{}" headlabel="{}" style="{}"{}]
)raw"sv;

constexpr auto NODE_FMT =
    R"raw({} [shape=circle fillcolor={} label="{}\n{}\n {}:{}\nlength={}\ncoverage={}"]
)raw"sv;

auto NodeFillColor(Node const& node, NodeID node_id, DotOverlaySets const& highlight,
                   DotOverlaySets const& background) -> std::string_view {
  if (background.mNodes.contains(node_id)) return COLOR_BACKGROUND;
  if (highlight.mNodes.contains(node_id)) return COLOR_HIGHLIGHT;
  if (node.IsShared()) return COLOR_SHARED;
  if (node.IsCaseOnly()) return COLOR_CASE_ONLY;
  if (node.IsCtrlOnly()) return COLOR_CTRL_ONLY;
  return COLOR_DEFAULT;
}

auto SignChar(Kmer::Sign sign) -> char {
  return sign == Kmer::Sign::PLUS ? '+' : '-';
}

void WriteEdgeDot(std::ofstream& out, Edge const& conn, DotOverlaySets const& highlight,
                  DotOverlaySets const& background) {
  fmt::print(out, EDGE_FMT, conn.SrcId(), conn.DstId(), SignChar(conn.SrcSign()),
             SignChar(conn.DstSign()),
             background.mEdges.contains(conn) ? STYLE_DOTTED : STYLE_SOLID,
             highlight.mEdges.contains(conn) ? ATTR_HIGHLIGHT_COLOR : ATTR_NO_COLOR);
}

void WriteNodeDot(std::ofstream& out, Node const& node, NodeID node_id,
                  DotOverlaySets const& highlight, DotOverlaySets const& background) {
  auto const dflt_seq = node.SequenceFor(Kmer::Ordering::DEFAULT);
  auto const oppo_seq = node.SequenceFor(Kmer::Ordering::OPPOSITE);
  auto const rev_oppo_seq = std::string(oppo_seq.crbegin(), oppo_seq.crend());
  auto const fill_color = NodeFillColor(node, node_id, highlight, background);

  fmt::print(out, NODE_FMT, node_id, fill_color, dflt_seq, rev_oppo_seq, node_id,
             SignChar(node.SignFor(Kmer::Ordering::DEFAULT)), node.Length(),
             node.TotalReadSupport());

  for (Edge const& conn : node) {
    WriteEdgeDot(out, conn, highlight, background);
  }
}

}  // namespace

void SerializeToDot(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                    std::filesystem::path const& out_path, usize const comp_id,
                    DotOverlaySets const& highlight, DotOverlaySets const& background) {
  std::ofstream out_handle(out_path, std::ios::trunc);
  out_handle << DOT_PREAMBLE;
  fmt::print(out_handle, "subgraph {} {{\n", out_path.stem().string());

  for (auto const& [node_id, node_ptr] : graph) {
    if (node_ptr->GetComponentId() != comp_id) continue;
    WriteNodeDot(out_handle, *node_ptr, node_id, highlight, background);
  }

  out_handle << DOT_FOOTER;
  out_handle.close();
}

}  // namespace lancet::cbdg
