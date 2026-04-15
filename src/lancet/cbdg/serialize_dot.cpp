#include "lancet/cbdg/serialize_dot.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ostream.h"

#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <string_view>

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)
void SerializeToDot(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                    std::filesystem::path const& out_path, usize const comp_id,
                    absl::flat_hash_set<NodeID> const& nodes_highlight,
                    absl::flat_hash_set<Edge> const& edges_highlight,
                    absl::flat_hash_set<NodeID> const& nodes_background,
                    absl::flat_hash_set<Edge> const& edges_background) {
  std::ofstream out_handle(out_path, std::ios::trunc);
  using namespace std::string_view_literals;

  out_handle << R"raw(strict digraph G {
graph [layout=neato,bgcolor=black,size="120,180",ratio=compress,rankdir=LR,overlap=vpsc,overlap_shrink=true,start=self];
node [style=filled,fontsize=2,width=2,height=2,fixedsize=false];
edge [color=gray,fontsize=8,fontcolor=floralwhite,len=3,fixedsize=false,headclip=true,tailclip=true];
)raw"sv;

  fmt::print(out_handle, "subgraph {} {{\n", out_path.stem().string());

  for (auto const& item : graph) {
    if (item.second->GetComponentId() != comp_id) continue;

    auto const dflt_seq = item.second->SequenceFor(Kmer::Ordering::DEFAULT);
    auto const oppo_seq = item.second->SequenceFor(Kmer::Ordering::OPPOSITE);
    auto const rev_oppo_seq = std::string(oppo_seq.crbegin(), oppo_seq.crend());
    auto const sign_dflt =
        item.second->SignFor(Kmer::Ordering::DEFAULT) == Kmer::Sign::PLUS ? '+' : '-';
    auto const is_background_node = nodes_background.contains(item.first);
    // NOLINTBEGIN(readability-avoid-nested-conditional-operator)
    auto const fill_color = is_background_node                     ? "darkgray"sv
                            : nodes_highlight.contains(item.first) ? "orchid"sv
                            : item.second->IsShared()              ? "steelblue"sv
                            : item.second->IsCaseOnly()            ? "indianred"sv
                            : item.second->IsCtrlOnly()            ? "mediumseagreen"sv
                                                                   : "lightblue"sv;
    // NOLINTEND(readability-avoid-nested-conditional-operator)

    fmt::print(out_handle,
               R"raw({} [shape=circle fillcolor={} label="{}\n{}\n {}:{}\nlength={}\ncoverage={}"]
)raw",
               item.first, fill_color, dflt_seq, rev_oppo_seq, item.first, sign_dflt,
               item.second->Length(), item.second->TotalReadSupport());

    for (Edge const& conn : *item.second) {
      auto const src_sign = conn.SrcSign() == Kmer::Sign::PLUS ? '+' : '-';
      auto const dst_sign = conn.DstSign() == Kmer::Sign::PLUS ? '+' : '-';
      auto const is_background_edge = edges_background.contains(conn);
      auto const is_highlight_edge = edges_highlight.contains(conn);
      fmt::print(out_handle, R"raw({} -> {} [taillabel="{}" headlabel="{}" style="{}"{}]
)raw",
                 conn.SrcId(), conn.DstId(), src_sign, dst_sign,
                 is_background_edge ? "dotted"sv : "solid"sv,
                 is_highlight_edge ? R"raw( color="goldenrod")raw"sv : ""sv);
    }
  }

  out_handle << "}\n}\n"sv;
  out_handle.close();
}

}  // namespace lancet::cbdg
