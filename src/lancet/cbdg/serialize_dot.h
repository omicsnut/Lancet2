#ifndef SRC_LANCET_CBDG_SERIALIZE_DOT_H_
#define SRC_LANCET_CBDG_SERIALIZE_DOT_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include <filesystem>
#include <memory>
#include <string>

namespace lancet::cbdg {

/// Serialize a graph component to a Graphviz DOT file for debug visualization.
void SerializeToDot(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                    std::filesystem::path const& out_path, usize comp_id = 0,
                    absl::flat_hash_set<NodeID> const& nodes_highlight = {},
                    absl::flat_hash_set<Edge> const& edges_highlight = {},
                    absl::flat_hash_set<NodeID> const& nodes_background = {},
                    absl::flat_hash_set<Edge> const& edges_background = {});

/// Pipeline stage tag for DOT graph debug snapshots.
enum class GraphState : u8 {
  FIRST_LOW_COV_REMOVAL = 0,
  FOUND_REF_ANCHORS = 1,
  FIRST_COMPRESSION = 2,
  SECOND_LOW_COV_REMOVAL = 3,
  SECOND_COMPRESSION = 4,
  SHORT_TIP_REMOVAL = 5,
  SHORT_LINK_REMOVAL = 6,
  FULLY_PRUNED_GRAPH = 7
};

#ifdef LANCET_DEVELOP_MODE
[[nodiscard]] auto ToString(GraphState state) -> std::string;
#endif

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_SERIALIZE_DOT_H_
