#ifndef SRC_LANCET_CBDG_SERIALIZE_DOT_H_
#define SRC_LANCET_CBDG_SERIALIZE_DOT_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

namespace lancet::cbdg {

// ============================================================================
// Overlay sets for debug highlighting
// ============================================================================
// Bundles node and edge sets used to color-code DOT graph nodes/edges.
// Highlight sets are shown in accent colors; background sets are dimmed.
struct DotOverlaySets {
  // ── 8B Align ────────────────────────────────────────────────────────────
  absl::flat_hash_set<NodeID> mNodes;
  absl::flat_hash_set<Edge> mEdges;
};

/// Serialize a graph component to a Graphviz DOT file for debug visualization.
/// If `walks` is non-empty, appends per-walk edge color override stanzas after
/// the normal graph rendering. Each walk's edges are wrapped in an anonymous
/// subgraph block with `edge [color=X penwidth=2]` defaults. Walks are iterated
/// in reverse order so DOT's last-writer-wins gives highest-confidence walk
/// priority on shared edges.
void SerializeToDot(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph,
                    std::filesystem::path const& out_path, usize comp_id = 0,
                    DotOverlaySets const& highlight = {}, DotOverlaySets const& background = {},
                    absl::Span<Path const> walks = {});

/// Render a graph component to a DOT-formatted string for in-memory buffering.
/// Same content as SerializeToDot but returns the body instead of writing to
/// disk; the caller (typically DotSnapshotBuffer) handles the file I/O on
/// commit. `subgraph_name` becomes the `subgraph` identifier — usually the
/// stem of the eventual output filename.
[[nodiscard]] auto SerializeToDotString(
    absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph, std::string_view subgraph_name,
    usize comp_id = 0, DotOverlaySets const& highlight = {}, DotOverlaySets const& background = {},
    absl::Span<Path const> walks = {}) -> std::string;

/// Pipeline stage tag for DOT graph debug snapshots.
enum class GraphState : u8 {
  FIRST_LOW_COV_REMOVAL = 0,
  FOUND_REF_ANCHORS = 1,
  FIRST_COMPRESSION = 2,
  SECOND_LOW_COV_REMOVAL = 3,
  SECOND_COMPRESSION = 4,
  SHORT_TIP_REMOVAL = 5,
  FULLY_PRUNED_GRAPH = 6,
  ENUMERATED_WALKS = 7
};

/// Return a maximally-distinct hex color (#RRGGBB) for walk `walk_index`.
/// 64-entry palette pre-computed via k-means in CIE L*a*b* with farthest-first
/// ordering. Indices beyond 64 cycle via modulo.
[[nodiscard]] auto WalkColor(usize walk_index) -> std::string_view;

#ifdef LANCET_DEVELOP_MODE
[[nodiscard]] auto ToString(GraphState state) -> std::string;
#endif

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_SERIALIZE_DOT_H_
