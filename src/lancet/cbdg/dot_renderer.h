#ifndef SRC_LANCET_CBDG_DOT_RENDERER_H_
#define SRC_LANCET_CBDG_DOT_RENDERER_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/dot_plan.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"

#include <memory>
#include <string>

namespace lancet::cbdg {

/// Render one component of the graph to a DOT string per the given plan.
/// Pure function over its inputs — no I/O, no globals. The caller
/// (typically `DotSnapshotBuffer` via `Graph::BufferFinalSnapshot` /
/// `Graph::BufferStageSnapshot`) handles the file write at commit time.
///
/// Layer composition rules:
///   • Node layers in ascending `mZOrder` override prior layers' fields per
///     attribute slot (fillcolor, border, peripheries, stripe accent). Role
///     fillcolor is a default that any layer can override.
///   • Edge layers' colors are concatenated into a graphviz `colorList`
///     (e.g. `color="#AAA:#BBB:#CCC"`) so multi-walk shared edges render as
///     parallel-color stripes with no information loss.
///   • Every overlayed `LogicalEdge` is emitted in BOTH directional DOT
///     statements (`lo->hi` and `hi->lo`) with the same style — bidirected
///     mirror coverage is a renderer invariant, not a per-call concern.
[[nodiscard]] auto SerializeToDotString(
    absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& graph, DotPlan const& plan)
    -> std::string;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_RENDERER_H_
