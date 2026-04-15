#ifndef SRC_LANCET_CBDG_CYCLE_FINDER_H_
#define SRC_LANCET_CBDG_CYCLE_FINDER_H_

#include "lancet/cbdg/traversal_index.h"

namespace lancet::cbdg {

/// O(V+E) cycle detection using three-color DFS on the flat adjacency list.
/// Tracks (NodeIdx, Sign) state pairs — a node reached via '+' and via '-' are
/// different states, correctly handling bidirected sign continuity (BCALM2).
/// Returns true if any cycle is reachable from the source node.
[[nodiscard]] auto HasCycle(TraversalIndex const& idx) -> bool;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_CYCLE_FINDER_H_
