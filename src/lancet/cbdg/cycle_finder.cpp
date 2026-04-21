#include "lancet/cbdg/cycle_finder.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/traversal_index.h"

#include <vector>

namespace lancet::cbdg {

// ============================================================================
//  HasCycle — O(V+E) Three-Color DFS on the Flat Adjacency List
// ============================================================================
//
// ALGORITHM
// ============================================================================
// Standard directed-graph cycle detection using three colors:
//   WHITE (0) = unvisited
//    GRAY (1) = on the current DFS stack (ancestor in the current path)
//   BLACK (2) = fully explored (all descendants visited)
//
// A cycle exists iff DFS encounters a "back edge" — an edge leading to a
// GRAY state. GRAY states are ancestors on the current DFS path, so a back
// edge forms a cycle.
//
// WHY THIS REPLACES THE OLD APPROACH
// ============================================================================
// The old HasCycle used backtracking (erase from visited set on return),
// which explored ALL paths from source — exponential in high-branching
// graphs. Three-color DFS visits each state exactly once: O(V+E).
// Profile data showed ~51.6s in the old HasCycle; this should be <1ms.
//
// BIDIRECTED SIGN CONTINUITY
// ============================================================================
// In the BCALM2 bidirected model, walks must satisfy sign continuity:
//   edge.DstSign must match the next edge.SrcSign
//
// We track state as (node_flat_idx, sign), so a node visited via '+' and
// via '-' are different DFS states. The TraversalIndex adjacency list is
// already partitioned by (node, sign), so edge iteration naturally respects
// sign continuity.
//
//   Example: DFS from source(+)
//   ┌──────────────────────────────────────────────────────┐
//   │  Visit state (A,+) → GRAY                            │
//   │    Edge to (B,+) → WHITE → push (B,+)                │
//   │      Edge to (C,-) → WHITE → push (C,-)              │
//   │        Edge to (A,+) → GRAY! → CYCLE FOUND           │
//   │      Edge to (A,-) → WHITE → push (A,-)              │ ← same node, different sign
//   │        ...no back edge → BLACK                       │
//   └──────────────────────────────────────────────────────┘
//
// FUTURE: Tarjan's SCC could provide cycle sizes and weakest edges for
// smarter retry-vs-skip decisions. For now, just detect presence/absence.
//
auto HasCycle(TraversalIndex const& idx) -> bool {
  enum class Color : u8 { WHITE = 0, GRAY = 1, BLACK = 2 };

  // Flat color array indexed by state_idx = node_flat_idx * 2 + sign_offset.
  // O(1) access per state — no hash lookups, cache-line-friendly.
  std::vector<Color> color(idx.NumStates(), Color::WHITE);

  // Iterative DFS stack frame: current state + position in its outgoing edge list.
  // Using an explicit stack avoids recursion depth limits on large graphs.
  struct DfsFrame {
    u32 mStateIdx;
    u32 mEdgePos;  // next edge to explore within this state's adjacency range
  };

  std::vector<DfsFrame> stack;
  stack.reserve(idx.NumNodes());

  color[idx.mSrcState] = Color::GRAY;
  stack.push_back({.mStateIdx = idx.mSrcState, .mEdgePos = 0});

  while (!stack.empty()) {
    auto& frame = stack.back();
    auto const& range = idx.mAdjRanges[frame.mStateIdx];

    // All children of this state explored → mark BLACK and backtrack
    if (frame.mEdgePos >= range.mCount) {
      color[frame.mStateIdx] = Color::BLACK;
      stack.pop_back();
      continue;
    }

    // Examine the next outgoing edge from this state
    auto const& out = idx.mAdjList[range.mStart + frame.mEdgePos];
    frame.mEdgePos++;

    if (color[out.mDstState] == Color::GRAY) return true;  // Back edge → cycle!

    if (color[out.mDstState] != Color::WHITE) continue;  // BLACK → already finished, skip

    // WHITE → unvisited. Push child onto DFS stack.
    color[out.mDstState] = Color::GRAY;
    stack.push_back({.mStateIdx = out.mDstState, .mEdgePos = 0});
  }

  return false;
}

}  // namespace lancet::cbdg
