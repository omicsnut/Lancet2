#include "lancet/cbdg/traversal_index.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"

#include <memory>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
//  BuildTraversalIndex — Construct Flat CSR Adjacency List
// ============================================================================
//
// Converts the hash-map-based NodeTable into a contiguous, integer-indexed
// adjacency list for a single connected component. This is built ONCE after
// all graph mutations (pruning, compression, tip removal) are complete.
//
// CONSTRUCTION PHASES
// --------------------
//  Phase 1: Assign contiguous u32 indices to nodes in this component.
//  Phase 2: Count outgoing edges per state (for CSR range sizing).
//  Phase 3: Compute prefix-sum offsets for each state's edge range.
//  Phase 4: Fill the adjacency list and assign edge ordinals.
//  Phase 5: Set source and sink state indices.
//
// COST: O(V + E) time and memory. The nid_to_flat hash map is only used
// during construction; all subsequent operations are flat-array-only.
//
auto BuildTraversalIndex(absl::flat_hash_map<NodeID, std::unique_ptr<Node>> const& nodes,
                         NodeIDPair const& source_and_sink_ids, usize const component_id)
    -> TraversalIndex {
  TraversalIndex traversal_index;

  // Phase 1: Assign contiguous u32 indices to nodes in this component
  absl::flat_hash_map<NodeID, u32> nid_to_flat;
  nid_to_flat.reserve(nodes.size());

  for (auto const& [node_id, node_ptr] : nodes) {
    if (node_ptr->GetComponentId() != component_id) continue;

    auto const flat = static_cast<u32>(traversal_index.mNodes.size());
    traversal_index.mNodes.push_back(node_ptr.get());
    traversal_index.mNodeIds.push_back(node_id);
    nid_to_flat.emplace(node_id, flat);
  }

  u32 const num_nodes = traversal_index.NumNodes();
  u32 const num_states = num_nodes * 2;
  traversal_index.mAdjRanges.resize(num_states, {.mStart = 0, .mCount = 0});

  // Phase 2: Count outgoing edges per state (for range sizing)
  for (u32 node_idx = 0; node_idx < num_nodes; node_idx++) {
    Node const* node = traversal_index.mNodes[node_idx];
    for (Edge const& edge : *node) {
      // Only count edges whose destination is in this component
      if (!nid_to_flat.contains(edge.DstId())) continue;

      u32 const state = TraversalIndex::MakeState(node_idx, edge.SrcSign());
      traversal_index.mAdjRanges[state].mCount++;
    }
  }

  // Phase 3: Compute starting offsets via prefix sum
  u32 offset = 0;
  for (u32 state_idx = 0; state_idx < num_states; state_idx++) {
    traversal_index.mAdjRanges[state_idx].mStart = offset;
    offset += traversal_index.mAdjRanges[state_idx].mCount;
    traversal_index.mAdjRanges[state_idx].mCount = 0;  // reset count; re-filled in phase 4
  }
  traversal_index.mAdjList.resize(offset);

  // Phase 4: Fill adjacency list entries and assign edge ordinals
  absl::flat_hash_map<Edge, u32> edge_to_ordinal;
  edge_to_ordinal.reserve(offset);

  for (u32 node_idx = 0; node_idx < num_nodes; node_idx++) {
    Node const* node = traversal_index.mNodes[node_idx];
    for (Edge const& edge : *node) {
      auto const dst_it = nid_to_flat.find(edge.DstId());
      if (dst_it == nid_to_flat.end()) continue;

      u32 const src_state = TraversalIndex::MakeState(node_idx, edge.SrcSign());
      u32 const dst_state = TraversalIndex::MakeState(dst_it->second, edge.DstSign());

      // Assign or reuse edge ordinal (edges appear at both endpoints as forward + mirror)
      u32 ordinal = 0;
      auto [iter, inserted] =
          edge_to_ordinal.emplace(edge, static_cast<u32>(traversal_index.mOrigEdges.size()));

      if (inserted) {
        ordinal = iter->second;
        traversal_index.mOrigEdges.push_back(edge);
      } else {
        ordinal = iter->second;
      }

      auto& range = traversal_index.mAdjRanges[src_state];
      traversal_index.mAdjList[range.mStart + range.mCount] = {.mDstState = dst_state,
                                                               .mEdgeOrdinal = ordinal};

      range.mCount++;
    }
  }

  // Phase 5: Set source and sink states
  auto const [source_id, sink_id] = source_and_sink_ids;
  auto const src_flat = nid_to_flat.at(source_id);
  auto const snk_flat = nid_to_flat.at(sink_id);
  auto const src_sign = traversal_index.mNodes[src_flat]->SignFor(Kmer::Ordering::DEFAULT);
  traversal_index.mSrcState = TraversalIndex::MakeState(src_flat, src_sign);
  traversal_index.mSnkNodeIdx = snk_flat;

  return traversal_index;
}

}  // namespace lancet::cbdg
