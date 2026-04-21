#include "lancet/cbdg/graph.h"

#include "lancet/base/assert.h"
#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/cycle_finder.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/max_flow.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/traversal_index.h"
#include "lancet/hts/phred_quality.h"

#include "absl/container/chunked_queue.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>

namespace lancet::cbdg {

/// Pipeline architecture for haplotype assembly from a colored de Bruijn graph:
///
///  ┌─────────────┐
///  │ Outer loop: │  Iterate k from min_k to max_k in steps of mKmerStepLen.
///  │ k-value scan│  If haplotypes are found at any k, stop.
///  │             │  If a cycle is detected or graph is too complex, abandon
///  │             │  this k and continue to next (via should_retry_kmer flag).
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐
///  │  BuildGraph │  Build the bidirected de Bruijn graph from reads + reference.
///  │  + prune    │  Remove low-coverage nodes, mark connected components.
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐   For each connected component with valid source/sink anchors:
///  │ Inner loop: │
///  │per-component│   1. Compress linear chains, remove low-cov nodes, remove tips
///  │             │   2. Build TraversalIndex (flat adjacency list) on frozen graph
///  │             │   3. HasCycle via O(V+E) three-color DFS on flat arrays
///  │             │   4. Log graph complexity metrics (cyclomatic, density, etc.)
///  │             │   5. Enumerate all source→sink walks via BFS walk-tree
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐
///  │  Assemble   │  Deduplicate haplotype sequences, prepend reference anchor.
///  │  + return   │  Return assembled haplotypes for genotyping / variant calling.
///  └─────────────┘
///
/// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
auto Graph::BuildComponentHaplotypes(RegionPtr region, ReadList reads) -> Result {
  mReads = reads;
  mRegion = std::move(region);

  lancet::base::Timer timer;
  GraphHaps per_comp_haps;
  std::string_view ref_anchor_seq;
  std::vector<usize> anchor_start_idxs;
  std::vector<GraphComplexity> component_metrics;
  absl::flat_hash_set<MateMer> mate_mers;

  static constexpr usize DEFAULT_EST_NUM_NODES = 32'768;
  static constexpr usize DEFAULT_MIN_ANCHOR_LENGTH = 150;
  static constexpr f64 DEFAULT_PCT_NODES_NEEDED = 10.0;

  // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
  auto const reg_str = mRegion->ToSamtoolsRegion();
  mCurrK = mParams.mMinKmerLen - mParams.mKmerStepLen;

  // Outer loop: increment k and retry until haplotypes are found or k is exhausted.
  // Replaces the old `goto IncrementKmerAndRetry` with structured control flow.
  while (per_comp_haps.empty() && (mCurrK + mParams.mKmerStepLen) <= mParams.mMaxKmerLen) {
    mCurrK += mParams.mKmerStepLen;
    timer.Reset();
    mSourceAndSinkIds = {0, 0};
    mNodes.reserve(DEFAULT_EST_NUM_NODES);

    // Skip this k if the reference itself has a repeated k-mer — the de Bruijn
    // graph would contain a cycle by construction, making assembly pointless.
    if (HasExactOrApproxRepeat(mRegion->SeqView(), mCurrK)) continue;

    mNodes.clear();
    BuildGraph(mate_mers);
    LOG_TRACE("Done building graph for {} with k={}, nodes={}, reads={}", reg_str, mCurrK,
              mNodes.size(), mReads.size())

    RemoveLowCovNodes(0);
    WriteDotDevelop(GraphState::FIRST_LOW_COV_REMOVAL, 0);

    auto const components = MarkConnectedComponents();
    per_comp_haps.reserve(components.size());
    anchor_start_idxs.reserve(components.size());
    LOG_TRACE("Found {} connected components in graph for {} with k={}", components.size(), reg_str,
              mCurrK)

    // Inner loop: process each connected component with valid source/sink anchors.
    // The should_retry_kmer flag is set to true by cycle detection to abandon all
    // remaining components at this k and retry at a higher k value.
    bool should_retry_kmer = false;
    for (auto const& cinfo : components) {
      if (should_retry_kmer) break;
      if (cinfo.mPctNodes < DEFAULT_PCT_NODES_NEEDED) continue;

      auto const comp_id = cinfo.mCompId;
      auto const source = FindSource(comp_id);
      auto const sink = FindSink(comp_id);

      if (!source.mFoundAnchor || !sink.mFoundAnchor || source.mAnchorId == sink.mAnchorId) {
        LOG_TRACE("Skipping comp{} in graph for {} because source/sink was not found", comp_id,
                  reg_str)
        continue;
      }

      auto const current_anchor_length = RefAnchorLength(source, sink, mCurrK);
      if (current_anchor_length < DEFAULT_MIN_ANCHOR_LENGTH) continue;

      LOG_TRACE("Found {}bp ref anchor for {} comp={} with k={}", current_anchor_length, reg_str,
                comp_id, mCurrK)

      mSourceAndSinkIds = NodeIDPair{source.mAnchorId, sink.mAnchorId};
      ref_anchor_seq = mRegion->SeqView().substr(source.mRefOffset, current_anchor_length);
      WriteDotDevelop(GraphState::FOUND_REF_ANCHORS, comp_id);

      PruneComponent(comp_id);

      // Build the flat traversal index on the frozen (fully-pruned) graph.
      // This maps NodeID -> contiguous u32 and constructs the CSR adjacency list.
      // Both HasCycle and MaxFlow operate on this flat structure for O(1) state tracking.
      auto const trav_idx = BuildTraversalIndex(mNodes, mSourceAndSinkIds, comp_id);

      // O(V+E) cycle detection using three-color DFS on the flat adjacency list.
      // See HasCycle() implementation for bidirected sign-continuity handling.
      if (HasCycle(trav_idx)) {
        LOG_TRACE("Cycle found in pruned graph for {} comp={} with k={}", reg_str, comp_id, mCurrK)
        should_retry_kmer = true;
        break;
      }

      // Log graph complexity metrics for debugging / correlating with runtime.
      // All metrics are O(V+E) to compute and help identify pathological windows.
      // Skip walk enumeration on pathological graphs — retry with larger k to
      // collapse branches. Same control flow as the HasCycle guard above.
      auto const cplx = ComputeComponentComplexity(comp_id);
      if (cplx.IsComplex()) {
        LOG_DEBUG("Graph too complex for {} comp={} k={}: CC={} branches={}", reg_str, comp_id,
                  mCurrK, cplx.CyclomaticComplexity(), cplx.NumBranchPoints())
        should_retry_kmer = true;
        break;
      }

      WriteDot(GraphState::FULLY_PRUNED_GRAPH, comp_id);

      auto haplotypes = EnumerateAndSortHaplotypes(comp_id, trav_idx, ref_anchor_seq);
      if (!haplotypes.empty()) {
        per_comp_haps.emplace_back(std::move(haplotypes));
        anchor_start_idxs.emplace_back(source.mRefOffset);
        component_metrics.emplace_back(cplx);
      }
    }

    // If any component triggered a retry, discard partial results and try next k
    if (should_retry_kmer) {
      per_comp_haps.clear();
      anchor_start_idxs.clear();
      component_metrics.clear();
      continue;
    }
  }

  // Count ALT haplotypes per component (excluding the leading reference path at index 0).
  // NOLINTNEXTLINE(clang-analyzer-deadcode.DeadStores)
  auto const num_haps = std::accumulate(
      per_comp_haps.cbegin(), per_comp_haps.cend(), u64{0},
      [](u64 const sum, auto const& comp_haps) -> u64 { return sum + comp_haps.size() - 1; });

  // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
  auto const human_rt = timer.HumanRuntime();

  LOG_TRACE("Assembled {} haplotypes for {} with k={} in {}", num_haps, reg_str, mCurrK, human_rt)

  return {.mGraphHaplotypes = per_comp_haps,
          .mAnchorStartIdxs = anchor_start_idxs,
          .mComponentMetrics = std::move(component_metrics)};
}

// ============================================================================
// Phase 1: Graph Construction
// ============================================================================

void Graph::BuildGraph(absl::flat_hash_set<MateMer>& mate_mers) {
  mRefNodeIds.clear();
  auto const ref_nodes = AddNodes(mRegion->SeqView(), Label(Label::REFERENCE));
  mRefNodeIds.reserve(ref_nodes.size());
  std::ranges::transform(ref_nodes, std::back_inserter(mRefNodeIds),
                         [](Node const* node) -> NodeID { return node->Identifier(); });

  mate_mers.clear();
  std::vector<f64> per_base_error_probs;
  static constexpr usize EXPECTED_READ_LENGTH = 150;
  per_base_error_probs.reserve(mReads.empty() ? EXPECTED_READ_LENGTH : mReads[0].Length());

  for (auto const& read : mReads) {
    if (!read.PassesAlnFilters()) continue;

    // Amortize Phred → error-prob look up to once per read instead of k times per kmer.
    // Prefix sum enables O(1) expected-error range queries:
    // expected_errors(i, i+k) = prefix[i+k] − prefix[i].
    per_base_error_probs.clear();
    std::ranges::transform(read.QualView(), std::back_inserter(per_base_error_probs),
                           [](u8 const qual) -> f64 { return hts::PhredToErrorProb(qual); });
    std::vector<f64> prefix_sum(per_base_error_probs.size() + 1, 0.0);
    std::partial_sum(per_base_error_probs.cbegin(), per_base_error_probs.cend(),
                     prefix_sum.begin() + 1);

    usize offset = 0;
    auto const added_nodes = AddNodes(read.SeqView(), read.SrcLabel());

    for (auto* node : added_nodes) {
      MateMer mm_pair{
          .mQname = read.QnameView(), .mKmerHash = node->Identifier(), .mTagKind = read.TagKind()};

      // Filter out low-quality kmers by expected error count (floor of summed Phred error
      // probabilities). Kmers with ≥1 expected error get no read support,
      // ensuring they are removed during the subsequent low-coverage pruning pass.
      // See https://www.drive5.com/usearch/manual/exp_errs.html
      // See https://doi.org/10.1093/bioinformatics/btv401 for proof on expected errors
      // O(1) expected-error check for kmer at [offset, offset+k)
      auto const raw_expected_error = prefix_sum[offset + mCurrK] - prefix_sum[offset];
      auto const expected_error = static_cast<i64>(std::floor(raw_expected_error));
      offset++;

      if (expected_error > 0 || mate_mers.contains(mm_pair)) continue;
      node->IncrementReadSupport(read.SampleIndex(), read.TagKind());
      mate_mers.emplace(mm_pair);
    }
  }
}

auto Graph::AddNodes(std::string_view sequence, Label const label) -> std::vector<Node*> {
  std::vector<Node*> result;
  auto const kplus_ones = lancet::base::SlidingView(sequence, mCurrK + 1);
  result.reserve(kplus_ones.size() + 1);

  for (usize mer_idx = 0; mer_idx < kplus_ones.size(); ++mer_idx) {
    auto const seq1 = absl::ClippedSubstr(kplus_ones[mer_idx], 0, mCurrK);
    auto const seq2 = absl::ClippedSubstr(kplus_ones[mer_idx], 1, mCurrK);

    auto left_mer = Kmer(seq1);
    auto right_mer = Kmer(seq2);
    auto const left_id = left_mer.Identifier();
    auto const right_id = right_mer.Identifier();

    mNodes.try_emplace(left_id, std::make_unique<Node>(std::move(left_mer), label));
    mNodes.try_emplace(right_id, std::make_unique<Node>(std::move(right_mer), label));

    auto& first = mNodes.at(left_id);
    auto& second = mNodes.at(right_id);

    if (mer_idx == 0) result.emplace_back(first.get());

    static constexpr auto DEFAULT_ORDER = Kmer::Ordering::DEFAULT;
    auto const fwd_edge =
        MakeFwdEdgeKind({first->SignFor(DEFAULT_ORDER), second->SignFor(DEFAULT_ORDER)});
    first->EmplaceEdge(NodeIDPair{left_id, right_id}, fwd_edge);
    second->EmplaceEdge(NodeIDPair{right_id, left_id}, RevEdgeKind(fwd_edge));

    result.emplace_back(second.get());
  }

  return result;
}

// ============================================================================
// Phase 2: Node Removal + Connected Components
// ============================================================================

void Graph::RemoveNode(NodeTable::iterator itr) {
  if (itr == mNodes.end()) return;

  // remove all incoming edges to the node first
  for (Edge const& conn : *itr->second) {
    if (conn.IsSelfLoop()) continue;

    auto nbour_itr = mNodes.find(conn.DstId());
    if (nbour_itr != mNodes.end()) nbour_itr->second->EraseEdge(conn.MirrorEdge());
  }

  mNodes.erase(itr);
}

void Graph::RemoveLowCovNodes(usize const component_id) {
  std::vector<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());

  auto const [source_id, sink_id] = mSourceAndSinkIds;
  for (auto const& [nid, node_ptr] : mNodes) {
    if (node_ptr->GetComponentId() != component_id) continue;
    if (nid == source_id || nid == sink_id) continue;

    // Remove nodes where every sample has at most 1 read (unreliable k-mers)
    // or total coverage is below the user-configured minimum.
    auto const all_singletons = node_ptr->IsAllSingletons();
    auto const total_sample_cov = node_ptr->TotalReadSupport();

    if (all_singletons || total_sample_cov < mParams.mMinNodeCov) {
      remove_nids.emplace_back(nid);
    }
  }

  if (!remove_nids.empty()) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    auto const region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Removing {:.4f}% (or) {} low cov nodes for {} in comp{} with k={}",
              100.0 * (static_cast<f64>(remove_nids.size()) / static_cast<f64>(mNodes.size())),
              remove_nids.size(), region_str, component_id, mCurrK)

    RemoveNodes(absl::MakeConstSpan(remove_nids));
  }
}

auto Graph::MarkConnectedComponents() -> std::vector<ComponentInfo> {
  usize current_component = 0;
  std::vector<ComponentInfo> results_info;

#ifdef LANCET_DEVELOP_MODE
  static auto const IS_UNASSIGNED = [](NodeTable::const_reference item) {
    return item.second->GetComponentId() == 0;
  };
#endif

  // Check that all nodes are component zero before we start
  LANCET_ASSERT(static_cast<usize>(std::ranges::count_if(mNodes, IS_UNASSIGNED)) == mNodes.size())

  for (NodeTable::reference item : mNodes) {
    if (item.second->GetComponentId() != 0) continue;

    current_component++;
    results_info.emplace_back(ComponentInfo{.mCompId = current_component, .mNumNodes = 0});

    absl::chunked_queue<Node*, 128, 1024> connected_nodes;
    connected_nodes.push_back(item.second.get());

    while (!connected_nodes.empty()) {
      auto* current_node = connected_nodes.front();
      LANCET_ASSERT(current_node != nullptr)

      if (current_node->GetComponentId() != 0) {
        connected_nodes.pop_front();
        continue;
      }

      current_node->SetComponentId(current_component);
      results_info[current_component - 1].mNumNodes += 1;
      for (Edge const& edge : *current_node) {
        auto const neighbour_itr = mNodes.find(edge.DstId());
        LANCET_ASSERT(neighbour_itr != mNodes.end())
        LANCET_ASSERT(neighbour_itr->second != nullptr)
        connected_nodes.push_back(neighbour_itr->second.get());
      }

      connected_nodes.pop_front();
    }
  }

  auto const total_num_nodes = static_cast<f64>(mNodes.size());
  for (auto& cinfo : results_info) {
    cinfo.mPctNodes = 100.0 * (static_cast<f64>(cinfo.mNumNodes) / total_num_nodes);
  }

  std::ranges::sort(results_info, [](ComponentInfo const& lhs, ComponentInfo const& rhs) -> bool {
    return lhs.mNumNodes > rhs.mNumNodes;
  });

  // Check that none of the nodes are component zero after we are done
  LANCET_ASSERT(static_cast<usize>(std::ranges::count_if(mNodes, IS_UNASSIGNED)) == 0)
  return results_info;
}

// ============================================================================
// Phase 3: Anchor Discovery
// ============================================================================

auto Graph::FindSource(usize const component_id) const -> RefAnchor {
  RefAnchor result{.mAnchorId = 0, .mRefOffset = 0, .mFoundAnchor = false};

  for (usize ref_idx = 0; ref_idx < mRefNodeIds.size(); ++ref_idx) {
    auto const itr = mNodes.find(mRefNodeIds[ref_idx]);
    if (itr == mNodes.end()) continue;

    LANCET_ASSERT(itr->second != nullptr)
    if (itr->second->GetComponentId() != component_id ||
        itr->second->TotalReadSupport() < mParams.mMinAnchorCov) {
      continue;
    }

    result.mAnchorId = itr->first;
    result.mRefOffset = ref_idx;
    result.mFoundAnchor = true;
    break;
  }

  return result;
}

auto Graph::FindSink(usize const component_id) const -> RefAnchor {
  RefAnchor result{.mAnchorId = 0, .mRefOffset = 0, .mFoundAnchor = false};

  for (i64 ref_idx = static_cast<i64>(mRefNodeIds.size() - 1); ref_idx >= 0; --ref_idx) {
    auto const itr = mNodes.find(mRefNodeIds[ref_idx]);
    if (itr == mNodes.end()) continue;

    LANCET_ASSERT(itr->second != nullptr)
    if (itr->second->GetComponentId() != component_id ||
        itr->second->TotalReadSupport() < mParams.mMinAnchorCov) {
      continue;
    }

    result.mAnchorId = itr->first;
    result.mRefOffset = static_cast<usize>(ref_idx);
    result.mFoundAnchor = true;
    break;
  }

  return result;
}

// ============================================================================
// Phase 4: Pruning + Compression
// ============================================================================

void Graph::PruneComponent(usize const component_id) {
  // Pruning pipeline: compress linear chains, remove low-coverage nodes, remove tips.
  // NOTE: The pre-compression HasCycle call (old line 103) has been removed. Graph
  // compression only merges degree-2 linear chain nodes — it cannot introduce or
  // remove cycles. Therefore cycle detection before compression was redundant with
  // the cycle detection after compression. We now check only once, on the smaller
  // (compressed) graph, which is faster.
  CompressGraph(component_id);
  WriteDotDevelop(GraphState::FIRST_COMPRESSION, component_id);
  RemoveLowCovNodes(component_id);
  WriteDotDevelop(GraphState::SECOND_LOW_COV_REMOVAL, component_id);
  CompressGraph(component_id);
  WriteDotDevelop(GraphState::SECOND_COMPRESSION, component_id);
  RemoveTips(component_id);
  WriteDotDevelop(GraphState::SHORT_TIP_REMOVAL, component_id);
}

// ============================================================================
// CompressGraph — Unitig compaction on the bidirected de Bruijn graph.
//
// Merges maximal linear chains (unitigs) into single nodes, reducing graph
// size without changing the walk space. A node is compressible when it has
// exactly one edge in a given sign-direction and its neighbor satisfies the
// buddy constraints (see IsPotentialBuddyEdge).
//
// Each node is compressed in both sign-directions (DEFAULT and OPPOSITE)
// to capture unitigs that extend in either orientation of the bidirected
// graph. Absorbed nodes are collected in `remove_nids` and deleted after
// the full pass to avoid iterator invalidation.
//
// Reference: BCALM2 bidirected graph compaction model
// https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
// ============================================================================
void Graph::CompressGraph(usize const component_id) {
  absl::flat_hash_set<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());

  for (NodeTable::const_reference item : mNodes) {
    if (item.second->GetComponentId() != component_id) continue;
    if (remove_nids.contains(item.first)) continue;

    CompressNode(item.first, Kmer::Ordering::DEFAULT, remove_nids);
    CompressNode(item.first, Kmer::Ordering::OPPOSITE, remove_nids);
  }

  if (!remove_nids.empty()) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    auto const region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Compressed {} nodes for {} in comp{} with k={}", remove_nids.size(), region_str,
              component_id, mCurrK)
    for (auto const nid : remove_nids) RemoveNode(mNodes.find(nid));
  }
}

// ============================================================================
// CompressNode — Extend a single node by absorbing unitig-interior neighbors.
//
// Starting from `nid`, repeatedly finds a compressible neighbor in the
// given sign-direction (`ord`) and merges it into `nid`:
//
//    BEFORE:  nid ──e₁──▶ buddy ──e₂──▶ next
//    AFTER:   nid′ ──────────────e₂′──▶ next    (buddy absorbed, edges rewired)
//
// Merge appends the buddy's unique suffix (or prefix, depending on edge
// signs) to nid's sequence and sums their coverage. The absorbed buddy's
// edges to its other neighbors are rewired to point directly to nid.
//
// Edge sign propagation during rewiring:
//   When buddy is absorbed, each of buddy's outgoing edges (buddy→next)
//   must be rewritten as (nid→next). The new edge's SrcSign depends on
//   whether the buddy's internal sign-continuity was preserved or flipped:
//     - If e₁.DstSign == e₂.SrcSign: signs are continuous → keep e₁.SrcSign
//     - If e₁.DstSign != e₂.SrcSign: signs are flipped   → flip e₁.SrcSign
//   This follows from the BCALM2 walk sign-continuity rule:
//   e_i.toSign must equal e_{i+1}.fromSign for a valid walk.
// ============================================================================
void Graph::CompressNode(NodeID nid, Kmer::Ordering const ord, NodeIdSet& compressed_ids) const {
  auto const node_itr = mNodes.find(nid);
  LANCET_ASSERT(node_itr != mNodes.end())
  LANCET_ASSERT(node_itr->second != nullptr)

  auto compressible_edge = FindCompressibleEdge(*node_itr->second, ord);
  while (compressible_edge.has_value()) {
    Edge const src2obdy = compressible_edge.value();
    LANCET_ASSERT(src2obdy.SrcId() == nid)
    auto const obdy_itr = mNodes.find(src2obdy.DstId());
    LANCET_ASSERT(obdy_itr != mNodes.end())
    LANCET_ASSERT(obdy_itr->second != nullptr)

    node_itr->second->Merge(*(obdy_itr->second), src2obdy.Kind(), mCurrK);
    node_itr->second->EraseEdge(src2obdy);  // src -->X--> old_buddy

    // Rewire buddy's outgoing edges to point at src (the absorbing node).
    // Sign propagation: if buddy's internal signs are continuous with the
    // merge edge, keep the original SrcSign; if flipped, reverse it.
    auto const rev_src2obdy_src_sign = Kmer::RevSign(src2obdy.SrcSign());
    for (Edge const& obdy2nbdy : *(obdy_itr->second)) {
      if (obdy2nbdy == src2obdy.MirrorEdge()) continue;

      LANCET_ASSERT(!obdy2nbdy.IsSelfLoop())
      LANCET_ASSERT(obdy2nbdy.DstId() != node_itr->second->Identifier())

      auto const nbdy_itr = mNodes.find(obdy2nbdy.DstId());
      LANCET_ASSERT(nbdy_itr != mNodes.end())
      LANCET_ASSERT(nbdy_itr->second != nullptr)

      auto const ne_src_sign =
          src2obdy.DstSign() != obdy2nbdy.SrcSign() ? rev_src2obdy_src_sign : src2obdy.SrcSign();
      auto const src2nbdy =
          Edge({nid, obdy2nbdy.DstId()}, MakeFwdEdgeKind({ne_src_sign, obdy2nbdy.DstSign()}));

      node_itr->second->EmplaceEdge(src2nbdy);               // src --> new_buddy
      nbdy_itr->second->EmplaceEdge(src2nbdy.MirrorEdge());  // new_buddy --> src
      nbdy_itr->second->EraseEdge(obdy2nbdy.MirrorEdge());   // new_buddy -->X--> old_buddy
    }

    compressed_ids.insert(src2obdy.DstId());
    compressible_edge = FindCompressibleEdge(*node_itr->second, ord);
  }
}

// ============================================================================
// FindCompressibleEdge — Find a unitig-interior neighbor to absorb.
//
// In the bidirected de Bruijn graph, each edge carries signs at both
// endpoints: (src, buddy, SrcSign, DstSign). The `ord` parameter selects
// which sign-direction to probe from `src`:
//   DEFAULT  → edges where SrcSign == src.SignFor(DEFAULT)   (canonical)
//   OPPOSITE → edges where SrcSign == src.SignFor(OPPOSITE)  (flipped)
//
// Returns the edge (src→buddy) suitable for merging, or std::nullopt.
//
//    Bidirected topology around src (probing in DEFAULT direction):
//
//      ┌───────────┐       merge edge       ┌───────┐
//      │ opp_nbour │ ◄──(∓)── src ──(+)──►  │ buddy │
//      └───────────┘   opposite   DEFAULT   └───────┘
//                       arm       direction
//
//   UNITIG MERGE CONDITIONS:
//
//    (1) src has 1–2 out-edges and no self-loops.
//    (2) Exactly one edge from src in the `ord` direction → merge candidate.
//    (3) src is NOT the source or sink anchor.
//    (4) buddy is NOT the source or sink anchor.
//    (5) IsPotentialBuddyEdge(src, merge_edge) passes.
//    (6) If an opposite-arm edge exists, it must also pass
//        IsPotentialBuddyEdge — both arms must be well-formed.
//
//   ANCHOR PROTECTION ((3) and (4)):
//
//   Source and sink nodes define the reference coordinate boundaries.
//   Neither may participate in compression in any direction:
//
//     ─ If src IS an anchor: absorbing neighbors extends the anchor's
//       sequence past the reference boundary. MaxFlow walks terminate
//       at the anchor but emit its full (inflated) sequence, producing
//       false InDels in the MSA.
//
//     ─ If buddy IS an anchor: absorbing it destroys the walk
//       termination point for MaxFlow's source→sink path discovery.
// ============================================================================
auto Graph::FindCompressibleEdge(Node const& src, Kmer::Ordering const ord) const
    -> std::optional<Edge> {
  if (src.NumOutEdges() > 2 || src.NumOutEdges() == 0 || src.HasSelfLoop()) return std::nullopt;

  // Anchor guard: prevent source/sink from absorbing any neighbors.
  auto const [source_id, sink_id] = mSourceAndSinkIds;
  if (src.Identifier() == source_id || src.Identifier() == sink_id) return std::nullopt;

  auto const mergeable_edges = src.FindEdgesInDirection(ord);
  if (mergeable_edges.size() != 1) return std::nullopt;

  auto const potential_result_edge = mergeable_edges[0];

  // Anchor guard: prevent any node from absorbing the source/sink.
  if (potential_result_edge.DstId() == source_id || potential_result_edge.DstId() == sink_id) {
    return std::nullopt;
  }

  if (!IsPotentialBuddyEdge(src, potential_result_edge)) return std::nullopt;

  // If src has a second edge in the opposite direction, verify its buddy
  // constraints too — both arms of the unitig must be well-formed.
  auto const opp_dir_edges = src.FindEdgesInDirection(Kmer::RevOrdering(ord));
  if (opp_dir_edges.empty()) return potential_result_edge;
  if (opp_dir_edges.size() > 1) return std::nullopt;

  if (!IsPotentialBuddyEdge(src, opp_dir_edges[0])) return std::nullopt;

  return potential_result_edge;
}

// ============================================================================
// IsPotentialBuddyEdge — Verify that a neighbor qualifies as a unitig interior.
//
// Given the edge (src→nbour), checks whether `nbour` is a valid unitig-
// interior node that can be absorbed into `src`.
//
// Topology being validated:
//
//     ┌─────┐  conn=(s₁,d₁)  ┌───────┐  (s₂,d₂)    ┌─────┐
//     │ src │ ──────────────►│ nbour │ ──────────► │ nnb │
//     └─────┘                └───────┘             └─────┘
//                             ◄──────
//                            mirror of conn
//                            (d₁,s₁) = nbour→src
//
//   s₁,d₁ = SrcSign,DstSign of conn   (determines orientation at each end)
//   s₂,d₂ = signs of nbour's onward edge to nnb (if it exists)
//
// UNITIG-INTERIOR CONDITIONS:
//
//   (1) nbour has 1–2 out-edges and no self-loops.
//   (2) nbour has exactly one edge in the d₁ sign-direction, and it must
//       be the mirror of conn (i.e., nbour→src). This confirms the
//       BCALM2 bidirected mirror constraint.
//   (3) If nbour has a second edge (in the opposite sign-direction), it
//       continues the chain to a next neighbor `nnb`.
//       - That edge must NOT loop back to src (cycle rejection).
//       - nnb must have at most 2 out-edges (not a branch point).
//
// DEGENERATE PAIR REJECTION:
//
//     ┌─────┐                ┌───────┐
//     │ src │ ◄───────────── │ nbour │
//     │     │ ──────────────►│       │
//     └─────┘  both degree-1 └───────┘
//
//   If src and nbour each have exactly one edge pointing at each other,
//   merging would produce a degenerate zero-edge node. Rejected.
// ============================================================================
auto Graph::IsPotentialBuddyEdge(Node const& src, Edge const& conn) const -> bool {
  auto const nbour_itr = mNodes.find(conn.DstId());
  LANCET_ASSERT(nbour_itr != mNodes.end())
  LANCET_ASSERT(nbour_itr->second != nullptr)
  Node const& nbour = *nbour_itr->second;

  // Degenerate pair: src ↔ nbour with both degree-1. Reject.
  if (src.NumOutEdges() == 1 && nbour.NumOutEdges() == 1) {
    auto const& edge_from_src = *src.cbegin();
    auto const& edge_from_nbour = *nbour.cbegin();
    if (edge_from_src.DstId() == nbour.Identifier() &&
        edge_from_nbour.DstId() == src.Identifier()) {
      return false;
    }
  }

  if (nbour.NumOutEdges() > 2 || nbour.NumOutEdges() == 0 || nbour.HasSelfLoop()) return false;

  // (2) Verify the mirror edge (nbour→src) exists and is the sole edge
  //     in the d₁ sign-direction from nbour.
  auto const expected_nbour2src = conn.MirrorEdge();
  auto const start_sign_nbour2src = expected_nbour2src.SrcSign();
  auto const dir_nbour2src = start_sign_nbour2src == nbour.SignFor(Kmer::Ordering::DEFAULT)
                                 ? Kmer::Ordering::DEFAULT
                                 : Kmer::Ordering::OPPOSITE;
  auto const nb_edges_in_nbour2src_dir = nbour.FindEdgesInDirection(dir_nbour2src);
  if (nb_edges_in_nbour2src_dir.size() != 1 || nb_edges_in_nbour2src_dir[0] != expected_nbour2src) {
    return false;
  }

  // (3) Check the onward edge in the opposite sign-direction.
  //     Must not loop back to src (cycle) and nnb must have degree ≤ 2.
  auto const nb_edges_in_opp_dir = nbour.FindEdgesInDirection(Kmer::RevOrdering(dir_nbour2src));
  if (nb_edges_in_opp_dir.size() != 1 || nb_edges_in_opp_dir[0].DstId() == conn.SrcId()) {
    return false;
  }

  auto const nnb_itr = mNodes.find(nb_edges_in_opp_dir[0].DstId());
  LANCET_ASSERT(nnb_itr != mNodes.end())
  LANCET_ASSERT(nnb_itr->second != nullptr)
  return nnb_itr->second->NumOutEdges() <= 2;
}

void Graph::RemoveTips(usize const component_id) {
  usize total_tips = 0;
  usize curr_tips = 1;

  std::vector<NodeID> remove_nids;
  remove_nids.reserve(mNodes.size());
  // remove tips and compress at least once. compression after tip removal
  // can produce new tips in the graph, so recursively remove tips from
  // the graph until there are no longer any tips left
  while (curr_tips > 0) {
    remove_nids.clear();

    auto const [source_id, sink_id] = mSourceAndSinkIds;
    for (auto const& [nid, node_ptr] : mNodes) {
      if (node_ptr->GetComponentId() != component_id || node_ptr->NumOutEdges() > 1) continue;
      if (nid == source_id || nid == sink_id) continue;

      auto const uniq_sequence_length = node_ptr->SeqLength() - mCurrK + 1;
      if (uniq_sequence_length >= mCurrK) continue;

      remove_nids.emplace_back(nid);
    }

    if (!remove_nids.empty()) {
      total_tips += curr_tips;
      RemoveNodes(absl::MakeConstSpan(remove_nids));
      CompressGraph(component_id);
    }

    curr_tips = remove_nids.size();
  }

  if (total_tips > 0) {
    // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
    auto const region_str = mRegion->ToSamtoolsRegion();
    LOG_TRACE("Removed {} tips for {} in comp{} with k={}", total_tips, region_str, component_id,
              mCurrK);
  }
}

// ============================================================================
// Phase 5: Haplotype Enumeration
// ============================================================================

auto Graph::EnumerateAndSortHaplotypes(usize comp_id, TraversalIndex const& trav_idx,
                                       std::string_view ref_anchor_seq) const -> std::vector<Path> {
  std::vector<Path> haplotypes;
  auto const reg_str = mRegion->ToSamtoolsRegion();

  LOG_TRACE("Starting walk enumeration for {} with k={}, num_nodes={}", reg_str, mCurrK,
            mNodes.size())

  MaxFlow max_flow(&mNodes, mSourceAndSinkIds, mCurrK, &trav_idx);
  auto path_seq = max_flow.NextPath();

  while (path_seq) {
    LOG_DEBUG("Assembled {}bp path sequence for {} comp={} with k={}",
              path_seq->Sequence().length(), reg_str, comp_id, mCurrK)
    haplotypes.emplace_back(std::move(*path_seq));
    path_seq = max_flow.NextPath();
  }

  if (!haplotypes.empty()) {
    // Sort non-reference ALT haplotypes strictly by descending mean coverage.
    // This ensures downstream Greedy Insertion Bias in SPOA prioritizes dominant somatic signals.
    std::ranges::sort(haplotypes, [](Path const& lhs, Path const& rhs) -> bool {
      return lhs.MeanCoverage() > rhs.MeanCoverage();
    });

    // Deduplicate Sequence matches in O(N). Because the array is already sorted
    // by coverage, this retains the highest-coverage path for any duplicate sequence.
    absl::flat_hash_set<std::string_view> seen_seqs;
    std::erase_if(haplotypes, [&seen_seqs](Path const& path) -> bool {
      auto [_unused, inserted] = seen_seqs.insert(path.Sequence());
      return !inserted;
    });

    Path ref_path;
    ref_path.AppendSequence(ref_anchor_seq);
    haplotypes.emplace(haplotypes.begin(), std::move(ref_path));
  }

  return haplotypes;
}

// ============================================================================
// Debug DOT Visualization
// ============================================================================

void Graph::WriteDot([[maybe_unused]] GraphState state, usize comp_id) {
  if (mParams.mOutGraphsDir.empty()) return;

#ifdef LANCET_DEVELOP_MODE
  auto const graph_state = ToString(state);
#else
  auto const* graph_state = "fully_pruned";
#endif

  using namespace std::string_view_literals;
  auto const win_id =
      fmt::format("{}_{}_{}", mRegion->ChromName(), mRegion->StartPos1(), mRegion->EndPos1());
  auto const fname =
      fmt::format("dbg__{}__{}__k{}__comp{}.dot", win_id, graph_state, mCurrK, comp_id);

  auto const out_path = mParams.mOutGraphsDir / "dbg_graph" / fname;
  std::filesystem::create_directories(mParams.mOutGraphsDir / "dbg_graph");
  DotOverlaySets const highlight{
      .mNodes = {mSourceAndSinkIds[0], mSourceAndSinkIds[1]},
  };
  SerializeToDot(mNodes, out_path, comp_id, highlight);
}

}  // namespace lancet::cbdg
