#ifndef SRC_LANCET_CBDG_GRAPH_H_
#define SRC_LANCET_CBDG_GRAPH_H_

#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/component_result.h"
#include "lancet/cbdg/dot_snapshot_buffer.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph_complexity.h"
#include "lancet/cbdg/graph_params.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"
#include "lancet/cbdg/probe_tracker.h"
#include "lancet/cbdg/read.h"
#include "lancet/cbdg/traversal_index.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeTable = absl::flat_hash_map<NodeID, NodePtr>;
  using RegionPtr = std::shared_ptr<hts::Reference::Region const>;
  using ReadList = absl::Span<Read const>;
  using ComponentResults = std::vector<ComponentResult>;

  explicit Graph(GraphParams params) : mParams(params) {}

  /// Current kmer length being used for this assembly attempt.
  [[nodiscard]] auto CurrentK() const noexcept -> usize { return mCurrK; }

  /// Const access to the node table — used by graph_complexity.cpp free functions.
  [[nodiscard]] auto Nodes() const noexcept -> NodeTable const& { return mNodes; }

  /// Main entry point: build, prune, and enumerate haplotypes from reads + reference.
  /// Iterates kmer lengths from min to max, returning per-component results on success.
  [[nodiscard]] auto BuildComponentResults(RegionPtr region, ReadList reads) -> ComponentResults;

  /// Set the external ProbeTracker for truth variant k-mer tracing. Null
  /// disables tracing (zero overhead in production).
  void SetProbeTracker(ProbeTracker* tracker) { mProbeTrackerPtr = tracker; }

  /// Set the external per-worker TarGzWriter shard that buffered DOTs are
  /// appended to. Null disables snapshot emission entirely (used when
  /// `--out-graphs-tgz` is unset; zero overhead in production).
  void SetGraphShardWriter(base::TarGzWriter* shard_writer) noexcept {
    mGraphShardWriter = shard_writer;
  }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  usize mCurrK = 0;
  RegionPtr mRegion;
  ReadList mReads;
  NodeTable mNodes;
  GraphParams mParams;

  /// Non-owning pointer to the ProbeTracker owned by ProbeDiagnostics.
  /// Null when --probe-variants is not specified (zero overhead in production).
  ProbeTracker* mProbeTrackerPtr = nullptr;

  /// Non-owning pointer to the per-worker TarGzWriter shard owned by
  /// VariantBuilder. Null when `--out-graphs-tgz` is not specified (zero
  /// overhead in production). On Commit, mDotBuffer appends each
  /// buffered DOT into this shard as a regular-file TAR entry.
  base::TarGzWriter* mGraphShardWriter = nullptr;

  std::vector<NodeID> mRefNodeIds;
  NodeIDPair mSourceAndSinkIds = {0, 0};

  /// In-memory accumulator for the per-component DOT snapshots emitted
  /// during the current k-attempt. Discarded on retry. On a successful
  /// k-attempt, the buffer drains its contents into the per-worker
  /// tar.gz shard via DotSnapshotBuffer::Commit.
  DotSnapshotBuffer mDotBuffer;

  using EdgeSet = absl::flat_hash_set<Edge>;
  using NodeIdSet = absl::flat_hash_set<NodeID>;

  // ============================================================================
  // Phase 1: Graph Construction
  // ============================================================================

  // De-duplicate read support: a kmer should receive at most one increment per
  // (read_name, sample_tag, kmer_hash) triple. MateMer keys this dedup set.
  struct MateMer {
    // ── 8B Align ────────────────────────────────────────────────────────────
    std::string_view mQname;  // 16B
    u64 mKmerHash;            // 8B
    // ── 1B Align ────────────────────────────────────────────────────────────
    Label::Tag mTagKind;  // 1B

    auto operator==(MateMer const& rhs) const -> bool {
      return mQname == rhs.mQname && mTagKind == rhs.mTagKind && mKmerHash == rhs.mKmerHash;
    }

    template <typename H>
    friend auto AbslHashValue(H state, MateMer const& mmer) -> H {
      return H::combine(std::move(state), mmer.mQname, mmer.mTagKind, mmer.mKmerHash);
    }
  };

  /// Construct the de Bruijn graph from reference + read sequences at current k.
  void BuildGraph(absl::flat_hash_set<MateMer>& mate_mers);

  /// Insert overlapping k+1-mers from a sequence, creating nodes and edges.
  auto AddNodes(std::string_view sequence, Label label) -> std::vector<Node*>;

  /// True if the reference sequence contains a repeated k-mer (exact or approximate),
  /// which would create a cycle by construction — making assembly at this k pointless.
  [[nodiscard]] static auto HasExactOrApproxRepeat(std::string_view seq, usize window) -> bool {
    auto const klen_seqs = lancet::base::SlidingView(seq, window);
    static constexpr usize NUM_ALLOWED_MISMATCHES = 2;
    return lancet::base::HasRepeat(absl::MakeConstSpan(klen_seqs), NUM_ALLOWED_MISMATCHES);
  }

  // ============================================================================
  // Phase 2: Node Removal + Connected Components
  // ============================================================================

  /// Erase a single node and all its incoming edges from the graph.
  void RemoveNode(NodeTable::iterator itr);

  /// Batch-remove nodes by ID (calls RemoveNode for each).
  void RemoveNodes(absl::Span<NodeID const> node_ids) {
    for (auto const nid : node_ids) RemoveNode(mNodes.find(nid));
  }

  /// Remove nodes with singleton or below-threshold coverage (sequencing errors).
  void RemoveLowCovNodes(usize component_id);

  /// Per-component metadata for sorting and filtering after BFS labeling.
  struct ComponentInfo {
    // ── 8B Align ────────────────────────────────────────────────────────────
    f64 mPctNodes = 0.0;  // 8B
    usize mCompId = 0;    // 8B
    usize mNumNodes = 0;  // 8B
  };

  /// BFS-label all nodes into connected components, returning component metadata
  /// sorted by descending node count.
  [[nodiscard]] auto MarkConnectedComponents() -> std::vector<ComponentInfo>;

  // ============================================================================
  // Phase 3: Anchor Discovery
  // ============================================================================

  /// A reference k-mer that anchors assembly: provides the start (source) or
  /// end (sink) of the reference path through the graph component.
  struct RefAnchor {
    // ── 8B Align ────────────────────────────────────────────────────────────
    NodeID mAnchorId = 0;  // 8B
    usize mRefOffset = 0;  // 8B
    // ── 1B Align ────────────────────────────────────────────────────────────
    bool mFoundAnchor = false;  // 1B
  };

  /// Walk mRefNodeIds left-to-right, returning the first reference k-mer in the component
  /// with sufficient coverage. This is the source anchor for walk enumeration.
  [[nodiscard]] auto FindSource(usize component_id) const -> RefAnchor;
  /// Walk mRefNodeIds right-to-left, returning the last reference k-mer in the component
  /// with sufficient coverage. This is the sink anchor for walk enumeration.
  [[nodiscard]] auto FindSink(usize component_id) const -> RefAnchor;

  /// Length of the reference anchor region: sink_offset - source_offset + k.
  [[nodiscard]] static auto RefAnchorLength(RefAnchor const& source, RefAnchor const& sink,
                                            usize currk) -> usize {
    return sink.mRefOffset - source.mRefOffset + currk;
  }

  // ============================================================================
  // Phase 4: Pruning + Compression
  // ============================================================================

  /// Run the full pruning pipeline on a component: compress → low-cov remove → compress → tips.
  void PruneComponent(usize component_id);

  /// Merge degree-2 linear chain nodes into their neighbors, reducing graph size.
  void CompressGraph(usize component_id);

  /// Merge a single node's compressible neighbor(s) in the given ordering direction.
  void CompressNode(NodeID nid, Kmer::Ordering ord, NodeIdSet& compressed_ids) const;

  /// Find an edge from `src` that can be compressed (merged into src) in the given direction.
  [[nodiscard]] auto FindCompressibleEdge(Node const& src, Kmer::Ordering ord) const
      -> std::optional<Edge>;

  /// Check whether `conn`'s destination qualifies as a buddy node for compression.
  [[nodiscard]] auto IsPotentialBuddyEdge(Node const& src, Edge const& conn) const -> bool;

  /// Iteratively remove short dead-end branches (tips) and re-compress until stable.
  void RemoveTips(usize component_id);

  // ============================================================================
  // Phase 5: Complexity Analysis + Haplotype Enumeration
  // ============================================================================

  // ComputeGraphComplexity — O(V+E) Metrics for Debug Logging
  //
  // Computes lightweight graph topology metrics that correlate with runtime:
  //
  //  Cyclomatic complexity (M = E - V + 1): number of independent cycles.
  //    M=0 → linear chain, M=1 → single variant bubble, M>>1 → STR hairball.
  //
  //  Edge-to-node density (E/V): a clean graph has E/V ≈ 1.0. Values >1.5
  //    indicate branching that causes path enumeration explosion.
  //
  //  Max single-direction degree: maximum outgoing edges in any one sign
  //    direction. Hub nodes with high degree are direct BFS blowup predictors.
  //
  //  Branch points: nodes with ≥2 outgoing edges in some direction. These
  //    are the "decision points" that cause combinatorial path blowup.
  //
  [[nodiscard]] auto ComputeComponentComplexity(usize component_id) const -> GraphComplexity {
    return cbdg::ComputeGraphComplexity(*this, component_id);
  }

  /// Enumerate all source→sink walks via MaxFlow BFS, sort ALT haplotypes by
  /// descending MinWeight (weakest-link confidence), deduplicate, and prepend
  /// a confidence-weighted reference haplotype. Each result bundles the
  /// assembled `Path` with its underlying source→sink edge walk; the walk is
  /// empty for the REF haplotype when reconstruction fails.
  [[nodiscard]] auto BuildHaplotypes(usize comp_id, TraversalIndex const& trav_idx,
                                     std::string_view ref_anchor_seq,
                                     ProbeTracker::Context const& probe_ctx) const
      -> std::vector<EnumeratedHaplotype>;

  /// Build the reference haplotype, bundling the median-confidence-weighted
  /// `Path` with the surviving REF backbone walk traced through the
  /// post-prune compressed graph. Returned walk is empty if the surviving
  /// backbone is fragmented.
  [[nodiscard]] auto BuildRefHaplotype(usize comp_id, std::string_view ref_anchor_seq) const
      -> EnumeratedHaplotype;

  // ============================================================================
  // Debug DOT Visualization
  // ============================================================================

  /// Render the post-prune ("final") DOT for one component into mDotBuffer.
  /// `haplotypes` non-empty → file basename uses `enumerated_walks` and per-
  /// walk edge color overlays are appended; empty → basename uses
  /// `fully_pruned` and only the graph body is rendered. The actual disk
  /// write is deferred until mDotBuffer.Commit is called after the outer
  /// k-loop succeeds.
  void BufferFinalSnapshot(usize comp_id, absl::Span<EnumeratedHaplotype const> haplotypes);

  /// Render an intermediate-stage DOT snapshot (post-pruning boundary) into
  /// mDotBuffer. No-op unless `mParams.mSnapshotMode == VERBOSE` and an
  /// output directory is configured. Replaces the old `WriteDotDevelop`
  /// template + `LANCET_DEVELOP_MODE` compile-time gate with a runtime
  /// check so the verbose mode is opt-in via `--graph-snapshots=verbose`.
  void BufferStageSnapshot(PruneStage stage, usize comp_id);

  // ============================================================================
  // ProbeTracker Forwarding Helpers (null-safe)
  //
  // Each helper checks mProbeTrackerPtr for null before forwarding.
  // In production (no --probe-variants), the pointer is null and all
  // calls are no-ops with zero overhead beyond a single branch prediction.
  // ============================================================================

  using Context = ProbeTracker::Context;

  [[nodiscard]] auto HasProbeTracker() const -> bool { return mProbeTrackerPtr != nullptr; }

  // Graph construction: tag ALT k-mers and count them in raw reads.
  void ProbeGenerateAndTag(Context const& ctx);
  void ProbeCountInReads(Context const& ctx);
  void ProbeLogStatus(PruneStage stage, Context const& ctx);

  // Graph mutation: keep side-table synchronized with node removal/merging.
  void ProbeOnNodeRemove(NodeID node_id);
  void ProbeOnNodeMerge(NodeID absorbed_id, NodeID surviving_id) const;

  // Component analysis: record component membership and anchor failures.
  using ProbeCompInfo = ProbeTracker::ComponentInfo;
  void ProbeRecordComponentInfo(absl::Span<ProbeCompInfo const> probe_comps, Context const& ctx);
  void ProbeSetNoAnchor(Context const& ctx);
  void ProbeSetShortAnchor(Context const& ctx);
  void ProbeCheckAnchorOverlap(RefAnchor const& source, RefAnchor const& sink, Context const& ctx);

  // Graph structure: flag probes when assembly encounters cycles or complexity.
  void ProbeSetGraphCycle(Context const& ctx);
  void ProbeSetGraphComplex(Context const& ctx);

  // Path enumeration: flag traversal-limited probes and check path survival.
  void ProbeSetTraversalLimit(Context const& ctx) const;
  void ProbeCheckPaths(absl::Span<EnumeratedHaplotype const> haplotypes, Context const& ctx);
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_H_
