#ifndef SRC_LANCET_CBDG_GRAPH_H_
#define SRC_LANCET_CBDG_GRAPH_H_

#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph_complexity.h"
#include "lancet/cbdg/graph_params.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"
#include "lancet/cbdg/read.h"
#include "lancet/cbdg/serialize_dot.h"
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

  explicit Graph(GraphParams params) : mParams(std::move(params)) {}

  /// Current kmer length being used for this assembly attempt.
  [[nodiscard]] auto CurrentK() const noexcept -> usize { return mCurrK; }

  /// Const access to the node table — used by graph_complexity.cpp free functions.
  [[nodiscard]] auto Nodes() const noexcept -> NodeTable const& { return mNodes; }

  /// Output of the main assembly pipeline: assembled haplotypes, anchor offsets,
  /// and per-component complexity metrics for downstream variant calling.
  struct Result {
    GraphHaps mGraphHaplotypes;
    std::vector<usize> mAnchorStartIdxs;
    std::vector<GraphComplexity> mComponentMetrics;
  };

  /// Main entry point: build, prune, and enumerate haplotypes from reads + reference.
  /// Iterates kmer lengths from min to max, returning assembled haplotypes on success.
  [[nodiscard]] auto BuildComponentHaplotypes(RegionPtr region, ReadList reads) -> Result;

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  usize mCurrK = 0;
  RegionPtr mRegion;
  ReadList mReads;
  NodeTable mNodes;
  GraphParams mParams;

  std::vector<NodeID> mRefNodeIds;
  NodeIDPair mSourceAndSinkIds = {0, 0};

  using EdgeSet = absl::flat_hash_set<Edge>;
  using NodeIdSet = absl::flat_hash_set<NodeID>;

  // ============================================================================
  // Phase 1: Graph Construction
  // ============================================================================

  // De-duplicate read support: a kmer should receive at most one increment per
  // (read_name, sample_tag, kmer_hash) triple. MateMer keys this dedup set.
  struct MateMer {
    std::string_view mQname;  // 16B
    u64 mKmerHash;            // 8B
    Label::Tag mTagKind;      // 1B

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
    static constexpr usize NUM_ALLOWED_MISMATCHES = 3;
    return lancet::base::HasRepeat(absl::MakeConstSpan(klen_seqs), NUM_ALLOWED_MISMATCHES);
  }

  // ============================================================================
  // Phase 2: Node Removal + Connected Components
  // ============================================================================

  /// Erase a single node and all its incoming edges from the graph.
  void RemoveNode(NodeTable::iterator itr);

  /// Batch-remove nodes by ID (calls RemoveNode for each).
  void RemoveNodes(absl::Span<NodeID const> node_ids) {
    for (auto const nid : node_ids) {
      RemoveNode(mNodes.find(nid));
    }
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
    NodeID mAnchorId = 0;       // 8B
    usize mRefOffset = 0;       // 8B
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
  /// descending mean coverage, deduplicate, and prepend the reference path.
  [[nodiscard]] auto EnumerateAndSortHaplotypes(usize comp_id, TraversalIndex const& trav_idx,
                                                std::string_view ref_anchor_seq) const
      -> std::vector<Path>;

  // ============================================================================
  // Debug DOT Visualization
  // ============================================================================

  /// Write the current graph state to a Graphviz DOT file for visual debugging.
  void WriteDot([[maybe_unused]] GraphState state, usize comp_id);

#ifdef LANCET_DEVELOP_MODE
  template <class... Args>
  constexpr void WriteDotDevelop(Args&&... args) {
    WriteDot(std::forward<Args>(args)...);
  }
#else
  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  constexpr void WriteDotDevelop([[maybe_unused]] Args&&... /*unused*/) {}
#endif
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_H_
