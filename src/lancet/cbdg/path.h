#ifndef SRC_LANCET_CBDG_PATH_H_
#define SRC_LANCET_CBDG_PATH_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"

#include "absl/container/inlined_vector.h"

#include <string>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

class Path {
 public:
  /// Per-base SPOA weight vector. Each element is the confidence weight for
  /// one base position in the haplotype sequence string.
  using BaseWeights = std::vector<u32>;
  /// Collection of per-base weight vectors, one per haplotype path.
  using HapWeights = std::vector<BaseWeights>;

  Path() = default;

  /// Concatenate a sequence fragment (k-mer or overlap extension) to the path.
  void AppendSequence(std::string_view seq);
  /// Pre-allocate storage for the expected total sequence length.
  void ReserveSequence(usize len);
  /// Record one node's total read coverage for downstream statistics.
  void AddNodeCoverage(u32 cov);

  /// Record a node's confidence weight and the number of sequence bases it
  /// contributed to the haplotype string. Called alongside AddNodeCoverage
  /// during walk construction — one call per node in the walk.
  void AddNodeWeight(u32 weight, u32 num_bases);

  /// Compute summary statistics (mean, median, CV, QCV) from accumulated node coverages.
  void Finalize();

  [[nodiscard]] auto IsEmpty() const -> bool { return mSequence.empty(); }
  [[nodiscard]] auto Sequence() const -> std::string_view { return mSequence; }

  /// Expand per-node weights into per-base weights for SPOA.
  /// Returns a vector of size == Sequence().size().
  [[nodiscard]] auto PerBaseWeights() const -> BaseWeights;

  /// Weakest-link: minimum node confidence across the entire path.
  /// A path is only as trustworthy as its least-supported node.
  [[nodiscard]] auto MinWeight() const -> u32;

  /// Average read coverage across all nodes constituting this path
  [[nodiscard]] auto MeanCoverage() const -> f64 { return mMeanCov; }
  /// Median read coverage across all nodes constituting this path
  [[nodiscard]] auto MedianCoverage() const -> f64 { return mMedianCov; }
  /// Standard deviation of read coverage across nodes
  [[nodiscard]] auto StdDevCoverage() const -> f64 { return mStdDevCov; }
  /// Coefficient of Variation (StdDev / Mean) -> higher means more coverage fluctuation
  [[nodiscard]] auto CoefficientOfVariationCoverage() const -> f64 { return mCvCov; }
  /// Quartile Coefficient of Variation ((Q3 - Q1) / (Q3 + Q1)) -> robust against outliers
  [[nodiscard]] auto QuartileCoefficientOfVariation() const -> f64 { return mQCvCov; }
  /// Aggregate total of all node coverages on this path
  [[nodiscard]] auto TotalCoverage() const -> f64 { return mTotalCov; }

 private:
  /// Per-node confidence weight entry. Stored during walk construction,
  /// expanded lazily into per-base weights when SPOA needs them.
  struct NodeWeightEntry {
    // ── 4B Align ──────────────────────────────────────────────────────────
    u32 mWeight;    // 4B — Node::Confidence() value
    u32 mNumBases;  // 4B — sequence bases this node contributed
  };

  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string mSequence;
  absl::InlinedVector<u32, 256> mNodeCoverages;
  absl::InlinedVector<NodeWeightEntry, 256> mNodeWeights;
  f64 mMeanCov = 0.0;
  f64 mMedianCov = 0.0;
  f64 mStdDevCov = 0.0;
  f64 mCvCov = 0.0;
  f64 mQCvCov = 0.0;
  f64 mTotalCov = 0.0;
};

/// Bundles an assembled haplotype `Path` (DNA sequence + per-node weights and
/// coverage statistics) with the underlying source→sink edge walk that
/// produced it. The walk is the bidirected-graph edge sequence that MaxFlow
/// traversed; downstream visualization consumes it to overlay walk colors on
/// the rendered DOT graph. Empty walk denotes the reference haplotype that
/// was not produced by a MaxFlow traversal (see `Graph::BuildRefHaplotypePath`).
struct EnumeratedHaplotype {
  // ── 8B Align ────────────────────────────────────────────────────────────
  Path mPath;
  std::vector<Edge> mWalk;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PATH_H_
