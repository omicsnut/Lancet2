#ifndef SRC_LANCET_CBDG_PATH_H_
#define SRC_LANCET_CBDG_PATH_H_

#include "lancet/base/types.h"

#include "absl/container/inlined_vector.h"

#include <string>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

class Path {
 public:
  Path() = default;

  /// Concatenate a sequence fragment (k-mer or overlap extension) to the path.
  void AppendSequence(std::string_view seq);
  /// Pre-allocate storage for the expected total sequence length.
  void ReserveSequence(usize len) { mSequence.reserve(len); }
  /// Record one node's total read coverage for downstream statistics.
  void AddNodeCoverage(u32 cov);
  /// Compute summary statistics (mean, median, CV, QCV) from accumulated node coverages.
  void Finalize();

  [[nodiscard]] auto IsEmpty() const -> bool { return mSequence.empty(); }
  [[nodiscard]] auto Sequence() const -> std::string_view { return mSequence; }
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
  std::string mSequence;
  absl::InlinedVector<u32, 256> mNodeCoverages;
  f64 mMeanCov = 0.0;
  f64 mMedianCov = 0.0;
  f64 mStdDevCov = 0.0;
  f64 mCvCov = 0.0;
  f64 mQCvCov = 0.0;
  f64 mTotalCov = 0.0;
};

// First is always ref hap. Subsequent ALT haplotypes are sorted by descending
// MeanCoverage, establishing Greedy Insertion Bias in downstream SPOA MSA.
using CompHaps = std::vector<Path>;
using GraphHaps = std::vector<CompHaps>;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PATH_H_
