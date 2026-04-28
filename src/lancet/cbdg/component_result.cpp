#include "lancet/cbdg/component_result.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/path.h"

#include <algorithm>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

ComponentResult::ComponentResult(std::vector<EnumeratedHaplotype> haplotypes,
                                 GraphComplexity metrics, u32 const anchor_start_offset)
    : mMetrics(metrics), mAnchorStartOffset(anchor_start_offset) {
  mPaths.reserve(haplotypes.size());
  mWalks.reserve(haplotypes.size());
  for (auto& hap : haplotypes) {
    mPaths.emplace_back(std::move(hap.mPath));
    mWalks.emplace_back(std::move(hap.mWalk));
  }
}

auto ComponentResult::HaplotypeSequenceViews() const -> std::vector<std::string_view> {
  std::vector<std::string_view> views;
  views.reserve(mPaths.size());
  std::ranges::transform(mPaths, std::back_inserter(views),
                         [](auto const& path) -> std::string_view { return path.Sequence(); });
  return views;
}

auto ComponentResult::HaplotypeSequences() const -> std::vector<std::string> {
  std::vector<std::string> seqs;
  seqs.reserve(mPaths.size());
  std::ranges::transform(mPaths, std::back_inserter(seqs), [](auto const& path) -> std::string {
    return std::string(path.Sequence());
  });
  return seqs;
}

auto ComponentResult::HaplotypeWeights() const -> Path::HapWeights {
  Path::HapWeights weights;
  weights.reserve(mPaths.size());
  for (auto const& path : mPaths) {
    weights.emplace_back(path.PerBaseWeights());
  }
  return weights;
}

auto ComponentResult::MaxAltPathCv() const -> std::optional<f64> {
  std::optional<f64> max_cv;
  for (usize hap = 1; hap < mPaths.size(); ++hap) {
    auto const cv_val = mPaths[hap].CoefficientOfVariationCoverage();
    max_cv = max_cv.has_value() ? std::max(*max_cv, cv_val) : cv_val;
  }
  return max_cv;
}

}  // namespace lancet::cbdg
