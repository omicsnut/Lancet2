#include "lancet/cbdg/path.h"

#include "lancet/base/compute_stats.h"
#include "lancet/base/types.h"

#include "absl/container/inlined_vector.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"

#include <algorithm>
#include <functional>
#include <ranges>
#include <string_view>

namespace lancet::cbdg {

void Path::AppendSequence(std::string_view seq) {
  absl::StrAppend(&mSequence, seq);
}

void Path::ReserveSequence(usize const len) {
  mSequence.reserve(len);
  mNodeCoverages.reserve(len);
  mNodeWeights.reserve(len);
}

void Path::AddNodeCoverage(u32 const cov) {
  mNodeCoverages.push_back(cov);
}

void Path::AddNodeWeight(u32 const weight, u32 const num_bases) {
  mNodeWeights.push_back({.mWeight = weight, .mNumBases = num_bases});
}

void Path::AddWalkNodeId(NodeID const node_id) {
  mWalkNodeIds.push_back(node_id);
}

auto Path::WalkNodeIds() const -> absl::Span<NodeID const> {
  return absl::MakeConstSpan(mWalkNodeIds);
}

auto Path::PerBaseWeights() const -> BaseWeights {
  BaseWeights per_base;
  per_base.reserve(mSequence.size());
  for (auto const& [weight, num_bases] : mNodeWeights) {
    per_base.insert(per_base.end(), num_bases, weight);
  }
  return per_base;
}

auto Path::MinWeight() const -> u32 {
  if (mNodeWeights.empty()) return 0;
  return std::ranges::min(mNodeWeights, {}, &NodeWeightEntry::mWeight).mWeight;
}

void Path::Finalize() {
  if (mNodeCoverages.empty()) return;

  lancet::base::OnlineStats stats;
  for (auto const cov : mNodeCoverages) {
    stats.Add(cov);
  }

  mMeanCov = stats.Mean();
  mStdDevCov = stats.StdDev();
  mTotalCov = mMeanCov * static_cast<f64>(stats.Count());
  if (mMeanCov > 0.0) mCvCov = mStdDevCov / mMeanCov;

  mMedianCov = static_cast<f64>(lancet::base::Median(absl::MakeConstSpan(mNodeCoverages)));

  if (mNodeCoverages.size() >= 4) {
    std::ranges::sort(mNodeCoverages);
    auto const quart1 = static_cast<f64>(mNodeCoverages[mNodeCoverages.size() / 4]);
    auto const quart3 = static_cast<f64>(mNodeCoverages[(mNodeCoverages.size() * 3) / 4]);
    // Quartile CV = (Q3 - Q1) / (Q3 + Q1)
    if (quart3 + quart1 > 0.0) mQCvCov = (quart3 - quart1) / (quart3 + quart1);
  }
}

}  // namespace lancet::cbdg
