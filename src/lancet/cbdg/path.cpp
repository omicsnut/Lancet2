#include "lancet/cbdg/path.h"

#include "lancet/base/compute_stats.h"
#include "lancet/base/types.h"

#include "absl/strings/str_cat.h"
#include "absl/types/span.h"

#include <absl/container/inlined_vector.h>
#include <algorithm>
#include <string_view>

namespace lancet::cbdg {

void Path::AppendSequence(std::string_view seq) {
  absl::StrAppend(&mSequence, seq);
}

void Path::AddNodeCoverage(u32 cov) {
  mNodeCoverages.push_back(cov);
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
