#include "lancet/cbdg/node.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>

namespace lancet::cbdg {

void Node::AddLabel(Label const& label) {
  mLabel.Merge(label);
}

void Node::IncrementReadSupport(usize const sample_index, Label::Tag const tag) {
  if (sample_index >= mCounts.size()) {
    mCounts.resize(sample_index + 1, 0);
  }
  mCounts[sample_index] += 1;
  mRoleCounts[RoleIndex(tag)] += 1;
}

auto Node::ReadSupportForSample(usize const sample_index) const -> u32 {
  return sample_index < mCounts.size() ? mCounts[sample_index] : 0;
}

auto Node::TotalReadSupport() const noexcept -> u32 {
  return std::reduce(mCounts.begin(), mCounts.end(), u32{0});
}

auto Node::ReadSupportForRole(Label::Tag const role) const noexcept -> u32 {
  return mRoleCounts[RoleIndex(role)];
}

auto Node::IsAllSingletons() const noexcept -> bool {
  // Every sample with coverage must have exactly 1 read, and at least one must have coverage.
  return std::ranges::any_of(mCounts, [](u32 const cnt) { return cnt > 0; }) &&
         std::ranges::all_of(mCounts, [](u32 const cnt) { return cnt <= 1; });
}

void Node::Merge(Node const& other, EdgeKind const conn_kind, usize const currk) {
  mKmer.Merge(other.mKmer, conn_kind, currk);
  mLabel.Merge(other.mLabel);

  // Weighted average of per-sample read support counts across both nodes.
  // The longer node's count dominates because it represents more k-mer
  // positions, giving a more reliable coverage estimate after compression.
  auto const max_size = std::max(mCounts.size(), other.mCounts.size());
  mCounts.resize(max_size, 0);

  auto const this_len = mKmer.Length();
  auto const other_len = other.Length();
  auto const total_len = this_len + other_len;

  for (usize idx = 0; idx < max_size; ++idx) {
    auto const this_count = static_cast<u64>(mCounts[idx]);
    auto const other_count =
        idx < other.mCounts.size() ? static_cast<u64>(other.mCounts[idx]) : 0ULL;
    auto const weighted_sum = (this_count * this_len) + (other_count * other_len);
    mCounts[idx] = static_cast<u32>(weighted_sum / total_len);
  }

  // Apply the same weighted average to per-role aggregates.
  // Each role counter is independently smoothed by node length, matching
  // the per-sample weighted average applied to mCounts above.
  for (usize role = 0; role < mRoleCounts.size(); ++role) {
    auto const this_rc = static_cast<u64>(mRoleCounts[role]);
    auto const other_rc = static_cast<u64>(other.mRoleCounts[role]);
    mRoleCounts[role] =
        static_cast<u32>(((this_rc * this_len) + (other_rc * other_len)) / total_len);
  }
}

auto Node::HasSelfLoop() const -> bool {
  return std::ranges::any_of(mEdges, [](Edge const& conn) -> bool { return conn.IsSelfLoop(); });
}

auto Node::FindEdgesInDirection(Kmer::Ordering const ord) const -> std::vector<Edge> {
  std::vector<Edge> results;
  results.reserve(mEdges.size());
  auto const expected_src_sign = mKmer.SignFor(ord);
  std::ranges::copy_if(mEdges, std::back_inserter(results),
                       [&expected_src_sign](Edge const& conn) -> bool {
                         return conn.SrcSign() == expected_src_sign;
                       });
  return results;
}

}  // namespace lancet::cbdg
