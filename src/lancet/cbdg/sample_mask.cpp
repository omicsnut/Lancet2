#include "lancet/cbdg/sample_mask.h"

#include <absl/container/inlined_vector.h>
#include <algorithm>
#include <bit>
#include <functional>
#include <numeric>

namespace lancet::cbdg {

static constexpr usize BITS_PER_WORD = 64;

void SampleMask::SetBit(usize const bit_index) {
  EnsureCapacity(bit_index);
  auto const word_idx = bit_index / BITS_PER_WORD;
  auto const bit_pos = bit_index % BITS_PER_WORD;
  mWords[word_idx] |= (1ULL << bit_pos);
}

auto SampleMask::TestBit(usize const bit_index) const -> bool {
  auto const word_idx = bit_index / BITS_PER_WORD;
  if (word_idx >= mWords.size()) return false;

  auto const bit_pos = bit_index % BITS_PER_WORD;
  return (mWords[word_idx] & (1ULL << bit_pos)) != 0;
}

void SampleMask::Merge(SampleMask const& other) {
  if (other.mWords.size() > mWords.size()) {
    mWords.resize(other.mWords.size(), 0ULL);
  }

  // OR each source word into the corresponding target word.
  // std::transform with std::bit_or produces the set-union of both masks.
  std::transform(other.mWords.begin(), other.mWords.end(), mWords.begin(), mWords.begin(),
                 std::bit_or<u64>{});
}

auto SampleMask::PopCount() const -> usize {
  // Sum popcount across all words. std::transform_reduce fuses the
  // per-element transform (std::popcount) with the accumulation (addition)
  // into a single pass — more expressive than a manual loop.
  return std::transform_reduce(
      mWords.begin(), mWords.end(), usize{0}, std::plus<>{},
      [](u64 const word) -> usize { return static_cast<usize>(std::popcount(word)); });
}

void SampleMask::EnsureCapacity(usize const bit_index) {
  auto const needed = (bit_index / BITS_PER_WORD) + 1;
  if (needed > mWords.size()) {
    mWords.resize(needed, 0ULL);
  }
}

}  // namespace lancet::cbdg
