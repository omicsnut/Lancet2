#ifndef SRC_LANCET_BASE_SLIDING_H_
#define SRC_LANCET_BASE_SLIDING_H_

#include "lancet/base/assert.h"
#include "lancet/base/types.h"

#include "absl/container/fixed_array.h"
#include "absl/strings/string_view.h"

#include <algorithm>
#include <ranges>
#include <string_view>

namespace lancet::base {

using SeqMers = absl::FixedArray<std::string_view>;
[[nodiscard]] inline auto SlidingView(std::string_view seq, usize const window) -> SeqMers {
  if (seq.length() < window) return absl::FixedArray<std::string_view>(0);

  auto const end_position = seq.length() - window;
  absl::FixedArray<std::string_view> result(end_position + 1);

  auto const offsets = std::views::iota(usize{0}, end_position + 1);
  std::ranges::transform(offsets, result.begin(), [&](usize const offset) {
    return absl::ClippedSubstr(seq, offset, window);
  });

  LANCET_ASSERT(std::ranges::all_of(
      result, [&](std::string_view const view) { return view.length() == window; }));

  return result;
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_SLIDING_H_
