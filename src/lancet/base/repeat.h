#ifndef SRC_LANCET_BASE_REPEAT_H_
#define SRC_LANCET_BASE_REPEAT_H_

#include "lancet/base/types.h"

#include "absl/types/span.h"

#include <string_view>

namespace lancet::base {

/// Byte-level Hamming distance: count of positions where two equal-length strings differ.
/// Designed for auto-vectorization — the compiler emits vpcmpeqb + vpsadbw on AVX2.
[[nodiscard]] auto HammingDist(std::string_view first, std::string_view second) -> usize;

/// True if any two kmers in the span differ in at most `max_mismatches` positions.
/// For exact repeats (max_mismatches=0), uses an O(n) hash-set duplicate check.
/// For approximate repeats, uses O(n²) pairwise Hamming distance with short-circuit.
/// Upper-triangle indexing ensures each pair is checked exactly once.
[[nodiscard]] auto HasRepeat(absl::Span<std::string_view const> kmers, usize max_mismatches)
    -> bool;

/// True if any two kmers in the span are identical.
/// Delegates to HasRepeat(kmers, 0), which uses an O(n) hash-set duplicate check.
[[nodiscard]] auto HasExactRepeat(absl::Span<std::string_view const> kmers) -> bool;

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_REPEAT_H_
