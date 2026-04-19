#include "lancet/base/repeat.h"

#include "lancet/base/assert.h"
#include "lancet/base/types.h"

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"

#include <algorithm>
#include <string_view>
#include <utility>

namespace lancet::base {

// ============================================================================
// HammingDist — auto-vectorized byte-level mismatch count
//
// Counts positions where two equal-length byte strings differ.  The inner
// loop accumulates into u8 in batches of 255 (max before u8 overflow).
//
// Why u8?  A comparison result is one byte (0 or 1), but a usize accumulator
// is 8 bytes.  When the accumulator is wider than the data, the compiler must
// insert extra instructions to widen each 1-byte result to 8 bytes before
// adding — a 4.7× instruction bloat on GCC 15.2 with AVX2 (168 vs 36 SIMD
// instructions).  Capping the inner loop at 255 iterations proves to the
// compiler's value-range analysis that no byte lane can overflow, so it
// safely drops the widening and sums byte lanes directly in hardware.
//
// Pure C++ — no intrinsics.  Works on x86 AVX2/AVX-512 and ARM64 NEON.
// For Lancet2 k-mer lengths (13–127 bp), the outer loop executes once.
// ============================================================================
auto HammingDist(std::string_view first, std::string_view second) -> usize {
  LANCET_ASSERT(first.length() == second.length())

  usize result = 0;
  auto const length = first.length();
  usize idx = 0;

  // 255 is the largest count a u8 can hold without overflow (each iteration
  // adds 0 or 1).  One scalar widen per 255 bytes is negligible overhead.
  while (idx < length) {
    u8 batch_sum = 0;
    auto const batch_end = std::min(idx + 255, length);
    for (; idx < batch_end; ++idx) {
      batch_sum += static_cast<u8>(first[idx] != second[idx]);
    }
    result += batch_sum;
  }

  return result;
}

// ============================================================================
// HasRepeat — repeat detector with fast path for exact matches
//
// For exact repeats (max_mismatches=0), uses an O(n) hash-set lookup — the
// set detects duplicate string_views by content, short-circuiting as soon as
// a collision is found.  For approximate repeats, falls back to upper-triangle
// pairwise Hamming distance checks: O(n(n-1)/2) comparisons, returning true
// as soon as any pair differs in at most `max_mismatches` positions.
// ============================================================================
auto HasRepeat(absl::Span<std::string_view const> kmers, usize const max_mismatches) -> bool {
  // Exact repeat: O(n) hash-set duplicate detection
  if (max_mismatches == 0) {
    absl::flat_hash_set<std::string_view> seen;
    seen.reserve(kmers.size());
    for (auto const kmer : kmers) {
      if (auto [iter, inserted] = seen.insert(kmer); !inserted) return true;
    }
    return false;
  }

  // Approximate repeat: O(n²) pairwise scan with short-circuit
  auto const num_kmers = kmers.size();
  for (usize i = 0; i < num_kmers; ++i) {
    for (usize j = i + 1; j < num_kmers; ++j) {
      if (HammingDist(kmers[i], kmers[j]) <= max_mismatches) return true;
    }
  }
  return false;
}

auto HasExactRepeat(absl::Span<std::string_view const> kmers) -> bool {
  return HasRepeat(kmers, 0);
}

}  // namespace lancet::base
