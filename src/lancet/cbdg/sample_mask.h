#ifndef SRC_LANCET_CBDG_SAMPLE_MASK_H_
#define SRC_LANCET_CBDG_SAMPLE_MASK_H_

#include "lancet/base/types.h"

#include "absl/container/inlined_vector.h"

namespace lancet::cbdg {

// ============================================================================
// SampleMask — dynamic per-sample bitmask for graph k-mer coloring.
//
// Generalizes beyond the fixed 3-bit Label::Tag system (REFERENCE | CTRL | CASE)
// with a scalable design that supports N samples without code changes.
//
// Bit layout:
//   Bit 0:             REFERENCE (always reserved)
//   Bit 1..N:          sample indices (0-based → bit position = index + 1)
//
// Storage:
//   ≤63 samples:  single inline u64 (no heap allocation via InlinedVector)
//   >63 samples:  spills to additional heap-backed u64 words
//
// Each u64 word covers 64 bit positions:
//   Word 0 → bits  0–63
//   Word 1 → bits 64–127
//   ...
//
// Intended use: during graph construction, each k-mer node carries a
// SampleMask that tracks which samples contributed reads containing that
// k-mer. The mask is merged (OR'd) when the same k-mer is seen from
// multiple reads or samples, producing a set-union of sample identities.
//
// Thread safety: none. Graph construction is single-threaded per window.
// ============================================================================
class SampleMask {
 public:
  SampleMask() = default;

  /// Set the bit at `bit_index`. Grows storage if needed.
  /// Bit 0 = REFERENCE. Bit (sample_index + 1) = sample identity.
  void SetBit(usize bit_index);

  /// Test whether the bit at `bit_index` is set.
  /// Returns false for out-of-range indices (no allocation on read).
  [[nodiscard]] auto TestBit(usize bit_index) const -> bool;

  /// Merge another mask into this one (bitwise OR across all words).
  /// After merge, this mask contains the set-union of both masks.
  void Merge(SampleMask const& other);

  /// Count total set bits across all words (includes reference bit if set).
  [[nodiscard]] auto PopCount() const -> usize;

  auto operator==(SampleMask const& rhs) const -> bool = default;
  auto operator!=(SampleMask const& rhs) const -> bool = default;

 private:
  /// Each u64 word holds 64 bits. Word 0 covers bits 0–63, etc.
  /// InlinedVector<u64, 1> stores the first word inline (no heap)
  /// and spills to heap only when >63 samples are used.
  // ── 8B Align ────────────────────────────────────────────────────────────
  absl::InlinedVector<u64, 1> mWords{0ULL};

  /// Ensure mWords has enough entries to hold `bit_index`.
  void EnsureCapacity(usize bit_index);
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_SAMPLE_MASK_H_
