#ifndef SRC_LANCET_CALLER_ALLELE_SCORING_TYPES_H_
#define SRC_LANCET_CALLER_ALLELE_SCORING_TYPES_H_

#include "lancet/base/types.h"
#include "lancet/caller/variant_support.h"

#include "absl/types/span.h"

namespace lancet::caller {

// ============================================================================
// HapVariantBounds: result of resolving an alignment's haplotype index
// to the variant's physical coordinates on that specific haplotype.
//
// Each assembled haplotype carries the variant at a DIFFERENT position.
// The REF haplotype uses the variant's original reference position
// (mLocalRefStart0Idx). Each ALT haplotype records its own position
// via the per-ALT mLocalHapStart0Idxs map (populated during variant
// extraction from the SPOA consensus paths).
//
//   Haplotype 0 (REF): [ ... var at pos 120, len=3 (REF allele) ... ]
//   Haplotype 1 (ALT): [ ... var at pos 118, len=5 (ALT allele) ... ]
//   Haplotype 2 (ALT): [ ... var at pos 119, len=4 (ALT allele) ... ]
// ============================================================================
struct HapVariantBounds {
  // ── 4B Align ──────────────────────────────────────────────────────────
  i32 mVarStart = 0;  // variant start position on this haplotype
  i32 mVarLen = 0;    // variant allele length on this haplotype

  // ── 1B Align ──────────────────────────────────────────────────────────
  AlleleIndex mAllele = REF_ALLELE_IDX;  // which allele this haplotype maps to
};

// ============================================================================
// ReadAlnContext: per-read invariant data constructed once per read and
// threaded through all variant scoring calls.
//
// Bundles the three query-read parameters that are constant across all
// variants within a single call to AssignReadToAlleles.
// ============================================================================
struct ReadAlnContext {
  // ── 16B Align ─────────────────────────────────────────────────────────
  absl::Span<u8 const> mSeqEncoded;  // numeric-encoded read sequence
  absl::Span<u8 const> mBaseQuals;   // per-base Phred quality scores

  // ── 8B Align ──────────────────────────────────────────────────────────
  usize mReadLength = 0;  // original read length (for folded position)
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_ALLELE_SCORING_TYPES_H_
