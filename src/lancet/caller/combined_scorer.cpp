#include "lancet/caller/combined_scorer.h"

#include "lancet/base/types.h"
#include "lancet/caller/allele_scoring_types.h"
#include "lancet/caller/genotyper.h"
#include "lancet/caller/local_scorer.h"
#include "lancet/caller/scoring_constants.h"
#include "lancet/hts/cigar_utils.h"

#include "absl/types/span.h"

#include <algorithm>
#include <vector>

namespace lancet::caller {

// ============================================================================
// ComputeHaplotypeEditDistance: NM against a specific haplotype.
//
// Calculate edit distance (NM) against the specified haplotype.
// Per SAM tags spec: NM excludes clipping (soft/hard clips do NOT contribute).
// hts::ComputeEditDistance() is already spec-compliant.
// ============================================================================
auto ComputeHaplotypeEditDistance(std::vector<Mm2AlnResult> const& alns,
                                  absl::Span<u8 const> encoded_haplotype,
                                  absl::Span<u8 const> qry_seq_encoded, usize qry_read_length,
                                  usize hap_idx) -> u32 {
  for (auto const& aln : alns) {
    if (aln.mHapIdx != hap_idx || aln.mRefStart >= aln.mRefEnd) {
      continue;
    }
    auto const aln_len = static_cast<usize>(aln.mRefEnd - aln.mRefStart);
    auto const target = encoded_haplotype.subspan(static_cast<usize>(aln.mRefStart), aln_len);
    return hts::ComputeEditDistance(aln.mCigar, qry_seq_encoded, target);
  }
  // No valid alignment found — fallback to full read length as worst-case NM.
  return static_cast<u32>(qry_read_length);
}

// ============================================================================
// ScoreReadAtVariant: pure scoring of one read-haplotype alignment at a
// variant site.
//
// Caller (Genotyper::AssignReadToAlleles) is responsible for:
//   1. Resolving haplotype bounds (ExtractHapBounds)
//   2. Checking overlap (OverlapsAlignment)
//   3. Compare-and-swap ranking (CombinedScore())
//
// This function only computes scoring:
//   combined = (global_score - sc_penalty - local_raw_score) + (local_pbq * identity)
// ============================================================================
auto ScoreReadAtVariant(Mm2AlnResult const& aln, absl::Span<u8 const> encoded_haplotype,
                        ReadAlnContext const& read_ctx, HapVariantBounds const& bounds)
    -> ReadAlleleAssignment {
  auto const aln_len = static_cast<usize>(aln.mRefEnd - aln.mRefStart);
  auto const target = encoded_haplotype.subspan(static_cast<usize>(aln.mRefStart), aln_len);

  auto const local =
      ComputeLocalScore(aln.mCigar, read_ctx.mSeqEncoded, target, read_ctx.mBaseQuals,
                        aln.mRefStart, bounds.mVarStart, bounds.mVarLen, SCORING_MATRIX);

  // Subtract the soft-clip penalty from the global score.
  // This prevents supplementary alignments (where half the read is soft-clipped)
  // from incorrectly outscoring full contiguous alignments.
  f64 const sc_penalty = ComputeSoftClipPenalty(aln.mCigar);
  f64 const global_adjusted = static_cast<f64>(aln.mScore) - sc_penalty;

  ReadAlleleAssignment result;
  result.mAllele = bounds.mAllele;
  result.mGlobalScore = static_cast<i32>(global_adjusted - local.mRawScore);
  result.mLocalScore = local.mPbqScore;
  result.mLocalIdentity = local.mIdentity;
  result.mBaseQualAtVar = local.mBaseQual;

  // SPOA path ID for HSE: which haplotype did this read align to?
  result.mAssignedHaplotypeId = static_cast<u32>(aln.mHapIdx);

  // Edit distance against this read's assigned haplotype (for AHDD FORMAT
  // field). Computed for ALL reads — both REF and ALT. For REF reads,
  // aln.mHapIdx == 0 (the REF haplotype); for ALT reads, aln.mHapIdx
  // points to the winning ALT haplotype. This gives AHDD an apples-to-
  // apples baseline: mean(ALT NM vs own hap) − mean(REF NM vs own hap).
  result.mOwnHapNm = hts::ComputeEditDistance(aln.mCigar, read_ctx.mSeqEncoded, target);

  // ============================================================================
  // Folded read position
  // ============================================================================
  usize var_start_in_aln = 0;
  if (bounds.mVarStart > aln.mRefStart) {
    var_start_in_aln = static_cast<usize>(bounds.mVarStart - aln.mRefStart);
  }

  auto const qpos_at_var = hts::CigarRefPosToQueryPos(aln.mCigar, var_start_in_aln);
  auto const rel_pos = read_ctx.mReadLength > 0
                           ? static_cast<f64>(qpos_at_var) / static_cast<f64>(read_ctx.mReadLength)
                           : 0.5;
  result.mFoldedReadPos = std::min(rel_pos, 1.0 - rel_pos);

  return result;
}

}  // namespace lancet::caller
