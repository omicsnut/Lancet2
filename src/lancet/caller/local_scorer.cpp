#include "lancet/caller/local_scorer.h"

#include "lancet/base/types.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/phred_quality.h"

#include "absl/types/span.h"

#include <algorithm>
#include <array>
#include <string_view>
#include <vector>

namespace lancet::caller {

namespace {

/// Convert a Phred quality score to confidence weight: 1 - 10^(-Q/10).
/// Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99, Q=30 → 0.999, Q=40 → 0.9999
inline auto PhredToConfidence(u8 const qual) -> f64 {
  return 1.0 - hts::PhredToErrorProb(qual);
}

// ============================================================================
// RegionAccumulator: Local variant scoring abstraction state.
//
// Encapsulates scoring so the CIGAR walk loop stays clean. Evaluates alignment
// quality specifically within a variant's physical array boundaries via absolute
// haplotypic mapping rather than relative offsets.
// ============================================================================
struct RegionAccumulator {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  absl::Span<u8 const> const mQuery;         // 16B
  absl::Span<u8 const> const mTarget;        // 16B
  absl::Span<u8 const> const mBaseQuals;     // 16B
  std::array<i8, 25> const& mScoringMatrix;  // 8B
  f64 mPbqScore = 0.0;                       // 8B
  f64 mRawScore = 0.0;                       // 8B
  usize mMatches = 0;                        // 8B
  usize mAligned = 0;                        // 8B
  i32 const mAlnRefStart;                    // 4B
  i32 const mVarStartHap;                    // 4B
  i32 const mVarEndHap;                      // 4B
  u8 mMinBq = 255;                           // 1B
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

  // Checks if the current alignment position overlaps the variant being scored.
  //
  //   Haplotype Array :  [ A  B  C  D  E  F  G  H  I ]
  //   Variant Region  :           [var_start ... var_end)
  //   Alignment       :     [aln_start ... tpos_rel ... ]
  //
  // We translate `tpos_rel` (which starts at 0 for the alignment) into an
  // absolute position (`abs_pos`) on the haplotype to see if it falls inside
  // the [mVarStartHap, mVarEndHap) window.
  [[nodiscard]] auto InRegion(i32 tpos_rel) const -> bool {
    i32 const abs_pos = mAlnRefStart + tpos_rel;
    return abs_pos >= mVarStartHap && abs_pos < mVarEndHap;
  }

  // Scores a single query-target base pair and adds it to the running totals.
  //   1. Raw Score: The penalty from the substitution matrix (e.g. mismatch = -4).
  //   2. PBQ Score: The raw penalty scaled by the base quality confidence.
  //                 A low quality mismatch is penalized less than a high quality one.
  //   3. Identity : Tracks pure exact matches to compute the exact-match fraction.
  void ScoreAlignedPair(i32 tpos_rel, usize qpos) {
    ++mAligned;
    if (qpos >= mQuery.size() || static_cast<usize>(tpos_rel) >= mTarget.size()) {
      return;
    }

    auto const raw = mScoringMatrix[(mTarget[tpos_rel] * 5) + mQuery[qpos]];
    mRawScore += static_cast<f64>(raw);
    f64 const weight = (qpos < mBaseQuals.size()) ? PhredToConfidence(mBaseQuals[qpos]) : 1.0;
    mPbqScore += static_cast<f64>(raw) * weight;
    mMatches += static_cast<usize>(mQuery[qpos] == mTarget[tpos_rel]);
  }

  // Record the weakest base quality in the variant region (weakest-link confidence).
  void TrackBaseQual(usize qpos) {
    if (qpos < mBaseQuals.size()) {
      mMinBq = std::min(mMinBq, mBaseQuals[qpos]);
    }
  }

  // Deletions do not have their own quality scores since the bases are missing.
  // Instead, we estimate the deletion's confidence by checking the quality of
  // the adjacent bases immediately surrounding the deletion.
  void TrackDeletionBounds(usize qpos) {
    if (qpos > 0 && qpos - 1 < mBaseQuals.size()) {
      mMinBq = std::min(mMinBq, mBaseQuals[qpos - 1]);
    }
    if (qpos < mBaseQuals.size()) {
      mMinBq = std::min(mMinBq, mBaseQuals[qpos]);
    }
  }

  // Finalize stats. min_bq defaults to 0 if we lacked base quality evidence entirely.
  [[nodiscard]] auto ToResult() const -> LocalScoreResult {
    return {
        .mPbqScore = mPbqScore,
        .mRawScore = mRawScore,
        .mIdentity = mAligned > 0 ? static_cast<f64>(mMatches) / static_cast<f64>(mAligned) : 0.0,
        .mBaseQual = (mMinBq == 255) ? u8{0} : mMinBq,
    };
  }
};

}  // namespace

// ============================================================================
// EncodeSequence: ASCII DNA → numeric (0-4) for local scoring.
// Single pass using the constexpr lookup table. O(n), no branches.
// ============================================================================
auto EncodeSequence(std::string_view const raw_seq) -> std::vector<u8> {
  std::vector<u8> encoded(raw_seq.size());
  std::ranges::transform(raw_seq, encoded.begin(),
                         [](char base) -> u8 { return ENCODE_TABLE[static_cast<u8>(base)]; });
  return encoded;
}

// ============================================================================
// ComputeLocalScore: evaluate alignment quality in a variant's region.
//
// Given a read→haplotype CIGAR alignment, this function extracts specific metrics
// for the sub-region of the haplotype that contains the variant:
//
//   1. mPbqScore:  PBQ-weighted DP score. Each position's substitution matrix
//                  contribution is scaled by (1 - ε) where ε = 10^(-PBQ/10).
//                  This down-weights low-confidence bases and up-weights
//                  high-confidence ones, analogous to GATK's PairHMM which
//                  bakes PBQ directly into the per-read log-likelihood.
//
//   2. mIdentity:  fraction of aligned bases that are exact matches.
//
//   4. mBaseQual:  Minimum Phred base quality (weakest-link) or flanking boundary.
//
// CRITICAL: tpos coordinates in the CIGAR are relative to the alignment start
// (ref_start from mm_map), NOT position 0 of the haplotype. The caller must
// adjust var_start by subtracting ref_start before calling this function.
//
//   Haplotype:   |----[var_start..........var_end)------|
//                      ^ref_start
//   CIGAR tpos:  0  1  2 ...
//                      ^var_start_in_aln = var_start - ref_start
// ============================================================================
// NOLINTNEXTLINE(readability-function-size)
auto ComputeLocalScore(std::vector<hts::CigarUnit> const& qry_aln_cigar,
                       absl::Span<u8 const> qry_seq, absl::Span<u8 const> hap_seq,
                       absl::Span<u8 const> qry_quals, i32 aln_start_on_hap, i32 var_start_on_hap,
                       i32 var_len_on_hap, std::array<i8, 25> const& score_matrix)
    -> LocalScoreResult {
  if (qry_aln_cigar.empty() || var_len_on_hap == 0) {
    return {};
  }

  RegionAccumulator acc{
      .mQuery = qry_seq,
      .mTarget = hap_seq,
      .mBaseQuals = qry_quals,
      .mScoringMatrix = score_matrix,
      .mAlnRefStart = aln_start_on_hap,
      .mVarStartHap = var_start_on_hap,
      .mVarEndHap = var_start_on_hap + var_len_on_hap,
  };
  i32 tpos = 0;
  usize qpos = 0;

  for (auto const& unit : qry_aln_cigar) {
    if (aln_start_on_hap + tpos >= acc.mVarEndHap && unit.ConsumesReference()) {
      break;
    }

    auto const cigar_op = unit.Operation();
    u32 const len = unit.Length();

    switch (cigar_op) {
      // ── Substitution Mappings ──
      // Consumes both query and target bases. Maps bases through the substitution
      // matrix to compute alignment identity and mismatch penalties within the
      // variant boundary (not the flanking context).
      case hts::CigarOp::ALIGNMENT_MATCH:
      case hts::CigarOp::SEQUENCE_MATCH:
      case hts::CigarOp::SEQUENCE_MISMATCH: {
        for (u32 i = 0; i < len; ++i, ++tpos, ++qpos) {
          if (!acc.InRegion(tpos)) {
            continue;
          }
          acc.ScoreAlignedPair(tpos, qpos);
          acc.TrackBaseQual(qpos);
        }
        break;
      }

      // ── Insertion Geometry ──
      // Consumes query bases only. Penalizes alignment quality via gap extension
      // weights because inserted bases have no reference counterpart.
      case hts::CigarOp::INSERTION: {
        // Inserted bases don't advance tpos, but if we're inside the variant
        // region they count as aligned content and contribute PBQ.
        bool const in_region = acc.InRegion(tpos);
        for (u32 i = 0; i < len; ++i, ++qpos) {
          if (!in_region) {
            continue;
          }
          ++acc.mAligned;
          acc.TrackBaseQual(qpos);
          acc.mRawScore += static_cast<f64>(SCORING_GAP_EXTEND);
          acc.mPbqScore += static_cast<f64>(SCORING_GAP_EXTEND);
        }
        break;
      }

      // ── Deletion Geometry ──
      // Consumes target bases only. Penalizes via gap extensions. Deletions have no
      // query bases, so the quality score borrows confidence from the flanking neighbors.
      case hts::CigarOp::DELETION: {
        for (u32 i = 0; i < len; ++i, ++tpos) {
          if (acc.InRegion(tpos)) {
            ++acc.mAligned;
            acc.mRawScore += static_cast<f64>(SCORING_GAP_EXTEND);
            acc.mPbqScore += static_cast<f64>(SCORING_GAP_EXTEND);
          }
        }
        acc.TrackDeletionBounds(qpos);
        break;
      }

      // ── Unmapped Sequences ──
      // Advances position counters without scoring.
      // Soft-clip penalty is handled in ComputeSoftClipPenalty.
      case hts::CigarOp::SOFT_CLIP: {
        qpos += len;
        break;
      }
      case hts::CigarOp::REFERENCE_SKIP: {
        tpos += static_cast<i32>(len);
        break;
      }
      default:
        break;
    }
  }

  return acc.ToResult();
}

// ============================================================================
// ComputeSoftClipPenalty: penalty for unaligned read tails.
//
//   Read:  [Soft Clip] ====== Aligned Sequence ====== [Soft Clip]
//          5' End                                     3' End
//
// We treat every soft-clipped base as a mismatch to penalize
// reads that only partially align to the haplotype.
// ============================================================================
auto ComputeSoftClipPenalty(std::vector<hts::CigarUnit> const& cigar) -> f64 {
  if (cigar.empty()) {
    return 0.0;
  }

  auto const& first = cigar.front();
  auto const& last = cigar.back();

  auto const is_5p_clipped = first.Operation() == hts::CigarOp::SOFT_CLIP;
  auto const is_3p_clipped = cigar.size() > 1 && last.Operation() == hts::CigarOp::SOFT_CLIP;

  i32 const unaligned_5p = is_5p_clipped ? static_cast<i32>(first.Length()) : 0;
  i32 const unaligned_3p = is_3p_clipped ? static_cast<i32>(last.Length()) : 0;

  return static_cast<f64>(unaligned_5p + unaligned_3p) * SCORING_MISMATCH;
}

}  // namespace lancet::caller
