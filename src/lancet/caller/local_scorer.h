#ifndef SRC_LANCET_CALLER_LOCAL_SCORER_H_
#define SRC_LANCET_CALLER_LOCAL_SCORER_H_

#include "lancet/base/types.h"
#include "lancet/hts/cigar_unit.h"

#include "absl/types/span.h"

#include <array>
#include <string_view>
#include <vector>

namespace lancet::caller {

// ============================================================================
// Result of ComputeLocalScore
// ============================================================================
struct LocalScoreResult {
  // ── 8B Align ────────────────────────────────────────────────────────────
  f64 mPbqScore = 0.0;  // 8B — PBQ-weighted DP score within variant region
  f64 mRawScore = 0.0;  // 8B — Unweighted matrix score (same paths)
  f64 mIdentity = 0.0;  // 8B — Fraction of exact matches in variant region
  // ── 1B Align ────────────────────────────────────────────────────────────
  u8 mBaseQual = 0;  // 1B — Minimum Phred base quality (weakest-link)
};

// ============================================================================
// Public API
// ============================================================================

/// ASCII DNA → numeric (0-4) for local scoring.
/// Single pass using the constexpr lookup table. O(n), no branches.
[[nodiscard]] auto EncodeSequence(std::string_view raw_seq) -> std::vector<u8>;

/// Evaluate alignment quality in a variant's physical region on the haplotype.
/// Pure scoring math — zero knowledge of variants, alleles, or minimap2.
/// See local_scorer.cpp for the full CIGAR walk algorithm.
[[nodiscard]] auto ComputeLocalScore(std::vector<hts::CigarUnit> const& qry_aln_cigar,
                                     absl::Span<u8 const> qry_seq, absl::Span<u8 const> hap_seq,
                                     absl::Span<u8 const> qry_quals, i32 aln_start_on_hap,
                                     i32 var_start_on_hap, i32 var_len_on_hap,
                                     std::array<i8, 25> const& score_matrix) -> LocalScoreResult;

/// Penalty for unaligned (soft-clipped) read tails.
/// Treats each soft-clipped base as a mismatch.
[[nodiscard]] auto ComputeSoftClipPenalty(std::vector<hts::CigarUnit> const& cigar) -> f64;

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_LOCAL_SCORER_H_
