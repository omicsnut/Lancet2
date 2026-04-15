#ifndef SRC_LANCET_CALLER_POSTERIOR_BASE_QUAL_H_
#define SRC_LANCET_CALLER_POSTERIOR_BASE_QUAL_H_

#include "lancet/base/types.h"

#include "absl/types/span.h"

namespace lancet::caller {

// ============================================================================
// RawPosteriorBaseQual (used to compute NPBQ FORMAT field)
//
// Edgar & Flyvbjerg (2014): Bayesian aggregation of per-read error probs.
//
// Given N reads supporting allele `a`, each with Phred quality Q_i:
//
//   ε_i = 10^(-Q_i / 10)
//
//   log10(P_err_all)  = Σ log10(ε_i)           // all reads are wrong
//   log10(P_ok_all)   = Σ log10(1 - ε_i)       // all reads are correct
//
//   posterior_err = P_err_all / (P_err_all + P_ok_all)
//
// In log10 space (for numerical stability via log-sum-exp):
//   log10(posterior_err) = log_err - log10(10^log_err + 10^log_ok)
//
//   raw PBQ = -10 * log10(posterior_err)
//
// Returns the raw uncapped value. The caller (BuildFormatFields) divides
// by allele_depth to produce NPBQ (per-read quality) for VCF output.
// ============================================================================

/// Compute raw posterior base quality from forward and reverse strand base quals.
/// Pure Bayesian math — no variant/allele knowledge.
[[nodiscard]] auto ComputeRawPosteriorBaseQual(absl::Span<u8 const> fwd_base_quals,
                                               absl::Span<u8 const> rev_base_quals) -> f64;

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_POSTERIOR_BASE_QUAL_H_
