#ifndef SRC_LANCET_CALLER_GENOTYPE_LIKELIHOOD_H_
#define SRC_LANCET_CALLER_GENOTYPE_LIKELIHOOD_H_

#include "lancet/base/types.h"

#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <vector>

namespace lancet::caller {

// ============================================================================
// Genotype Likelihood Engine (Dirichlet-Multinomial Model)
//
// Replaces the legacy per-read binomial model with a count-based
// Dirichlet-Multinomial (DM) distribution that handles:
//
//   1. Multi-allelic sites (K alleles, not just bi-allelic)
//   2. Correlated errors at ultra-high depth (overdispersion via ρ)
//   3. Numerically stable computation via std::lgamma
//
// The DM generalizes the Beta-Binomial to K dimensions:
//   ln P(c | μ, M) = lnΓ(ΣMμ_i) − lnΓ(N + ΣMμ_i)
//                   + Σ[ lnΓ(c_i + Mμ_i) − lnΓ(Mμ_i) ]
//
// See docs/guides/variant_discovery_genotyping.md for the full derivation.
// ============================================================================

/// Compute Phred-scaled genotype likelihoods (PLs) from allele counts.
/// Uses the DM model with overdispersion to handle correlated sequencing errors.
///
/// For K alleles, returns K*(K+1)/2 PLs in VCF-standard unphased ordering:
///   (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...
/// Best genotype PL = 0; others are scaled relative to it.
[[nodiscard]] auto ComputeGenotypePLs(absl::Span<int const> allele_counts)
    -> absl::InlinedVector<u32, 6>;

/// Genotype Quality: second-smallest PL value, capped at 99.
/// Standard GATK convention: GQ = min2 - min1. After normalization, min1 = 0.
/// See: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
[[nodiscard]] auto ComputeGenotypeQuality(absl::Span<u32 const> phred_likelihoods) -> u32;

// ============================================================================
// Continuous Mixture Log-Odds (CMLOD)
//
// Per-ALT-allele LOD comparing the observed pileup likelihood under MLE
// allele fractions vs. a null hypothesis where the target ALT fraction
// is forced to zero. Integrates exact per-read base qualities.
//
//   CMLOD[alt] = max(0, LL(f_MLE) − LL(f_null))   (log10 scale)
// ============================================================================

/// Per-allele strand-split base qualities needed by the mixture LOD engine.
/// Caller must provide one entry per allele (REF at index 0, ALTs at 1..K-1).
struct AlleleBaseQuals {
  absl::Span<u8 const> mFwdBaseQuals;  // 16B
  absl::Span<u8 const> mRevBaseQuals;  // 16B
};

/// Compute per-ALT CMLOD scores. Index 0 (REF) is always 0.0.
/// Each ALT score measures evidence that the allele exists at its observed
/// frequency vs. being absent.
[[nodiscard]] auto ComputeContinuousMixtureLods(absl::Span<AlleleBaseQuals const> allele_quals,
                                                absl::Span<usize const> allele_coverages)
    -> std::vector<f64>;

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_GENOTYPE_LIKELIHOOD_H_
