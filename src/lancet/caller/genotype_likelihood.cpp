#include "lancet/caller/genotype_likelihood.h"

#include "lancet/base/types.h"
#include "lancet/hts/phred_quality.h"

#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <algorithm>
#include <limits>
#include <numbers>
#include <vector>

#include <cmath>

namespace lancet::caller {

namespace {

// ============================================================================
// DM Model Constants
// ============================================================================
// These are hardcoded defaults tuned for Illumina short-read sequencing.
// Both values interact: raising ρ widens the overdispersion tail (absorbs
// more correlated noise but reduces sensitivity to true low-VAF germline
// variants). Raising ε distributes more background mass (better numeric
// stability at the cost of slightly reduced genotype discrimination).

/// Background error rate: fraction of reads expected to be mismapped or
/// mis-called even under the correct genotype hypothesis.
///
/// Default: 0.005 (0.5%), matching typical Illumina error profiles.
///
/// When to change:
///   - ONT/PacBio CLR: raise to 0.01–0.03 (higher per-base error rate)
///   - HiFi/CCS:       keep at 0.005 (error correction makes it Illumina-like)
///   - Ultra-clean PCR-free libraries: could lower to 0.002–0.003
///
/// Effect of increasing: more probability mass on "background," reducing
/// PL separation between genotypes. The model becomes more conservative
/// (fewer high-confidence calls) but more robust to systematic noise.
/// Effect of decreasing: sharper genotype discrimination but more
/// vulnerable to correlated sequencing artifacts at ultra-high depth.
constexpr f64 DM_BACKGROUND_ERROR = 0.005;

/// Overdispersion parameter ρ. Controls the heavy-tail width of the DM
/// distribution. Higher ρ = more variance absorbed = PLs plateau earlier.
///
/// Precision M = (1 − ρ) / ρ.
///   ρ = 0.01  → M = 99   (default, moderate overdispersion)
///   ρ = 0.001 → M = 999  (nearly binomial, PLs scale ~linearly with depth)
///   ρ = 0.1   → M = 9    (aggressive, PLs plateau very early)
///
/// Default: 0.01 — empirically balances depth-robustness (2500x WGS) against
/// sensitivity to genuine germline heterozygotes at standard coverage (30×–60×).
///
/// When to change:
///   - Ultra-deep targeted panels (2000x+): consider 0.02–0.05 if PLs still
///     explode, indicating higher-than-expected error correlation
///   - Low-coverage WGS (10×–15×): consider 0.005 to preserve every quantum
///     of evidence (overdispersion barely matters at low N)
///   - ONT duplex: consider 0.02 (systematic strand-phase errors)
///
/// Effect of increasing: PLs plateau at lower depth. A true het at 30×
/// may produce PL[0/1]=0, PL[0/0]=45 instead of 60 — lower GQ but more
/// robust against correlated artifacts.
/// Effect of decreasing: PLs scale closer to linearly with depth.
/// A true het at 2000× could produce PL[0/0] > 50,000 — technically
/// correct but numerically extreme and potentially misleading for
/// downstream ML models that expect bounded features.
constexpr f64 DM_OVERDISPERSION = 0.01;

/// Minimum alpha floor to prevent lgamma(0) singularity.
/// lgamma(x) diverges as x → 0. This floor is small enough to never
/// affect genotype discrimination (α_i ≈ 0 means zero expected counts
/// for that allele under this genotype, which is the correct behavior).
constexpr f64 DM_ALPHA_FLOOR = 1e-6;

// ============================================================================
// LogDirichletMultinomial: DM log-likelihood computation.
//
// In plain terms: given an expected allele mix (dm_alphas, derived from the
// genotype hypothesis) and the actual read counts, this function returns a
// single number measuring "how likely is this data under this hypothesis?"
// More negative = worse fit. The log-gamma terms handle the combinatorial
// bookkeeping of how many ways the observed counts could arise.
//
//   ln P(c | α) = lnΓ(Σα_i) − lnΓ(N + Σα_i) + Σ[lnΓ(c_i + α_i) − lnΓ(α_i)]
//
// Uses std::lgamma for numerical stability. All values in natural-log space;
// conversion to Phred (log10) happens only at the final PL scaling step.
// ============================================================================
auto LogDirichletMultinomial(absl::Span<int const> allele_counts, std::vector<f64> const& dm_alphas)
    -> f64 {
  f64 log_prob = 0.0;
  f64 alpha_sum = 0.0;
  f64 count_alpha_sum = 0.0;

  for (usize index = 0; index < allele_counts.size(); ++index) {
    auto const count_val = static_cast<f64>(allele_counts[index]);
    f64 const alpha_val = dm_alphas[index];
    log_prob += std::lgamma(count_val + alpha_val) - std::lgamma(alpha_val);
    alpha_sum += alpha_val;
    count_alpha_sum += count_val + alpha_val;
  }

  log_prob += std::lgamma(alpha_sum) - std::lgamma(count_alpha_sum);
  return log_prob;
}

/// Normalize raw natural-log likelihoods to Phred-scaled PLs.
/// Best genotype gets PL=0; others are scaled relative to it.
/// PL = −10 × (LL − best_LL) / ln(10)
auto NormalizeToPLs(std::vector<f64> const& log_likelihoods) -> absl::InlinedVector<u32, 6> {
  static constexpr f64 PL_CAP = static_cast<f64>(std::numeric_limits<u32>::max()) / 2.0;

  f64 const best_log_lk = *std::ranges::max_element(log_likelihoods);

  f64 const ln_ten = std::numbers::ln10;
  absl::InlinedVector<u32, 6> phred_likelihoods(log_likelihoods.size());
  for (usize index = 0; index < log_likelihoods.size(); ++index) {
    f64 const raw_pl = -10.0 * (log_likelihoods[index] - best_log_lk) / ln_ten;
    phred_likelihoods[index] = static_cast<u32>(std::round(std::min(raw_pl, PL_CAP)));
  }

  return phred_likelihoods;
}

/// Compute the log10-probability of a single read under the mixture model.
///
/// In plain terms: given a read that was called as allele 's' with a certain
/// base quality, this function asks "how likely is this read under the current
/// allele frequency mix?" A high-quality read matching a common allele is very
/// likely (high probability); a high-quality read matching a rare allele is
/// less likely; a low-quality read is always somewhat ambiguous.
///
/// Sums over all possible true allele origins using the Law of Total Probability.
///
///   P(read called as s | f) = Σ_t f[t] × P(s | true_origin = t)
///   P(s|t) = (1−ε) if s==t, else ε/(K−1)
auto ReadMixtureProbLog10(int const called_as, u8 const base_qual, absl::Span<f64 const> frac_vec,
                          int const num_alleles) -> f64 {
  f64 const error_prob = hts::PhredToErrorProb(base_qual);
  f64 const mismatch_prob = error_prob / std::max(1, num_alleles - 1);
  f64 const match_bonus = (1.0 - error_prob) - mismatch_prob;

  // Start with the uniform mismatch contribution from all alleles:
  //   Σ_t f[t] × mismatch_prob = mismatch_prob  (since Σf = 1)
  // Then add the extra mass from the matching allele.
  f64 const prob_read = mismatch_prob + (frac_vec[called_as] * match_bonus);
  return std::log10(std::max(1e-15, prob_read));
}

/// Accumulate log10-likelihoods from a vector of base qualities for one allele.
auto AccumulateStrandLogLk(absl::Span<u8 const> base_quals, int const called_as,
                           absl::Span<f64 const> frac_vec, int const num_alleles) -> f64 {
  f64 strand_log_lk = 0.0;
  for (auto const base_qual : base_quals) {
    strand_log_lk += ReadMixtureProbLog10(called_as, base_qual, frac_vec, num_alleles);
  }
  return strand_log_lk;
}

/// Compute the full-pileup log10-likelihood under a K-allele mixture model.
/// Sums per-read contributions (fwd + rev strands) across all alleles.
auto PileupLogLikelihood(absl::Span<AlleleBaseQuals const> allele_quals,
                         absl::Span<f64 const> frac_vec, int const num_alleles) -> f64 {
  f64 log_lk = 0.0;
  for (int called_as = 0; called_as < num_alleles; ++called_as) {
    auto const& per_allele = allele_quals[called_as];
    log_lk += AccumulateStrandLogLk(per_allele.mFwdBaseQuals, called_as, frac_vec, num_alleles);
    log_lk += AccumulateStrandLogLk(per_allele.mRevBaseQuals, called_as, frac_vec, num_alleles);
  }
  return log_lk;
}

/// Build the null-hypothesis fraction vector: redistribute
/// the target ALT's mass proportionally among remaining alleles.
auto BuildNullFractions(absl::Span<f64 const> frac_mle, int const target_alt, int const num_alleles)
    -> std::vector<f64> {
  std::vector<f64> frac_null(frac_mle.begin(), frac_mle.end());
  f64 const null_mass = frac_null[target_alt];
  frac_null[target_alt] = 0.0;

  f64 const remaining_sum = 1.0 - null_mass;
  if (remaining_sum <= 0.0) {
    frac_null[0] = 1.0;
    return frac_null;
  }

  for (int index = 0; index < num_alleles; ++index) {
    frac_null[index] /= remaining_sum;
  }

  return frac_null;
}

}  // namespace

// ============================================================================
// ComputeGenotypePLs: Phred-scaled likelihoods from allele counts.
//
// For K alleles, evaluates K*(K+1)/2 diploid genotype hypotheses.
// Each hypothesis predicts an allele fraction vector μ, and the DM model
// scores how well the observed counts match that prediction.
// ============================================================================
auto ComputeGenotypePLs(absl::Span<int const> allele_counts) -> absl::InlinedVector<u32, 6> {
  auto const num_alleles = static_cast<int>(allele_counts.size());
  if (num_alleles == 0) return {};

  auto const num_genotypes = num_alleles * (num_alleles + 1) / 2;

  // DM precision: M = (1 − ρ) / ρ.
  // Higher M → closer to multinomial (less overdispersion).
  // ρ=0.01 → M=99   ρ=0.001 → M=999   ρ=0.1 → M=9
  f64 const precision = (1.0 - DM_OVERDISPERSION) / DM_OVERDISPERSION;

  std::vector<f64> genotype_log_lks(num_genotypes, 0.0);
  int genotype_idx = 0;

  // VCF-standard unphased genotype ordering:
  //   (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...
  for (int allele_b = 0; allele_b < num_alleles; ++allele_b) {
    for (int allele_a = 0; allele_a <= allele_b; ++allele_a) {
      // Build expected allele fraction vector for this genotype.
      // Base: distribute ε uniformly across all alleles as background noise.
      std::vector<f64> expected_mu(num_alleles, DM_BACKGROUND_ERROR / num_alleles);
      f64 const main_mass = 1.0 - DM_BACKGROUND_ERROR;

      if (allele_a == allele_b) {
        // Homozygous: all main mass on one allele
        expected_mu[allele_a] += main_mass;
      } else {
        // Heterozygous: split main mass equally between two alleles
        expected_mu[allele_a] += main_mass / 2.0;
        expected_mu[allele_b] += main_mass / 2.0;
      }

      // Map expected fractions to Dirichlet concentration parameters: α_i = M × μ_i
      std::vector<f64> dm_alphas(num_alleles);
      for (int kidx = 0; kidx < num_alleles; ++kidx) {
        dm_alphas[kidx] = std::max(DM_ALPHA_FLOOR, expected_mu[kidx] * precision);
      }

      genotype_log_lks[genotype_idx] = LogDirichletMultinomial(allele_counts, dm_alphas);
      ++genotype_idx;
    }
  }

  return NormalizeToPLs(genotype_log_lks);
}

// ============================================================================
// ComputeGenotypeQuality: GQ from PLs.
//
// Standard GATK convention: GQ is the second-smallest PL value.
// After normalization, the smallest PL is always 0 (the called genotype).
// GQ = second-smallest PL, capped at 99.
// ============================================================================
auto ComputeGenotypeQuality(absl::Span<u32 const> phred_likelihoods) -> u32 {
  if (phred_likelihoods.size() < 2) return 0;

  // Find minimum and second-minimum PL values
  u32 min1 = std::numeric_limits<u32>::max();
  u32 min2 = std::numeric_limits<u32>::max();
  for (u32 const val : phred_likelihoods) {
    if (val < min1) {
      min2 = min1;
      min1 = val;
    } else if (val < min2) {
      min2 = val;
    }
  }

  // After normalization, min1 should be 0. GQ = min2 - min1 = min2.
  static constexpr u32 MAX_GQ = 99;
  return std::min(min2 - min1, MAX_GQ);
}

// ============================================================================
// ComputeContinuousMixtureLods: per-ALT CMLOD scores.
//
// In plain terms: CMLOD measures "how much evidence is there that this
// variant is real?" For each ALT allele, it compares two scenarios:
//   1. The variant exists at its observed frequency (e.g., 15% VAF)
//   2. The variant doesn't exist — those ALT reads are just errors
// The score is the log-odds between these scenarios. Higher = more
// evidence the variant is real. Each read's contribution is weighted by
// its base quality, so 10 high-quality reads matter far more than 10
// low-quality reads at the same allele frequency.
//
// For mosaic/somatic calling, rigid diploid genotype states are inappropriate.
// A subclonal variant at 2% VAF is a continuous frequency, not a discrete
// diploid state (0/0, 0/1, or 1/1).
//
// CMLOD integrates exact per-read base qualities into a K-dimensional mixture
// model, breaking the count-based identifiability paradox:
//
//   P(read called as s | f) = Σ_t f[t] × P(s | true_origin = t)
//   P(s|t) = (1−ε) if s==t, else ε/(K−1)
//
//   LL(f) = Σ_reads log10( P(read | f) )
//   CMLOD[alt] = max(0, LL(f_MLE) − LL(f_null))
//
// Complexity: O(N × K) per ALT allele, where N = total read count.
// ============================================================================
auto ComputeContinuousMixtureLods(absl::Span<AlleleBaseQuals const> allele_quals,
                                  absl::Span<usize const> allele_coverages) -> std::vector<f64> {
  auto const num_alleles = static_cast<int>(allele_quals.size());
  std::vector<f64> lod_scores(num_alleles, 0.0);
  if (num_alleles < 2) return lod_scores;

  // Compute total depth
  usize total_depth = 0;
  for (auto const coverage : allele_coverages) {
    total_depth += coverage;
  }
  if (total_depth == 0) return lod_scores;

  auto const total_depth_f64 = static_cast<f64>(total_depth);

  // Empirical allele fractions (MLE): f_i = count_i / total_depth
  std::vector<f64> frac_mle(num_alleles, 0.0);
  for (int index = 0; index < num_alleles; ++index) {
    frac_mle[index] = static_cast<f64>(allele_coverages[index]) / total_depth_f64;
  }

  f64 const log_lk_mle =
      PileupLogLikelihood(allele_quals, absl::MakeConstSpan(frac_mle), num_alleles);

  // Compute CMLOD for each ALT allele independently
  for (int target_alt = 1; target_alt < num_alleles; ++target_alt) {
    if (allele_coverages[target_alt] == 0) continue;

    auto const frac_null =
        BuildNullFractions(absl::MakeConstSpan(frac_mle), target_alt, num_alleles);

    f64 const log_lk_null =
        PileupLogLikelihood(allele_quals, absl::MakeConstSpan(frac_null), num_alleles);

    lod_scores[target_alt] = std::max(0.0, log_lk_mle - log_lk_null);
  }

  return lod_scores;
}

}  // namespace lancet::caller
