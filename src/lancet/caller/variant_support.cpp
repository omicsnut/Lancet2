#include "lancet/caller/variant_support.h"

#include "lancet/base/mann_whitney.h"
#include "lancet/base/types.h"
#include "lancet/hts/phred_quality.h"

#include <absl/types/span.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>

namespace lancet::caller {

// ============================================================================
// AddEvidence
// ============================================================================
void VariantSupport::AddEvidence(ReadEvidence const& evidence) {
  EnsureAlleleSlot(evidence.mAllele);
  auto& data = mAlleleData[evidence.mAllele];

  // Deduplicate by read name hash: only first-seen strand counts.
  auto const [iter, inserted] = data.mNameHashes.try_emplace(evidence.mRnameHash, evidence.mStrand);
  if (!inserted) return;

  // Each read contributes exactly one representative base quality value.
  // For indels, this is already collapsed to min(PBQ) across the variant region
  // by the genotyper (see ReadAlleleAssignment::mBaseQualAtVar comment).
  if (evidence.mStrand == Strand::FWD) {
    data.mFwdBaseQuals.push_back(evidence.mBaseQual);
  } else {
    data.mRevBaseQuals.push_back(evidence.mBaseQual);
  }

  data.mMapQuals.push_back(evidence.mMapQual);
  data.mAlnScores.push_back(evidence.mAlnScore);

  // Track soft-clip status from original alignment (for SCA FORMAT tag)
  if (evidence.mIsSoftClipped) ++data.mSoftClipCount;

  // Track insert sizes from properly-paired reads (for FLD FORMAT tag).
  // Signed insert size preserved: negative = ALT shorter (cfDNA biology).
  if (evidence.mIsProperPair && evidence.mInsertSize != 0) {
    data.mProperPairIsizes.push_back(static_cast<f64>(evidence.mInsertSize));
  }

  // Track folded read position (for RPCD FORMAT tag)
  data.mFoldedReadPositions.push_back(evidence.mFoldedReadPos);

  // Track edit distance to REF haplotype (for ASMD FORMAT tag)
  data.mRefNmValues.push_back(static_cast<f64>(evidence.mRefNm));

  // Track fragment start position (for FSSE: PCR optical duplicate detection)
  data.mAlignmentStarts.push_back(evidence.mAlignmentStart);

  // Track edit distance against assigned ALT haplotype (for AHDD)
  data.mAltNmValues.push_back(static_cast<f64>(evidence.mAltNm));

  // Track assigned SPOA haplotype ID (for HSE: path co-segregation)
  data.mHaplotypeIds.push_back(evidence.mAssignedHaplotypeId);
}

// ============================================================================
// MergeAlleleFrom — copy one allele's read data from src into dst_allele slot.
//
// Used by the multi-allelic VariantCall constructor to build a unified
// VariantSupport from N bi-allelic per-variant supports:
//
//   variant[0]: {REF=allele0, ALT=allele1}  → merge REF→0, ALT→1
//   variant[1]: {REF=allele0, ALT=allele1}  → merge ALT→2
//   variant[2]: {REF=allele0, ALT=allele1}  → merge ALT→3
//
// Deduplicates by read name hash (same read can't count twice).
// ============================================================================
void VariantSupport::MergeAlleleFrom(VariantSupport const& src, AlleleIndex const src_allele,
                                     AlleleIndex const dst_allele) {
  if (src_allele >= src.mAlleleData.size()) return;

  EnsureAlleleSlot(dst_allele);
  auto const& src_data = src.mAlleleData[src_allele];
  auto& dst_data = mAlleleData[dst_allele];

  for (auto const& [rname_hash, strand] : src_data.mNameHashes) {
    auto const [iter, inserted] = dst_data.mNameHashes.try_emplace(rname_hash, strand);
    // already seen this read in dst
    if (!inserted) continue;
  }

  // Append PBQ, RMQ, and alignment scores (reads already deduped above by mNameHashes,
  // but the qual vectors are append-only during AddEvidence, so we append here too).
  dst_data.mFwdBaseQuals.insert(dst_data.mFwdBaseQuals.end(), src_data.mFwdBaseQuals.begin(),
                                src_data.mFwdBaseQuals.end());

  dst_data.mRevBaseQuals.insert(dst_data.mRevBaseQuals.end(), src_data.mRevBaseQuals.begin(),
                                src_data.mRevBaseQuals.end());

  dst_data.mMapQuals.insert(dst_data.mMapQuals.end(), src_data.mMapQuals.begin(),
                            src_data.mMapQuals.end());

  dst_data.mAlnScores.insert(dst_data.mAlnScores.end(), src_data.mAlnScores.begin(),
                             src_data.mAlnScores.end());

  dst_data.mSoftClipCount += src_data.mSoftClipCount;

  dst_data.mProperPairIsizes.insert(dst_data.mProperPairIsizes.end(),
                                    src_data.mProperPairIsizes.begin(),
                                    src_data.mProperPairIsizes.end());

  dst_data.mFoldedReadPositions.insert(dst_data.mFoldedReadPositions.end(),
                                       src_data.mFoldedReadPositions.begin(),
                                       src_data.mFoldedReadPositions.end());

  dst_data.mRefNmValues.insert(dst_data.mRefNmValues.end(), src_data.mRefNmValues.begin(),
                               src_data.mRefNmValues.end());

  // Merge new artifact metric vectors
  dst_data.mAlignmentStarts.insert(dst_data.mAlignmentStarts.end(),
                                   src_data.mAlignmentStarts.begin(),
                                   src_data.mAlignmentStarts.end());

  dst_data.mAltNmValues.insert(dst_data.mAltNmValues.end(), src_data.mAltNmValues.begin(),
                               src_data.mAltNmValues.end());

  dst_data.mHaplotypeIds.insert(dst_data.mHaplotypeIds.end(), src_data.mHaplotypeIds.begin(),
                                src_data.mHaplotypeIds.end());
}

// ============================================================================
// Per-Allele Accessors
// ============================================================================
auto VariantSupport::FwdCount(AlleleIndex const idx) const -> usize {
  if (idx >= mAlleleData.size()) return 0;
  return mAlleleData[idx].mFwdBaseQuals.size();
}

auto VariantSupport::RevCount(AlleleIndex const idx) const -> usize {
  if (idx >= mAlleleData.size()) return 0;
  return mAlleleData[idx].mRevBaseQuals.size();
}

auto VariantSupport::TotalAlleleCov(AlleleIndex const idx) const -> usize {
  return FwdCount(idx) + RevCount(idx);
}

auto VariantSupport::TotalSampleCov() const noexcept -> usize {
  usize total = 0;
  for (auto const& entry : mAlleleData) {
    total += entry.mFwdBaseQuals.size() + entry.mRevBaseQuals.size();
  }
  return total;
}

auto VariantSupport::TotalAltCov() const -> usize {
  usize total = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }
  return total;
}

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
auto VariantSupport::RawPosteriorBaseQual(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size()) return 0.0;

  auto const& data = mAlleleData[idx];
  // Per-read representative base quality at the variant position, split by
  // strand. For indels, this is the MINIMUM PBQ across the variant region
  // (one entry per read, NOT per base position).
  auto const& fwd_bq = data.mFwdBaseQuals;
  auto const& rev_bq = data.mRevBaseQuals;

  if (fwd_bq.empty() && rev_bq.empty()) return 0.0;

  f64 log_err = 0.0;  // Σ log10(ε_i)
  f64 log_ok = 0.0;   // Σ log10(1 - ε_i)

  auto const accumulate_bq = [&log_err, &log_ok](std::vector<u8> const& quals) -> void {
    for (auto const qual : quals) {
      f64 const eps = hts::PhredToErrorProb(qual);
      log_err += std::log10(std::max(eps, 1e-300));
      log_ok += std::log10(std::max(1.0 - eps, 1e-300));
    }
  };

  accumulate_bq(fwd_bq);
  accumulate_bq(rev_bq);

  // log-sum-exp: log10(10^a + 10^b) = max(a,b) + log10(1 + 10^(min-max))
  f64 const max_log = std::max(log_err, log_ok);
  f64 const log_sum =
      max_log + std::log10(1.0 + std::pow(10.0, std::min(log_err, log_ok) - max_log));
  f64 const log_posterior_err = log_err - log_sum;
  f64 const posterior_bq = -10.0 * log_posterior_err;

  return posterior_bq;
}

// ============================================================================
// RmsMappingQual (RMQ FORMAT field)
//
// Root-mean-square of the mapping qualities of all reads assigned to this
// allele. This follows the standard INFO/RMQ convention (samtools, GATK).
//
//   RMQ = sqrt( Σ(mapq_i²) / N )
// ============================================================================
auto VariantSupport::RmsMappingQual(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].mMapQuals.empty()) return 0.0;

  auto const& mqs = mAlleleData[idx].mMapQuals;
  auto const square_summer = [](f64 const acc, u8 const mapq) -> f64 {
    return acc + (static_cast<f64>(mapq) * static_cast<f64>(mapq));
  };

  f64 const sum_sq = std::accumulate(mqs.cbegin(), mqs.cend(), 0.0, square_summer);
  return std::sqrt(sum_sq / static_cast<f64>(mqs.size()));
}

// ============================================================================
// StrandBiasLogOR (SB FORMAT field)
//
// Natural log odds ratio for strand bias with Haldane correction (+1 to
// all cells). Coverage-invariant: measures the effect size of strand
// imbalance, not statistical significance.
//
//        |  FWD  |  REV  |
//   -----|-------|-------|
//   REF  | rf    | rr    |
//   ALT  | af    | ar    |      (ALT = sum across all non-REF alleles)
//
//   SB = ln( ((rf+1)(ar+1)) / ((rr+1)(af+1)) )
//
// Positive: ALT enriched on forward strand relative to REF.
// Negative: ALT enriched on reverse strand relative to REF.
// Near 0.0: balanced strands (no bias).
//
// Coverage stability: a 3:1 strand ratio produces SB ≈ 1.1 at any depth
// from 20× to 2000× (vs Fisher SB ranging 3.5→290 over the same range).
// ============================================================================
auto VariantSupport::StrandBiasLogOR() const -> f64 {
  // REF strand counts
  auto const ref_fwd = static_cast<int>(FwdCount(REF_ALLELE_IDX));
  auto const ref_rev = static_cast<int>(RevCount(REF_ALLELE_IDX));

  // ALT strand counts (summed across all non-REF alleles)
  int alt_fwd = 0;
  int alt_rev = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    alt_fwd += static_cast<int>(mAlleleData[i].mFwdBaseQuals.size());
    alt_rev += static_cast<int>(mAlleleData[i].mRevBaseQuals.size());
  }

  // Haldane correction: +1 to all cells handles zero-count edge cases
  // without sentinels. When all counts are zero, returns 0.0 (no bias).
  auto const rf1 = static_cast<f64>(ref_fwd + 1);
  auto const rr1 = static_cast<f64>(ref_rev + 1);
  auto const af1 = static_cast<f64>(alt_fwd + 1);
  auto const ar1 = static_cast<f64>(alt_rev + 1);
  return std::log((rf1 * ar1) / (rr1 * af1));
}

// ============================================================================
// MeanAlnScore
// ============================================================================
auto VariantSupport::MeanAlnScore(AlleleIndex const idx) const -> f64 {
  if (idx >= mAlleleData.size() || mAlleleData[idx].mAlnScores.empty()) return 0.0;

  auto const& scores = mAlleleData[idx].mAlnScores;
  f64 const sum = std::accumulate(scores.cbegin(), scores.cend(), 0.0);
  return sum / static_cast<f64>(scores.size());
}

// ============================================================================
// SoftClipAsymmetry (SCA FORMAT field)
//
// Fraction of soft-clipped reads in ALT minus fraction in REF.
//   SCA = (alt_sc/alt_total) - (ref_sc/ref_total)
//
// Positive values indicate ALT reads are more often soft-clipped than REF,
// which may flag unresolved larger structural events masquerading as smaller
// local variants. Soft-clip status comes from the original whole-genome
// alignment CIGAR (>= 6% of read length), not the genotyper re-alignment.
// ============================================================================
auto VariantSupport::SoftClipAsymmetry() const -> f64 {
  // ALT soft-clip fraction (summed across all non-REF alleles)
  usize alt_sc = 0;
  usize alt_total = 0;

  for (usize i = 1; i < mAlleleData.size(); ++i) {
    alt_sc += mAlleleData[i].mSoftClipCount;
    alt_total += TotalAlleleCov(static_cast<AlleleIndex>(i));
  }

  // REF soft-clip fraction
  auto const has_ref = REF_ALLELE_IDX < mAlleleData.size();
  auto const ref_sc = has_ref ? mAlleleData[REF_ALLELE_IDX].mSoftClipCount : usize{0};
  auto const ref_total = TotalAlleleCov(REF_ALLELE_IDX);

  f64 const alt_frac = alt_total > 0 ? static_cast<f64>(alt_sc) / static_cast<f64>(alt_total) : 0.0;
  f64 const ref_frac = ref_total > 0 ? static_cast<f64>(ref_sc) / static_cast<f64>(ref_total) : 0.0;
  return alt_frac - ref_frac;
}

// ============================================================================
// FragLengthDelta (FLD FORMAT field)
//
// Signed difference in mean insert sizes between ALT-supporting and
// REF-supporting properly-paired reads. Directionality is preserved:
//   FLD = mean_alt_isize − mean_ref_isize
//   Negative → ALT fragments shorter (e.g., cfDNA tumor biology)
//   Positive → ALT fragments longer (e.g., chimeric library artifact)
//
// Insert sizes come from the original whole-genome alignment
// (bam1_t::core.isize), not the genotyper re-alignment.
// Only non-zero insert sizes from properly-paired reads are included.
// ============================================================================
auto VariantSupport::FragLengthDelta() const -> std::optional<f64> {
  // ALT mean insert size (summed across all non-REF alleles)
  f64 alt_isize_sum = 0.0;
  usize alt_pairs = 0;

  for (usize i = 1; i < mAlleleData.size(); ++i) {
    for (auto const isize : mAlleleData[i].mProperPairIsizes) {
      alt_isize_sum += isize;
      ++alt_pairs;
    }
  }

  // REF mean insert size
  f64 ref_isize_sum = 0.0;
  usize ref_pairs = 0;

  if (REF_ALLELE_IDX < mAlleleData.size()) {
    for (auto const isize : mAlleleData[REF_ALLELE_IDX].mProperPairIsizes) {
      ref_isize_sum += isize;
      ++ref_pairs;
    }
  }

  // Untestable if either group has no proper pairs
  if (alt_pairs == 0 || ref_pairs == 0) return std::nullopt;
  f64 const alt_mean = alt_isize_sum / static_cast<f64>(alt_pairs);
  f64 const ref_mean = ref_isize_sum / static_cast<f64>(ref_pairs);
  return alt_mean - ref_mean;
}

// ============================================================================
// MappingQualCohenD (MQCD FORMAT field)
//
// Coverage-normalized effect size comparing mapping qualities of REF vs
// ALT reads. Uses Mann-Whitney U Z/√N to remove √N power amplification.
//
// A mild ALT MAPQ depression (2 units lower) produces MQCD ≈ −0.34 at any
// depth from 20× to 2000× (vs raw Z-score ranging −1.5 to −14.9).
//
// Pools all ALT alleles into a single group (REF vs all-ALT) because the
// test is about whether ALT reads as a class are mismapped.
//
// Returns std::nullopt if either group is empty (untestable).
// Returns 0.0 when test ran but found no bias (genuine zero).
// ============================================================================
auto VariantSupport::MappingQualCohenD() const -> std::optional<f64> {
  // Collect REF mapping qualities
  absl::Span<u8 const> ref_mqs;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    ref_mqs = absl::MakeConstSpan(mAlleleData[REF_ALLELE_IDX].mMapQuals);
  }

  // Pool all ALT mapping qualities into a single group
  std::vector<u8> alt_mqs;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    auto const& mqs = mAlleleData[i].mMapQuals;
    alt_mqs.insert(alt_mqs.end(), mqs.begin(), mqs.end());
  }

  // MannWhitneyEffectSize now returns std::optional<f64>:
  //   nullopt → one or both groups empty (untestable → NaN in VCF)
  //   0.0     → test ran, no bias detected (genuine zero)
  return base::MannWhitneyEffectSize<u8>(ref_mqs, absl::MakeConstSpan(alt_mqs));
}

// ============================================================================
// ReadPosCohenD (RPCD FORMAT field)
//
// Coverage-normalized effect size comparing folded read positions of REF vs
// ALT reads. Folded position = min(p, 1−p) where p = variant_query_pos / read_len.
// Maps both 5' and 3' read edges to 0.0, centers to 0.5.
//
// True variants should be uniformly distributed across read positions.
// Artifacts from 3' quality degradation cluster at read edges (low folded
// position), producing a negative effect size for ALT.
//
// Returns std::nullopt if either group is empty (untestable).
// Returns 0.0 when test ran but found no bias (genuine zero).
// ============================================================================
auto VariantSupport::ReadPosCohenD() const -> std::optional<f64> {
  absl::Span<f64 const> ref_positions;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    ref_positions = absl::MakeConstSpan(mAlleleData[REF_ALLELE_IDX].mFoldedReadPositions);
  }

  std::vector<f64> alt_positions;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    auto const& pos = mAlleleData[i].mFoldedReadPositions;
    alt_positions.insert(alt_positions.end(), pos.begin(), pos.end());
  }

  return base::MannWhitneyEffectSize<f64>(ref_positions, absl::MakeConstSpan(alt_positions));
}

// ============================================================================
// BaseQualCohenD (BQCD FORMAT field)
//
// Coverage-normalized effect size comparing base qualities of REF vs ALT reads.
// Concatenates forward and reverse strand base qualities per allele.
// Detects 8-oxoguanine oxidation artifacts where the miscalled base has
// characteristically low Phred confidence.
//
// Returns std::nullopt if either group is empty (untestable).
// Returns 0.0 when test ran but found no bias (genuine zero).
// ============================================================================
auto VariantSupport::BaseQualCohenD() const -> std::optional<f64> {
  // Collect REF base qualities (fwd + rev concatenated)
  std::vector<u8> ref_bqs;
  if (REF_ALLELE_IDX < mAlleleData.size()) {
    auto const& ref = mAlleleData[REF_ALLELE_IDX];
    ref_bqs.reserve(ref.mFwdBaseQuals.size() + ref.mRevBaseQuals.size());
    ref_bqs.insert(ref_bqs.end(), ref.mFwdBaseQuals.begin(), ref.mFwdBaseQuals.end());
    ref_bqs.insert(ref_bqs.end(), ref.mRevBaseQuals.begin(), ref.mRevBaseQuals.end());
  }

  // Pool all ALT base qualities
  std::vector<u8> alt_bqs;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    auto const& alt = mAlleleData[i];
    alt_bqs.insert(alt_bqs.end(), alt.mFwdBaseQuals.begin(), alt.mFwdBaseQuals.end());
    alt_bqs.insert(alt_bqs.end(), alt.mRevBaseQuals.begin(), alt.mRevBaseQuals.end());
  }

  return base::MannWhitneyEffectSize<u8>(absl::MakeConstSpan(ref_bqs),
                                         absl::MakeConstSpan(alt_bqs));
}

// ============================================================================
// AlleleMismatchDelta (ASMD FORMAT field)
//
// mean(ALT NM) − mean(REF NM) − variant_length, where NM is the edit
// distance of each read against the REF haplotype. Both groups share the
// variant's own edit contribution, but ALT reads also carry the variant's
// inherent edit distance (e.g., a clean 50bp deletion = NM +50 against REF).
// Subtracting variant_length isolates excess noise beyond the expected
// structural difference.
//
// Positive ASMD → ALT reads have more mismatches (chimeric/paralogous signal)
// Near zero     → expected for true variants
// Returns std::nullopt if either group is empty (untestable).
// ============================================================================
auto VariantSupport::AlleleMismatchDelta(usize const variant_length) const -> std::optional<f64> {
  if (REF_ALLELE_IDX >= mAlleleData.size()) return std::nullopt;
  auto const& ref_nms = mAlleleData[REF_ALLELE_IDX].mRefNmValues;
  if (ref_nms.empty()) return std::nullopt;

  auto const ref_sum = std::accumulate(ref_nms.begin(), ref_nms.end(), 0.0);
  auto const ref_mean = ref_sum / static_cast<f64>(ref_nms.size());

  // Pool all ALT NM values
  f64 alt_sum = 0.0;
  usize alt_count = 0;
  for (usize i = 1; i < mAlleleData.size(); ++i) {
    auto const& nms = mAlleleData[i].mRefNmValues;
    alt_sum += std::accumulate(nms.begin(), nms.end(), 0.0);
    alt_count += nms.size();
  }

  if (alt_count == 0) return std::nullopt;
  auto const alt_mean = alt_sum / static_cast<f64>(alt_count);
  // Subtract variant length: a clean indel contributes its own length
  // to NM against REF, but that's the variant itself, not noise.
  auto const adjusted_alt_mean = alt_mean - static_cast<f64>(variant_length);
  return adjusted_alt_mean - ref_mean;
}

// ============================================================================
// Dirichlet-Multinomial Genotype Likelihood Engine
//
// In plain terms: the DM model scores each possible genotype by asking
// "how well do the actual read counts match what this genotype predicts?"
// A true heterozygous site should have roughly 50/50 REF/ALT reads;
// a true homozygous site should be nearly all one allele. The model
// measures this agreement while accounting for the fact that sequencing
// reads are not perfectly independent — correlated errors from PCR,
// flow-cell artifacts, and mismapping mean that extra depth gives
// diminishing returns in confidence, not infinite certainty.
//
// Replaces the legacy per-read binomial model with a count-based
// Dirichlet-Multinomial (DM) distribution that natively handles:
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

namespace {

// ── DM Model Constants ──────────────────────────────────────────────────
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

/// Compute the log-likelihood of observing allele counts under the DM model
/// with concentration parameters dm_alphas.
///
/// In plain terms: given an expected allele mix (dm_alphas, derived from the
/// genotype hypothesis) and the actual read counts, this function returns a
/// single number measuring "how likely is this data under this hypothesis?"
/// More negative = worse fit. The log-gamma terms handle the combinatorial
/// bookkeeping of how many ways the observed counts could arise.
///
///   ln P(c | α) = lnΓ(Σα_i) − lnΓ(N + Σα_i) + Σ[lnΓ(c_i + α_i) − lnΓ(α_i)]
///
/// Uses std::lgamma for numerical stability. All values in natural-log space;
/// conversion to Phred (log10) happens only at the final PL scaling step.
auto LogDirichletMultinomial(std::vector<int> const& allele_counts,
                             std::vector<f64> const& dm_alphas) -> f64 {
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
auto NormalizeToPLs(std::vector<f64> const& log_likelihoods) -> std::vector<int> {
  static constexpr f64 PL_CAP = static_cast<f64>(std::numeric_limits<int>::max()) / 2.0;

  f64 const best_log_lk = *std::ranges::max_element(log_likelihoods);

  f64 const ln_ten = std::numbers::ln10;
  std::vector<int> phred_likelihoods(log_likelihoods.size());
  for (usize index = 0; index < log_likelihoods.size(); ++index) {
    f64 const raw_pl = -10.0 * (log_likelihoods[index] - best_log_lk) / ln_ten;
    phred_likelihoods[index] = static_cast<int>(std::round(std::min(raw_pl, PL_CAP)));
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
/// Templated to accept the private PerAlleleData vector without naming the type.
template <typename AlleleDataVec>
auto PileupLogLikelihood(AlleleDataVec const& allele_data_vec, absl::Span<f64 const> frac_vec,
                         int const num_alleles) -> f64 {
  f64 log_lk = 0.0;
  for (int called_as = 0; called_as < num_alleles; ++called_as) {
    auto const& per_allele = allele_data_vec[called_as];
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

auto VariantSupport::ComputePLs() const -> std::vector<int> {
  auto const num_alleles = static_cast<int>(mAlleleData.size());
  if (num_alleles == 0) return {};

  auto const num_genotypes = num_alleles * (num_alleles + 1) / 2;

  // Build allele count vector: count[i] = total reads assigned to allele i
  std::vector<int> allele_counts(num_alleles, 0);
  for (int allele_idx = 0; allele_idx < num_alleles; ++allele_idx) {
    auto const aidx = static_cast<AlleleIndex>(allele_idx);
    allele_counts[allele_idx] = static_cast<int>(TotalAlleleCov(aidx));
  }

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
// ComputeGQ — Genotype Quality
//
// Standard GATK convention: GQ is the second-smallest PL value.
// After normalization, the smallest PL is always 0 (the called genotype).
// GQ = second-smallest PL, capped at 99.
//
// See: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
// ============================================================================
auto VariantSupport::ComputeGQ(std::vector<int> const& pls) -> int {
  if (pls.size() < 2) return 0;

  // Find minimum and second-minimum PL values
  int min1 = std::numeric_limits<int>::max();
  int min2 = std::numeric_limits<int>::max();
  for (int const val : pls) {
    if (val < min1) {
      min2 = min1;
      min1 = val;
    } else if (val < min2) {
      min2 = val;
    }
  }

  // After normalization, min1 should be 0. GQ = min2 - min1 = min2.
  static constexpr int MAX_GQ = 99;
  return std::min(min2 - min1, MAX_GQ);
}

// ============================================================================
// ComputeContinuousMixtureLods — Continuous Mixture LOD (CMLOD FORMAT field)
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
// Where f_MLE = empirical allele fractions, f_null = redistribute the target
// ALT's mass proportionally among remaining alleles.
//
// A Q40 read supporting ALT contributes ~4.0 to the LOD; a Q10 read
// contributes ~0.5. This natively separates high-confidence signal from
// sequencer noise without requiring any external error model.
//
// Complexity: O(N × K) per ALT allele, where N = total read count.
// ============================================================================
auto VariantSupport::ComputeContinuousMixtureLods() const -> std::vector<f64> {
  auto const num_alleles = static_cast<int>(mAlleleData.size());
  std::vector<f64> lod_scores(num_alleles, 0.0);
  if (num_alleles < 2) return lod_scores;

  // Empirical allele fractions (MLE): f_i = count_i / total_depth
  auto const total_depth = static_cast<f64>(TotalSampleCov());
  if (total_depth == 0.0) return lod_scores;

  std::vector<f64> frac_mle(num_alleles, 0.0);
  for (int index = 0; index < num_alleles; ++index) {
    auto const allele_depth = TotalAlleleCov(static_cast<AlleleIndex>(index));
    frac_mle[index] = static_cast<f64>(allele_depth) / total_depth;
  }

  f64 const log_lk_mle =
      PileupLogLikelihood(mAlleleData, absl::MakeConstSpan(frac_mle), num_alleles);

  // Compute CMLOD for each ALT allele independently
  for (int target_alt = 1; target_alt < num_alleles; ++target_alt) {
    if (TotalAlleleCov(static_cast<AlleleIndex>(target_alt)) == 0) continue;

    auto const frac_null =
        BuildNullFractions(absl::MakeConstSpan(frac_mle), target_alt, num_alleles);

    f64 const log_lk_null =
        PileupLogLikelihood(mAlleleData, absl::MakeConstSpan(frac_null), num_alleles);

    lod_scores[target_alt] = std::max(0.0, log_lk_mle - log_lk_null);
  }

  return lod_scores;
}

// ============================================================================
// EnsureAlleleSlot
// ============================================================================
void VariantSupport::EnsureAlleleSlot(AlleleIndex const idx) {
  if (idx >= mAlleleData.size()) mAlleleData.resize(idx + 1);
}

// ============================================================================
// SupportArray Accessors
// ============================================================================
auto SupportArray::Find(std::string_view sample_name) const -> VariantSupport const* {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto const* const iter = std::ranges::find_if(mItems, pred);
  return iter != mItems.end() ? iter->mData.get() : nullptr;
}

auto SupportArray::FindOrCreate(std::string_view sample_name) -> VariantSupport& {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto* const iter = std::ranges::find_if(mItems, pred);
  if (iter != mItems.end()) return *iter->mData;

  mItems.push_back(
      NamedSupport{.mSampleName = sample_name, .mData = std::make_unique<VariantSupport>()});

  return *mItems.back().mData;
}

auto SupportArray::Extract(std::string_view sample_name) -> std::unique_ptr<VariantSupport> {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto* const iter = std::ranges::find_if(mItems, pred);
  if (iter != mItems.end()) {
    auto extracted_data = std::move(iter->mData);
    mItems.erase(iter);
    return extracted_data;
  }

  return nullptr;
}

// SupportArray inline definitions placed in the header for performance.

}  // namespace lancet::caller
