#ifndef SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
#define SRC_LANCET_CALLER_VARIANT_SUPPORT_H_

#include "lancet/base/mann_whitney.h"
#include "lancet/base/types.h"
#include "lancet/caller/per_allele_data.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <functional>
#include <numeric>
#include <optional>
#include <vector>

#include <cmath>

namespace lancet::caller {

// ============================================================================
// AlleleIndex: integer-based allele identifier for multi-allelic support
//
// Unlike the previous binary Allele enum {REF, ALT}, AlleleIndex supports
// arbitrary numbers of alleles:
//   0 = REF, 1 = ALT1, 2 = ALT2, ...
//
// This is the foundation for multi-allelic VCF output (comma-separated ALTs)
// and for Phase 2 graph alignment where a read's path through the POA graph
// directly identifies which allele it supports at each variant position.
// ============================================================================
using AlleleIndex = u8;
static constexpr AlleleIndex REF_ALLELE_IDX = 0;

// ============================================================================
// VariantSupport: per-sample allele evidence aggregator
//
// Collects read-level evidence for each allele at a variant site and computes
// aggregate metrics for VCF FORMAT fields. Designed for multi-allelic sites.
//
//  ┌──────────────────┐
//  │  Read Evidence   │  ← one per aligned read at this variant
//  │  {allele, quals} │
//  └────────┬─────────┘
//           ▼
//  ┌──────────────────┐
//  │  PerAlleleData   │  ← one per allele (REF, ALT1, ALT2...)
//  │  [fwd_bq, rev_bq │     stored in dense vector indexed by AlleleIndex
//  │   map_quals, ...]│
//  └────────┬─────────┘
//           ▼
//  ┌──────────────────┐
//  │ Aggregation      │  → PL, NPBQ, RMQ, SB, GQ for VCF FORMAT
//  └──────────────────┘
//
// The PerAlleleData is stored as std::vector<PerAlleleData> indexed directly
// by AlleleIndex. This is efficient because AlleleIndex is a dense, zero-based
// integer (typically 0-3 for most sites, never more than ~8).
// ============================================================================
class VariantSupport {
 public:
  VariantSupport() = default;

  struct ReadEvidence {
    // ── 8B Align ────────────────────────────────────────────────────────────
    i64 mInsertSize;      // template length from original alignment (for FLD)
    i64 mAlignmentStart;  // fragment genomic start position (for FSSE)
    f64 mAlnScore;        // normalized alignment score to the assigned haplotype
    f64 mFoldedReadPos;   // 0.0=read edge, 0.5=read center (for RPCD)

    // ── 4B Align ────────────────────────────────────────────────────────────
    u32 mRnameHash;            // hash of read name (for dedup)
    u32 mRefNm;                // edit distance to REF haplotype (for ASMD)
    u32 mOwnHapNm;             // edit distance to assigned haplotype (for AHDD)
    u32 mAssignedHaplotypeId;  // SPOA path index this read was assigned to (for HSE)

    // ── 1B Align ────────────────────────────────────────────────────────────
    AlleleIndex mAllele;  // which allele this read supports
    Strand mStrand;       // forward or reverse strand
    u8 mBaseQual;         // representative PBQ (min across variant region for indels)
    u8 mMapQual;          // original mapping quality of the read
    bool mIsSoftClipped;  // soft-clip bases >= 6% of read length in original alignment
    bool mIsProperPair;   // properly paired in original alignment
  };

  void AddEvidence(ReadEvidence const& evidence);

  // ── Per-Allele Accessors ──
  [[nodiscard]] auto FwdCount(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto RevCount(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto TotalAlleleCov(AlleleIndex idx) const -> usize;
  [[nodiscard]] auto TotalSampleCov() const noexcept -> usize;
  [[nodiscard]] auto NumAlleles() const noexcept -> usize { return mAlleleData.size(); }

  // Convenience shims matching the old REF/ALT interface (for VariantCall)
  [[nodiscard]] auto TotalRefCov() const -> usize { return TotalAlleleCov(REF_ALLELE_IDX); }
  [[nodiscard]] auto TotalAltCov() const -> usize;

  // ── Aggregate Metrics for VCF FORMAT Fields ──

  // Normalized Posterior Base Quality (NPBQ FORMAT field):
  // Raw Bayesian posterior base quality divided by the number of supporting reads,
  // yielding the effective per-read quality contribution. Coverage-invariant:
  // returns ~30 for Q30 reads at any depth from 20× to 2000×.
  //
  // Given N reads with Phred qualities Q_i:
  //   ε_i = 10^(-Q_i/10)
  //   log_err  = Σ log10(ε_i)
  //   log_ok   = Σ log10(1 - ε_i)
  //   posterior = 10^log_err / (10^log_err + 10^log_ok)
  //   raw PBQ = -10 * log10(posterior)
  //
  // Returns the raw uncapped PBQ value. Caller divides by allele_depth
  // to produce NPBQ for VCF output. Not capped at 255.
  [[nodiscard]] auto RawPosteriorBaseQual(AlleleIndex idx) const -> f64;

  // RMS mapping quality: sqrt(mean(mapq_i^2)) for reads supporting this allele.
  // Per-allele FORMAT field (RMQ), follows the standard samtools/GATK convention.
  // Coverage stability: bounded by MAPQ ceiling (60). Converges rapidly;
  // effectively coverage-invariant above ~10 reads per allele.
  [[nodiscard]] auto RmsMappingQual(AlleleIndex idx) const -> f64;

  // Strand Bias Log Odds Ratio (SB FORMAT field):
  // Natural log of the odds ratio from a 2×2 REF/ALT × FWD/REV table,
  // with Haldane +1 correction to handle zero cells without sentinels.
  //
  //   SB = ln( ((ref_fwd+1)(alt_rev+1)) / ((ref_rev+1)(alt_fwd+1)) )
  //
  // Coverage-invariant: measures the *effect size* of strand imbalance,
  // not statistical significance. A 3:1 strand ratio produces SB ≈ 1.1
  // at any depth from 20× to 2000×. Preserves directionality:
  //   Positive = ALT enriched on forward strand relative to REF
  //   Negative = ALT enriched on reverse strand relative to REF
  [[nodiscard]] auto StrandBiasLogOR() const -> f64;

  // Mean normalized alignment score for reads supporting this allele.
  [[nodiscard]] auto MeanAlnScore(AlleleIndex idx) const -> f64;

  // Soft Clip Asymmetry (SCA FORMAT field): fraction of soft-clipped reads in ALT
  // minus fraction in REF. Detects unresolved larger variant events masquerading
  // as smaller local calls.
  //   SCA = (alt_sc / alt_total) - (ref_sc / ref_total)
  [[nodiscard]] auto SoftClipAsymmetry() const -> f64;

  // Fragment Length Delta (FLD FORMAT field): signed difference in mean insert
  // sizes between ALT-supporting and REF-supporting read pairs. Preserves
  // directionality: negative = ALT fragments shorter (cfDNA biology),
  // positive = ALT fragments longer (chimeric artifact).
  //   FLD = mean_alt_isize - mean_ref_isize
  //   Returns std::nullopt when no proper pairs exist in either group.
  [[nodiscard]] auto FragLengthDelta() const -> std::optional<f64>;

  // Mapping Quality Cohen's D (MQCD FORMAT field): coverage-normalized effect
  // size comparing mapping qualities of REF-supporting vs ALT-supporting reads.
  // Uses Mann-Whitney U test Z/√N to remove √N power amplification.
  //
  // Coverage-invariant: a mild ALT MAPQ depression produces MQCD ≈ −0.34
  // at any depth from 20× to 2000×.
  //   Negative → ALT reads are systematically lower MAPQ (paralogous mismapping)
  //   Near zero → no systematic difference (expected for true variants)
  //   Returns std::nullopt when either group is empty (untestable).
  //   Returns 0.0 when test ran but found no bias (genuine zero).
  [[nodiscard]] auto MappingQualCohenD() const -> std::optional<f64>;

  // Read Position Cohen's D (RPCD FORMAT field): coverage-normalized effect
  // size comparing folded read positions of REF vs ALT reads. Folded position
  // maps both 5' and 3' read edges to 0.0, centers to 0.5.
  //   Negative → ALT allele systematically at read edges (artifact signal)
  //   Near zero → uniform distribution (expected for true variants)
  //   Returns std::nullopt when either group is empty (untestable).
  //   Returns 0.0 when test ran but found no bias (genuine zero).
  [[nodiscard]] auto ReadPosCohenD() const -> std::optional<f64>;

  // Base Quality Cohen's D (BQCD FORMAT field): coverage-normalized effect
  // size comparing per-allele base qualities (REF vs ALT). Detects 8-oxoG
  // oxidation artifacts where the miscalled base has characteristically low
  // Phred confidence.
  //   Negative → ALT allele has systematically lower base quality
  //   Returns std::nullopt when either group is empty (untestable).
  //   Returns 0.0 when test ran but found no bias (genuine zero).
  [[nodiscard]] auto BaseQualCohenD() const -> std::optional<f64>;

  // Allele-Specific Mismatch Delta (ASMD FORMAT field):
  // mean(ALT NM) − mean(REF NM) − variant_length, where NM is edit distance
  // to the REF haplotype. ALT reads carry the variant's inherent edit distance
  // (e.g., a clean 50bp deletion = NM +50 against REF). Subtracting
  // variant_length removes this expected structural difference, isolating only
  // excess noise from chimeric or paralogous mismapping.
  //   Returns std::nullopt when either group is empty (untestable).
  //   Returns 0.0 when no excess mismatch detected (genuine zero).
  [[nodiscard]] auto AlleleMismatchDelta(usize variant_length = 0) const -> std::optional<f64>;

  // Fragment Start Shannon Entropy (FSSE FORMAT field):
  // Measures spatial diversity of ALT read alignment start positions.
  // Catches PCR duplicate jackpot artifacts that survive upstream
  // MarkDuplicates, which is binary, global-coordinate, and allele-blind.
  //
  // Five failure modes MarkDuplicates misses:
  //   1. Exonuclease fraying: 1–4 bp enzymatic nibbling staggers clone
  //      starts past exact-coordinate dedup. 3 bp binning absorbs this.
  //   2. Alignment jitter: BWA-MEM maps complex variants chaotically,
  //      giving clones divergent global POS. FSSE operates post-SPOA.
  //   3. Birthday Paradox: at ≥5,000× Poisson crowding, coordinate
  //      collisions among independent fragments are expected.
  //   4. Representative read roulette: Picard keeps the highest-BQ read,
  //      potentially amplifying a single error to 100% VAF.
  //   5. Continuous vs binary: Shannon entropy preserves distributional
  //      information that binary dedup destroys.
  //
  // Computation: bin ALT starts into 3 bp windows → H = −Σ pᵢ log₂(pᵢ),
  // normalized by log₂(min(N, 20)) to [0, 1].
  //   Returns std::nullopt when < 3 ALT reads (insufficient data).
  [[nodiscard]] auto ComputeFSSE() const -> std::optional<f64>;

  // ALT-Haplotype Discordance Delta (AHDD FORMAT field):
  // mean(ALT reads' NM against their assigned ALT haplotype) −
  // mean(REF reads' NM against REF haplotype).
  //
  // Both NM values come from the same mOwnHapNmValues vector, which
  // stores each read's edit distance against its assigned haplotype
  // (after the mAltNm → mOwnHapNm rename and guard removal — see
  // genotyper.cpp). This ensures an apples-to-apples comparison: both
  // groups use the same ComputeEditDistance() path.
  //
  // ASMD asks "do ALT reads fit the REF haplotype?" — AHDD asks "do ALT
  // reads fit their OWN assembled haplotype?" High AHDD means the SPOA
  // consensus doesn't match its supporting reads, signaling an assembly
  // hallucination.
  //   Returns std::nullopt when either group is empty (untestable).
  [[nodiscard]] auto ComputeAHDD() const -> std::optional<f64>;

  // Haplotype Segregation Entropy (HSE FORMAT field):
  // Shannon entropy of ALT reads' SPOA haplotype path assignments,
  // normalized to [0, 1].
  //
  // True somatic variants concentrate ALT reads on a single haplotype
  // (HSE ≈ 0). Random sequencer errors scatter across multiple assembly
  // paths (HSE > 0.5). With only one haplotype, segregation is
  // undefined — entropy requires at least two categories.
  //   Returns std::nullopt when < 3 ALT reads or total_haplotypes < 2.
  [[nodiscard]] auto ComputeHSE(usize total_haplotypes) const -> std::optional<f64>;

  // ── Multi-Allelic Genotype Likelihoods (Dirichlet-Multinomial Model) ──

  // Computes Phred-scaled genotype likelihoods (PLs) for all possible diploid
  // genotypes at a multi-allelic site using the Dirichlet-Multinomial (DM)
  // distribution instead of the standard per-read binomial model.
  //
  // Intuition: for each possible genotype (e.g., REF/REF, REF/ALT, ALT/ALT),
  // the DM model asks "given what this genotype predicts about allele counts,
  // how well do the observed read counts match?" Unlike the standard model
  // which treats each read as an independent coin flip, the DM model accounts
  // for the reality that reads from the same sequencing run share correlated
  // errors — so adding more reads gives diminishing returns in confidence,
  // preventing extreme PL values at ultra-high depth.
  //
  // The DM model generalizes the overdispersed Beta-Binomial to K dimensions,
  // absorbing correlated sequencing errors at ultra-high depths
  // (2500x+). PLs plateau without depth-division clamping.
  //
  // For K alleles (0=REF, 1..K-1=ALTs), there are K*(K+1)/2 diploid genotypes.
  // For each genotype (a1, a2), expected allele fractions μ_i are:
  //   Homozygous (a1==a2):  μ[a1] = 1 − ε,     μ[other] = ε/(K−1)
  //   Heterozygous:         μ[a1] = μ[a2] = (1−ε)/2, μ[other] = ε/(K−1)
  //
  // DM log-likelihood with precision M = (1−ρ)/ρ and α_i = M × μ_i:
  //   ln P(counts | μ, M) = lnΓ(Σα) − lnΓ(N + Σα) + Σ[lnΓ(c_i + α_i) − lnΓ(α_i)]
  //
  // Genotypes ordered per VCF 4.3 spec §1.6.2:
  //   GL_index(i,j) = j*(j+1)/2 + i   where i ≤ j
  //
  // Returns: vector of K*(K+1)/2 Phred-scaled likelihoods (best genotype PL=0).
  //
  // See docs/guides/variant_discovery_genotyping.md for the full DM derivation.
  [[nodiscard]] auto ComputePLs(usize num_alleles) const -> absl::InlinedVector<u32, 6>;

  // Genotype Quality (GQ): confidence in the called genotype.
  // Second-smallest PL value, capped at 99. Standard GATK convention.
  //   GQ = second_min(PLs) − min(PLs)   [min is always 0 after normalization]
  //
  // Value now derives from DM-based PLs, which plateau at high depth.
  // See: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
  [[nodiscard]] static auto ComputeGQ(absl::Span<u32 const> phred_likelihoods) -> u32;

  // ── Continuous Mixture Log-Odds (CMLOD FORMAT field) ──

  // Computes a per-ALT-allele Log-Odds score comparing the likelihood of the
  // observed read pileup under the MLE allele fractions vs. a null hypothesis
  // where the target ALT fraction is forced to zero.
  //
  // Intuition: CMLOD asks "is there any evidence this variant is real?" by
  // comparing two stories about the data — one where the variant exists at
  // its observed frequency, and one where it doesn't exist at all. Each
  // read's vote is weighted by its sequencing quality: a high-quality Q40
  // read is a strong vote (~4 points), while a noisy Q10 read barely counts
  // (~0.5 points). A higher CMLOD means more confident evidence that the
  // variant is real. This is especially useful for detecting low-frequency
  // variants (e.g., 2–15% VAF) where the standard PL model — which only
  // considers 0%, 50%, and 100% — cannot represent the true frequency.
  //
  // Unlike DM-based PLs (which evaluate discrete diploid genotype states), the
  // CMLOD operates over continuous frequency space and integrates exact per-read
  // base qualities. This separates Q40 reads from Q10 noise, solving
  // the identifiability problem where count-based models can't distinguish a 2%
  // true mosaic from a 2% systematic artifact.
  //
  // For each read assigned to allele s with error probability ε:
  //   P(read | f) = Σ_t f[t] × P(read called as s | true origin = t)
  //   P(s|t) = (1−ε) if s==t, else ε/(K−1)
  //
  //   CMLOD[alt] = max(0, LL(f_MLE) − LL(f_null))   (log10 scale)
  //
  // Returns: vector of K f64 values. Index 0 (REF) is always 0.0.
  //          Index i>0 is the CMLOD for ALT allele i.
  //
  // See docs/guides/variant_discovery_genotyping.md for the full derivation.
  [[nodiscard]] auto ComputeContinuousMixtureLods(usize num_alleles) const -> std::vector<f64>;

  // Copy allele data from `src` allele `src_allele` into `dst_allele` slot
  // in this object. Used for multi-allelic merging: each bi-allelic variant
  // has alleles {0=REF, 1=ALT}. When merging N variants at a locus, we remap
  // variant[i]'s ALT(1) → merged allele (i+1).
  void MergeAlleleFrom(VariantSupport const& src, AlleleIndex src_allele, AlleleIndex dst_allele);

 private:
  // Dense vector indexed by AlleleIndex: mAlleleData[0]=REF, [1]=ALT1, ...
  std::vector<PerAlleleData> mAlleleData;

  // Grow the vector to accommodate a new allele index.
  void EnsureAlleleSlot(AlleleIndex idx);

  // ── Template Helpers (defined in header for instantiation) ─────────────

  /// Pool values from all non-REF alleles into a single vector.
  template <typename ValueType, typename FieldAccessor>
  [[nodiscard]] auto PoolAltValues(FieldAccessor const& field_getter) const
      -> std::vector<ValueType> {
    std::vector<ValueType> pooled;
    for (usize aidx = 1; aidx < mAlleleData.size(); ++aidx) {
      auto const& field = field_getter(mAlleleData[aidx]);
      pooled.insert(pooled.end(), field.begin(), field.end());
    }
    return pooled;
  }

  /// Compute REF-vs-pooled-ALT effect size using a field accessor.
  /// Pattern A: wraps MannWhitneyEffectSize for Cohen's D metrics.
  template <typename ValueType, typename FieldAccessor>
  [[nodiscard]] auto RefVsAltEffectSize(FieldAccessor&& field_getter) const -> std::optional<f64> {
    if (REF_ALLELE_IDX >= mAlleleData.size()) return std::nullopt;
    auto const ref_vals = absl::MakeConstSpan(field_getter(mAlleleData[REF_ALLELE_IDX]));
    auto const alt_vals = PoolAltValues<ValueType>(std::forward<FieldAccessor>(field_getter));
    return base::MannWhitneyEffectSize<ValueType>(ref_vals, absl::MakeConstSpan(alt_vals));
  }

  /// Compute mean(pooled ALT values) − mean(REF values) with an optional offset.
  /// Pattern B: uses std::transform_reduce with explicit f64 promotion.
  template <typename ValueType, typename FieldAccessor>
  [[nodiscard]] auto MeanAltMinusRef(FieldAccessor&& field_getter, f64 offset = 0.0) const
      -> std::optional<f64> {
    if (REF_ALLELE_IDX >= mAlleleData.size()) return std::nullopt;
    auto const& ref_vals = field_getter(mAlleleData[REF_ALLELE_IDX]);
    if (ref_vals.empty()) return std::nullopt;

    // std::transform_reduce with explicit f64 cast avoids narrowing from u32→f64 implicit
    // promotion in std::accumulate. Both produce identical results for u32 inputs (f64 has
    // 52-bit mantissa, sufficient for exact integer sums up to 2^53), but transform_reduce
    // makes the widening intent explicit at the call site.
    auto const ref_sum =
        std::transform_reduce(ref_vals.begin(), ref_vals.end(), 0.0, std::plus<>{},
                              [](ValueType value) { return static_cast<f64>(value); });
    auto const ref_mean = ref_sum / static_cast<f64>(ref_vals.size());

    auto const alt_vals = PoolAltValues<ValueType>(std::forward<FieldAccessor>(field_getter));
    if (alt_vals.empty()) return std::nullopt;

    auto const alt_sum =
        std::transform_reduce(alt_vals.begin(), alt_vals.end(), 0.0, std::plus<>{},
                              [](ValueType value) { return static_cast<f64>(value); });
    auto const alt_mean = alt_sum / static_cast<f64>(alt_vals.size());
    return (alt_mean - offset) - ref_mean;
  }

  /// Compute normalized Shannon entropy of pooled ALT values.
  /// Pattern C: bins values via bin_func, normalizes by log2(max_bins).
  template <typename ValueType, typename BinKeyType, typename FieldAccessor, typename BinFunc>
  [[nodiscard]] auto AltPooledEntropy(FieldAccessor const& field_getter, BinFunc const& bin_func,
                                      f64 max_bins) const -> std::optional<f64> {
    auto const pooled = PoolAltValues<ValueType>(field_getter);
    if (pooled.size() < 3) return std::nullopt;

    absl::flat_hash_map<BinKeyType, usize> bin_counts;
    for (auto const& value : pooled) {
      bin_counts[bin_func(value)]++;
    }

    auto const total = static_cast<f64>(pooled.size());
    f64 entropy = 0.0;
    for (auto const& [_key, count] : bin_counts) {
      auto const prob = static_cast<f64>(count) / total;
      entropy -= prob * std::log2(prob);
    }

    auto const max_entropy = std::log2(std::min(total, max_bins));
    return max_entropy > 0.0 ? (entropy / max_entropy) : 0.0;
  }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SUPPORT_H_
