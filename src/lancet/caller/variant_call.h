#ifndef SRC_LANCET_CALLER_VARIANT_CALL_H_
#define SRC_LANCET_CALLER_VARIANT_CALL_H_

#include "lancet/base/sequence_complexity.h"
#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/sample_format_data.h"
#include "lancet/caller/support_array.h"
#include "lancet/caller/variant_support.h"
#include "lancet/core/sample_info.h"

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"

#include <algorithm>
#include <compare>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::caller {

using VariantID = u64;

// ============================================================================
// VariantCall: a finalized VCF record, potentially multi-allelic.
//
// Multi-allelic support:
//   The constructor takes a single multiallelic RawVariant from the Genotyper.
//   It unpacks the pre-mapped SupportArray layout into a unified VCF trace:
//
//     Input:  {chr1:100 A→T, A→G}  (One Multiallelic RawVariant)
//     Output: VCF record: chr1 100 . A T,G ... (comma-separated ALTs)
//
// FORMAT fields (per-sample): see vcf_formatter.h::FORMAT_HEADER for the
// authoritative single-source-of-truth list and field ordering rationale.
//
//   GT   - Genotype (e.g., 0/1, 1/2 for multi-allelic)
//   AD   - Number=R: read depth per allele (REF, ALT1, ALT2, ...)
//   ADF  - Number=R: forward strand depth per allele
//   ADR  - Number=R: reverse strand depth per allele
//   DP   - Total read depth
//   RMQ  - Number=R: RMS mapping quality per allele
//   NPBQ - Number=R: normalized posterior base quality per allele (PBQ/N)
//   SB   - Number=1: Strand bias log odds ratio (Haldane-corrected)
//   SCA  - Number=1: Soft Clip Asymmetry (ALT - REF soft-clip fraction)
//   FLD  - Number=1: Fragment Length Delta (mean ALT isize − mean REF isize, signed)
//   RPCD - Number=1: Read Position Cohen's D (folded position effect size)
//   BQCD - Number=1: Base Quality Cohen's D (base quality effect size)
//   MQCD - Number=1: Mapping Quality Cohen's D (MAPQ effect size)
//   ASMD - Number=1: Allele-Specific Mismatch Delta (mean ALT NM − mean REF NM − variant_length)
//   SDFC - Number=1: Site Depth Fold Change (sample DP / per-sample window mean coverage)
//   PRAD - Number=1: Polar Radius log10(1 + sqrt(AD_Ref² + AD_Alt²))
//   PANG - Number=1: Polar Angle atan2(AD_Alt, AD_Ref) in radians
//   CMLOD - Number=A: Continuous Mixture LOD per ALT (quality-weighted)
//   FSSE  - Number=1: Fragment Start Shannon Entropy [0,1] (ALT start position diversity)
//   AHDD  - Number=1: ALT-Haplotype Discordance Delta (ALT reads vs own haplotype)
//   HSE   - Number=1: Haplotype Segregation Entropy [0,1] (ALT path concentration)
//   PDCV  - Number=1: Path Depth Coefficient of Variation (graph coverage uniformity)
//   PL    - Number=G: Phred-scaled genotype likelihoods (Dirichlet-Multinomial)
//   GQ    - Genotype quality (second-lowest DM PL, capped at 99)
// ============================================================================
class VariantCall {
 public:
  using Samples = absl::Span<core::SampleInfo const>;

  // SampleFormatData: per-sample FORMAT field payload.
  // See sample_format_data.h for the full class definition.
  using SampleFormatData = caller::SampleFormatData;

  // Native multi-allelic constructor
  using SupportsByVariant = absl::flat_hash_map<RawVariant const*, SupportArray>;
  VariantCall(RawVariant const* var, SupportsByVariant const& all_supports, Samples samps,
              usize window_length);

  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIndex; }
  [[nodiscard]] auto ChromName() const -> std::string_view { return mChromName; }
  [[nodiscard]] auto StartPos1() const -> usize { return mStartPos1; }
  [[nodiscard]] auto RefAllele() const -> std::string_view { return mRefAllele; }
  [[nodiscard]] auto AltAlleles() const -> absl::Span<std::string const> { return mAltAlleles; }
  [[nodiscard]] auto NumAltAlleles() const -> usize { return mAltAlleles.size(); }
  [[nodiscard]] auto VariantLengths() const -> absl::Span<i64 const> { return mVariantLengths; }
  [[nodiscard]] auto Quality() const -> f64 { return mSiteQuality; }
  [[nodiscard]] auto State() const -> AlleleState { return mState; }
  [[nodiscard]] auto Categories() const -> absl::Span<AlleleType const> { return mCategories; }
  [[nodiscard]] auto IsMultiallelic() const -> bool { return mIsMultiallelic; }

  [[nodiscard]] auto NumSamples() const -> usize { return mSampleGenotypes.size(); }
  [[nodiscard]] auto Identifier() const -> VariantID { return mVariantId; }
  [[nodiscard]] auto TotalCoverage() const -> usize { return mTotalSampleCov; }

  /// Returns true if any sample has ALT allele support.
  /// Used by VariantStore to filter out zero-evidence calls.
  [[nodiscard]] auto HasAltSupport() const -> bool { return mHasAltSupport; }

  [[nodiscard]] auto AsVcfRecord() const -> std::string;

  // ============================================================================
  // VARIANT STORE EXTENSIONS (* ALLELE OVERLAPS)
  // ============================================================================
  [[nodiscard]] auto IsDeletion() const -> bool {
    return std::ranges::any_of(mCategories,
                               [](AlleleType type) { return type == AlleleType::DEL; });
  }

  [[nodiscard]] auto GetMaxDeletionLength() const -> i64 {
    static constexpr auto SELECT_MAX = [](i64 max_so_far, i64 curr_len) {
      return std::max(max_so_far, curr_len);
    };

    static constexpr auto FILTER_DEL = [](AlleleType var_type, i64 len_val) {
      return var_type == AlleleType::DEL ? len_val : 0;
    };

    return std::transform_reduce(mCategories.cbegin(), mCategories.cend(), mVariantLengths.cbegin(),
                                 i64{0}, SELECT_MAX, FILTER_DEL);
  }

  [[nodiscard]] auto RefLength() const -> usize { return mRefAllele.length(); }

  // ============================================================================
  // VARIANT IDENTITY & ORDERING DESIGN
  // ============================================================================
  // Identity (operator==) and ordering (operator<) both use the same
  // conceptual key: CHROM + POS + REF (locus-level, ALTs excluded).
  //
  // WHY NO ALTs IN IDENTITY:
  //   Overlapping genomic windows independently assemble the same locus.
  //   Window A might produce chr1:100 A→T at 80x, while Window B produces
  //   chr1:100 A→T,G at 120x. These are the same locus — the higher-coverage
  //   window assembled a more complete multi-allelic picture. With ALTs in the
  //   hash, they'd be treated as *different* variants, producing duplicate VCF
  //   records at the same locus (poor practice per VCF v4.5 spec). With
  //   CHROM+POS+REF identity, VariantStore dedup keeps only the higher-coverage
  //   call, ensuring at most one VCF record per locus.
  //
  // DOWNSTREAM INVARIANT (VariantStore):
  //   VariantStore uses Identifier() (== mVariantId) as flat_hash_map key.
  //   When a duplicate is found (same CHROM+POS+REF), the variant with higher
  //   TotalCoverage() replaces the existing one. This means after dedup, at
  //   most one variant per locus exists, so the TotalCov and mVariantId
  //   tiebreakers in operator< are purely for mathematical completeness.
  //
  // WARNING: Do NOT add ALT alleles to HashRawVariant or change the fields
  //   used by operator== / operator< without updating VariantStore's dedup
  //   logic. Mismatched identity semantics will produce duplicate VCF records
  //   and break tabix indexing.
  // ============================================================================
  friend auto operator==(VariantCall const& lhs, VariantCall const& rhs) -> bool {
    return lhs.mVariantId == rhs.mVariantId;
  }

  friend auto operator<(VariantCall const& lhs, VariantCall const& rhs) -> bool {
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;
    if (lhs.mStartPos1 != rhs.mStartPos1) return lhs.mStartPos1 < rhs.mStartPos1;
    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;

    if (lhs.mTotalSampleCov != rhs.mTotalSampleCov) {
      return lhs.mTotalSampleCov < rhs.mTotalSampleCov;
    }

    return lhs.mVariantId < rhs.mVariantId;
  }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  u64 mVariantId;
  usize mChromIndex;
  usize mStartPos1;
  usize mTotalSampleCov{0};
  f64 mSiteQuality{0};

  // Raw variant pointer — valid for the lifetime of the VariantBuilder's
  // btree_set (which outlives all VariantCall instances). Needed by
  // BuildFormatFields to access mNumTotalHaps (HSE) and mMaxPathCv (PDCV).
  RawVariant const* mRawVariant{nullptr};

  usize mWindowLength;

  std::string mChromName;
  std::string mRefAllele;
  std::string mInfoField;

  std::vector<i64> mVariantLengths;
  std::vector<std::string> mAltAlleles;
  std::vector<AlleleType> mCategories;
  std::vector<SampleFormatData> mSampleGenotypes;

  // ============================================================================
  // Sequence complexity (11 coverage-invariant features, from RawVariant)
  // ============================================================================
  // ============================================================================
  // Graph complexity metrics (from RawVariant, transcribed)
  // ============================================================================
  base::SequenceComplexity mSeqCx;
  GraphMetrics mGraphCx;

  // ── 1B Align ────────────────────────────────────────────────────────────
  AlleleState mState = AlleleState::NONE;
  bool mIsMultiallelic = false;
  bool mHasAltSupport = false;

  /// Site Depth Fold Change: sample DP / per-sample window mean coverage.
  /// Returns nullopt if window coverage is zero (renders as "." in VCF).
  [[nodiscard]] auto SiteDepthFoldChange(core::SampleInfo const& sinfo, usize sample_dp) const
      -> std::optional<f64> {
    auto const window_cov = sinfo.MeanCoverage(mWindowLength);
    if (window_cov <= 0.0) return std::nullopt;
    return static_cast<f64>(sample_dp) / window_cov;
  }

  // ============================================================================
  // Evidence collection (shared by both constructors)
  // ============================================================================

  /// Common finalization after evidence is assembled: builds FORMAT, state, and INFO fields.
  void Finalize(SupportArray const& evidence, Samples samps);

  // ============================================================================
  // Modular field builders
  // ============================================================================
  struct PerAlleleMetrics {
    // ── 8B Align ────────────────────────────────────────────────────────────
    std::string mAd;
    std::string mAdf;
    std::string mAdr;
    std::string mRmq;
    std::string mNpbq;
  };

  /// Build per-sample FORMAT components and track multi-allelic site qualities.
  void BuildFormatFields(SupportArray const& evidence, Samples samps, bool case_ctrl_mode);

  // Sets mGenotypeIndices from PL values.
  static void AssignGenotype(SampleFormatData& sample, absl::Span<u32 const> pls,
                             usize num_alleles);

  /// Convert a GL index back to genotype string (e.g., GL=4 with k=3 → "1/2")
  [[nodiscard]] static auto GenotypeFromGLIndex(usize gl_index, usize num_alleles) -> std::string;

  void UpdateSiteQuality(core::SampleInfo const& sinfo, VariantSupport const* support,
                         SupportArray const& evidence, Samples samps, bool case_ctrl_mode,
                         absl::Span<u32 const> pls);

  /// Somatic log odds ratio: case ALT enrichment vs control.
  [[nodiscard]] static auto SomaticLogOddsRatio(core::SampleInfo const& curr,
                                                SupportArray const& supports, Samples samps) -> f64;

  /// Extract per-allele metrics from VariantSupport.
  static void AssignPerAlleleMetrics(SampleFormatData& sample, VariantSupport const* support,
                                     usize num_alleles);

  /// Compute SHARED/CTRL/CASE/UNKNOWN state from evidence.
  /// In non-case-control mode (i.e. control-only), state is always UNKNOWN.
  void ComputeState(SupportArray const& evidence, Samples samps, bool case_ctrl_mode);

  /// Assemble the INFO field string (TYPE, LENGTH,
  /// optional state prefix, complexity annotations).
  void BuildInfoField(bool case_ctrl_mode);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_CALL_H_
