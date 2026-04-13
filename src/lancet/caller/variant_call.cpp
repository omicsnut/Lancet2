#include "lancet/caller/variant_call.h"

#include "lancet/base/assert.h"
#include "lancet/base/polar_coords.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/label.h"

#include "absl/hash/hash.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

// Variant identity is defined by CHROM+POS+REF (locus-level, no ALTs).
// Uses mChromIndex (integer, always correctly set) instead of mChromName for robustness.
// See the identity design note on VariantCall::operator< for full rationale.
[[nodiscard]] inline auto HashRawVariant(lancet::caller::RawVariant const* var) -> u64 {
  return static_cast<u64>(absl::HashOf(var->mChromIndex, var->mGenomeChromPos1, var->mRefAllele));
}

// Convert std::optional<f64> from VariantSupport to std::optional<f32> for SampleGenotypeData.
// Preserves nullopt (untestable → "." in VCF) and narrows f64 → f32 for storage.
[[nodiscard]] inline auto ToOptF32(std::optional<f64> const& val) -> std::optional<f32> {
  return val.has_value() ? std::optional(static_cast<f32>(val.value())) : std::nullopt;
}

// Format a std::optional<f32> for VCF output, emitting "." for missing values.
// VCF 4.5 spec: "." = missing value. This is the single point where the
// optional → dot conversion happens for all FORMAT fields.
// NOTE: Cannot use IEEE NaN here because -ffast-math (cmake/defaults.cmake)
// implies -ffinite-math-only, which optimizes away std::isnan().
[[nodiscard]] inline auto FormatOptional(std::optional<f32> const& val, char const* fmt_spec)
    -> std::string {
  if (!val.has_value()) return ".";
  return fmt::format(fmt::runtime(fmt_spec), val.value());
}

}  // namespace

namespace lancet::caller {

// ============================================================================
// Multi-allelic constructor
//
// Unspools a multi-allelic RawVariant into the VCF array mapping.
// ============================================================================
VariantCall::VariantCall(RawVariant const* var, SupportsByVariant const& all_supports,
                         Samples samps, PerSampleCov per_sample_cov)
    : mVariantId(HashRawVariant(var)),
      mChromIndex(var->mChromIndex),
      mStartPos1(var->mGenomeChromPos1),

      mRawVariant(var),
      mPerSampleCov(std::move(per_sample_cov)),
      mChromName(var->mChromName),
      mRefAllele(var->mRefAllele),
      mGraphCx(var->mGraphMetrics),
      mSeqCx(var->mSeqCx),
      mIsMultiallelic(mAltAlleles.size() > 1) {
  auto const num_alts = var->mAlts.size();
  mAltAlleles.reserve(num_alts);
  mCategories.reserve(num_alts);
  mVariantLengths.reserve(num_alts);

  std::ranges::for_each(var->mAlts, [this](auto const& alt) -> void {
    mAltAlleles.push_back(alt.mSequence);
    mCategories.push_back(alt.mType);
    mVariantLengths.push_back(alt.mLength);
  });

  if (auto const var_it = all_supports.find(var); var_it != all_supports.end()) {
    Finalize(var_it->second, samps);
  } else {
    Finalize(SupportArray(), samps);
  }
}

// ============================================================================
// Finalize: common post-construction steps for both constructors.
//
// Determines tumor-normal (somatic) vs normal-only mode, then delegates to
// three focused methods:
//   1. BuildFormatFields  — per-sample FORMAT strings + site quality
//   2. ComputeState       — SHARED/NORMAL/TUMOR classification
//   3. BuildInfoField     — INFO string assembly
// ============================================================================
void VariantCall::Finalize(SupportArray const& evidence, Samples samps) {
  static auto const IS_TUMOR = [](auto const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::TUMOR;
  };
  static auto const IS_NORMAL = [](auto const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::NORMAL;
  };
  auto const tumor_normal_mode =
      std::ranges::any_of(samps, IS_TUMOR) && std::ranges::any_of(samps, IS_NORMAL);

  BuildFormatFields(evidence, samps, tumor_normal_mode);
  ComputeState(evidence, samps, tumor_normal_mode);
  BuildInfoField(tumor_normal_mode);
}

// ============================================================================
// BuildFormatFields: per-sample FORMAT strings and site quality.
//
// clang-format off
// FORMAT: GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ
// clang-format on
//
//   GT    - Genotype derived from minimum PL
//   AD    - Number=R: allele depths (REF, ALT1, ALT2, ...)
//   ADF   - Number=R: forward strand allele depths
//   ADR   - Number=R: reverse strand allele depths
//   DP    - Total depth
//   RMQ   - Number=R: RMS mapping quality per allele
//   NPBQ  - Number=R: normalized posterior base quality per allele (PBQ/N)
//   SB    - Number=1: Strand bias log odds ratio (Haldane-corrected)
//   SCA   - Number=1: Soft Clip Asymmetry (ALT - REF soft-clip fraction)
//   FLD   - Number=1: Fragment Length Delta (mean ALT isize − mean REF isize, signed)
//   RPCD  - Number=1: Read Position Cohen's D (folded position effect size)
//   BQCD  - Number=1: Base Quality Cohen's D (base quality effect size)
//   MQCD  - Number=1: Mapping Quality Cohen's D (MAPQ effect size)
//   ASMD  - Number=1: Allele-Specific Mismatch Delta (mean ALT NM − mean REF NM − variant_length)
//   SDFC  - Number=1: Site Depth Fold Change (sample DP / per-sample window mean coverage)
//   PRAD  - Number=1: Polar Radius log10(1 + sqrt(AD_Ref² + AD_Alt²))
//   PANG  - Number=1: Polar Angle atan2(AD_Alt, AD_Ref) in radians
//   CMLOD - Number=A: Continuous Mixture LOD per ALT (quality-weighted)
//   FSSE  - Number=1: Fragment Start Shannon Entropy [0,1] (ALT start position diversity)
//   AHDD  - Number=1: ALT-Haplotype Discordance Delta (ALT reads vs own haplotype)
//   HSE   - Number=1: Haplotype Segregation Entropy [0,1] (ALT path concentration)
//   PDCV  - Number=1: Path Depth Coefficient of Variation (graph coverage uniformity)
//   PL    - Number=G: Phred-scaled genotype likelihoods (Dirichlet-Multinomial)
//   GQ    - Genotype quality (second-lowest PL from Dirichlet-Multinomial, capped at 99)
// ============================================================================
void VariantCall::BuildFormatFields(SupportArray const& evidence, Samples samps,
                                    bool const tumor_normal_mode) {
  auto const num_alleles = mAltAlleles.size() + 1;  // +1 for REF

  mSampleGenotypes.reserve(samps.size());
  for (auto const& sinfo : samps) {
    auto const* support = evidence.Find(sinfo.SampleName());
    if (support == nullptr) {
      mSampleGenotypes.emplace_back(SampleGenotypeData{.mIsMissingSupport = true});
      continue;
    }

    mTotalSampleCov += support->TotalSampleCov();
    auto const pls = support->ComputePLs();
    SampleGenotypeData sample;
    sample.mIsMissingSupport = false;

    AssignGenotype(sample, pls, num_alleles);
    AssignPerAlleleMetrics(sample, support, num_alleles);

    sample.mTotalDepth = support->TotalSampleCov();
    sample.mGenotypeQuality = VariantSupport::ComputeGQ(pls);

    UpdateSiteQuality(sinfo, support, evidence, samps, tumor_normal_mode, pls);
    mHasAltSupport = mHasAltSupport || (support->TotalAltCov() > 0);

    // Strand bias log odds ratio (Number=1, per-sample)
    sample.mStrandBias = static_cast<f32>(support->StrandBiasLogOR());

    // Alignment-derived per-sample annotations (coverage-normalized effect sizes).
    // ToOptF32 bridges std::optional<f64> (VariantSupport) → std::optional<f32> (VCF storage).
    sample.mSoftClipAsym = static_cast<f32>(support->SoftClipAsymmetry());
    sample.mFragLenDelta = ToOptF32(support->FragLengthDelta());
    sample.mReadPosCohenD = ToOptF32(support->ReadPosCohenD());
    sample.mBaseQualCohenD = ToOptF32(support->BaseQualCohenD());
    sample.mMapQualCohenD = ToOptF32(support->MappingQualCohenD());

    // ASMD: subtract max variant length so the variant's own edit distance
    // against REF doesn't inflate the mismatch delta (e.g., 50bp del = +50 NM).
    auto const max_var_len = std::transform_reduce(
        mVariantLengths.cbegin(), mVariantLengths.cend(), usize{0},
        [](usize acc, usize cur) { return std::max(acc, cur); },
        [](i64 len) { return static_cast<usize>(std::abs(len)); });
    sample.mAlleleMismatchDelta = ToOptF32(support->AlleleMismatchDelta(max_var_len));

    sample.mSiteDepthFoldChange =
        static_cast<f32>(SiteDepthFoldChange(sinfo.SampleName(), support->TotalSampleCov()));

    // Polar coordinate features for ML variant classification
    // PRAD/PANG separate allele identity from depth (see polar_coords.h)
    auto const ad_ref = static_cast<f64>(support->TotalRefCov());
    auto const ad_alt = static_cast<f64>(support->TotalAltCov());
    sample.mPolarRadius = static_cast<f32>(base::PolarRadius(ad_ref, ad_alt));
    sample.mPolarAngle = static_cast<f32>(base::PolarAngle(ad_alt, ad_ref));

    // Continuous Mixture LOD scores (CMLOD FORMAT field)
    auto const cmlod_scores = support->ComputeContinuousMixtureLods();
    sample.mContinuousMixtureLods.reserve(cmlod_scores.size());
    sample.mContinuousMixtureLods.insert(sample.mContinuousMixtureLods.end(), cmlod_scores.cbegin(),
                                         cmlod_scores.cend());

    // ── Artifact detection metrics ─────────────────────────────────────
    sample.mFragStartEntropy = ToOptF32(support->ComputeFSSE());
    sample.mAltHapDiscordDelta = ToOptF32(support->ComputeAHDD());

    if (mRawVariant != nullptr) {
      sample.mHaplotypeSegEntropy = ToOptF32(support->ComputeHSE(mRawVariant->mNumTotalHaps));
      sample.mPathDepthCv = ToOptF32(mRawVariant->mMaxPathCv);
    } else {
      sample.mHaplotypeSegEntropy = std::nullopt;
      sample.mPathDepthCv = std::nullopt;
    }

    sample.mPhredLikelihoods = pls;
    mSampleGenotypes.push_back(std::move(sample));
  }
}

// ============================================================================
// GenotypeFromGLIndex / AssignGenotype: convert a GL index back to a genotype tuple.
//
// VCF 4.3 §1.6.2 defines the GL index for a diploid genotype (i, jdx) as:
//
//   GL_index = jdx * (jdx + 1) / 2 + i     where i ≤ jdx
//
// This maps genotypes to a flat array indexed by triangular numbers:
//
//   Genotype:  0/0  0/1  1/1  0/2  1/2  2/2  0/3  1/3  2/3  3/3
//   GL index:   0    1    2    3    4    5    6    7    8    9
//                    └ jdx=1 ┘  └─ jdx=2 ─┘    └── jdx=3 ──┘
//
// Each "row" jdx starts at T(jdx) = jdx*(jdx+1)/2 and contains i = 0..jdx.
//
// To invert, we find which row jdx the GL index falls in, then recover i.
// We use the same integer-only algorithm as htslib's bcf_gt2alleles()
// (htslib/vcf.h): walk triangular numbers T(1), T(2), ... until T(dkmer) ≥ gl_index.
//
//   klen  = running triangular number T(dkmer-1)
//   dkmer = next row size (starts at 1, incremented each step)
//
//   After the loop: dkmer-1 = jdx, and i = gl_index - T(jdx) = gl_index - klen + jdx
//
// Example: gl_index = 4
//   dkmer=1, klen=0 → klen<4, dkmer=2, klen=2 → klen<4, dkmer=3, klen=5 → klen≥4, stop
//   jdx = dkmer-1 = 2,  i = 4 - 5 + 2 = 1  → "1/2" ✓
// ============================================================================
void VariantCall::AssignGenotype(SampleGenotypeData& sample, absl::Span<int const> pls,
                                 [[maybe_unused]] usize num_alleles) {
  if (pls.empty()) {
    sample.mGenotypeIndices = {-1, -1};
    return;
  }

  auto const best_gt_idx =
      static_cast<usize>(std::distance(pls.cbegin(), std::ranges::min_element(pls)));

  // Integer-only triangular number walk (matches htslib bcf_gt2alleles)
  usize klen = 0;
  usize dkmer = 1;
  while (klen < best_gt_idx) {
    dkmer++;
    klen += dkmer;
  }

  auto const jdx = dkmer - 1;
  auto const aidx = best_gt_idx - klen + jdx;
  // This should never trigger: gl_index comes from std::min_element over the
  // PL vector (size = num_alleles*(num_alleles+1)/2), so the inverted (i,j)
  // will always satisfy i <= j < num_alleles. If it does fire, ComputePLs or
  // its caller has a bug — don't silently return "0/0" and mask it.
  LANCET_ASSERT(aidx < num_alleles && jdx < num_alleles);
  sample.mGenotypeIndices = {static_cast<i16>(aidx), static_cast<i16>(jdx)};
}

// ============================================================================
// UpdateSiteQuality
//
// With DM-based PLs, the ref-hom PL naturally asymptotes at high depth
// (overdispersion absorbs correlated errors). The legacy `PL / SAMPLE_DP`
// normalization hack is fully deprecated — no artificial cap is needed.
//
// Germline mode: QUAL = ref-hom PL (PL[0/0]).
//   Directly measures confidence that the site is non-reference.
//   DM PLs asymptote via overdispersion — no artificial cap needed.
//
// Somatic mode:  QUAL = Somatic Log Odds Ratio (SOLOR).
//   Coverage-invariant tumor-vs-normal enrichment metric.
// ============================================================================
void VariantCall::UpdateSiteQuality(core::SampleInfo const& sinfo,
                                    [[maybe_unused]] VariantSupport const* support,
                                    SupportArray const& evidence, Samples samps,
                                    bool tumor_normal_mode, absl::Span<int const> pls) {
  if (tumor_normal_mode) {
    auto const somatic_lor = SomaticLogOddsRatio(sinfo, evidence, samps);
    mSiteQuality = std::max(mSiteQuality, somatic_lor);
  } else {
    // Ref-hom PL (index 0) = confidence the site is non-reference.
    // DM PLs asymptote via overdispersion — no artificial cap needed.
    auto const ref_hom_pl = pls.empty() ? 0 : pls[0];
    auto const germline_qual = static_cast<f64>(ref_hom_pl);
    mSiteQuality = std::max(mSiteQuality, germline_qual);
  }
}

// ============================================================================
// SomaticLogOddsRatio: somatic variant evidence as a coverage-invariant
// log odds ratio. Compares ALT/REF allele counts between the current tumor
// sample and the average across normal samples.
//
//   SOLOR = ln( ((tmr_alt+1)(nml_ref+1)) / ((tmr_ref+1)(nml_alt+1)) )
//
// Haldane correction (+1) handles zero-count edge cases without sentinels.
// A clean somatic produces SOLOR ≈ 5; germline produces SOLOR ≈ 0.
// Coverage-stable: once normal VAF stabilizes (≥60×), SOLOR varies < 3%.
// ============================================================================
auto VariantCall::SomaticLogOddsRatio(core::SampleInfo const& curr, SupportArray const& supports,
                                      Samples samps) -> f64 {
  if (curr.TagKind() != cbdg::Label::TUMOR) return 0.0;

  auto const* tmr_evidence = supports.Find(curr.SampleName());
  // Haldane correction (+1) mitigates undefined zero-division smoothly
  f64 const tmr_alt = tmr_evidence ? static_cast<f64>(tmr_evidence->TotalAltCov()) + 1.0 : 1.0;
  f64 const tmr_ref = tmr_evidence ? static_cast<f64>(tmr_evidence->TotalRefCov()) + 1.0 : 1.0;

  f64 sum_na = 0.0;
  f64 sum_nr = 0.0;
  f64 count_nml = 0.0;

  for (auto const& sinfo : samps) {
    auto const* evidence = supports.Find(sinfo.SampleName());
    if (sinfo.TagKind() != cbdg::Label::NORMAL || evidence == nullptr) continue;

    sum_na += static_cast<f64>(evidence->TotalAltCov());
    sum_nr += static_cast<f64>(evidence->TotalRefCov());
    count_nml += 1.0;
  }

  f64 const norm_val = std::max(count_nml, 1.0);
  f64 const nml_alt = (sum_na / norm_val) + 1.0;
  f64 const nml_ref = (sum_nr / norm_val) + 1.0;

  return std::log((tmr_alt * nml_ref) / (tmr_ref * nml_alt));
}

void VariantCall::AssignPerAlleleMetrics(SampleGenotypeData& sample, VariantSupport const* support,
                                         usize num_alleles) {
  sample.mAlleleDepths.reserve(num_alleles);
  sample.mFwdAlleleDepths.reserve(num_alleles);
  sample.mRevAlleleDepths.reserve(num_alleles);
  sample.mRmsMappingQualities.reserve(num_alleles);
  sample.mNormPosteriorBQs.reserve(num_alleles);

  for (usize allele = 0; allele < num_alleles; ++allele) {
    auto const idx = static_cast<AlleleIndex>(allele);
    sample.mAlleleDepths.push_back(support->TotalAlleleCov(idx));
    sample.mFwdAlleleDepths.push_back(support->FwdCount(idx));
    sample.mRevAlleleDepths.push_back(support->RevCount(idx));
    sample.mRmsMappingQualities.push_back(static_cast<f32>(support->RmsMappingQual(idx)));

    // NPBQ: raw posterior base quality divided by allele depth
    // Recovers the effective per-read quality (~30 for Q30 reads at any depth)
    auto const raw_pbq = support->RawPosteriorBaseQual(idx);
    auto const allele_cov = support->TotalAlleleCov(idx);
    auto const npbq = allele_cov > 0 ? raw_pbq / static_cast<f64>(allele_cov) : 0.0;
    sample.mNormPosteriorBQs.push_back(static_cast<f32>(npbq));
  }
}

// ============================================================================
// ComputeState: classify variant as SHARED, NORMAL, TUMOR, UNKNOWN, or NONE.
//
// In tumor-normal mode: SHARED/NORMAL/TUMOR based on ALT presence, NONE if no support.
// In normal-only mode: always UNKNOWN (not enough info to classify).
// ============================================================================
void VariantCall::ComputeState(SupportArray const& evidence, Samples samps,
                               bool const tumor_normal_mode) {
  if (!tumor_normal_mode) {
    mState = RawVariant::State::UNKNOWN;
    return;
  }

  auto const has_alt = [&evidence](core::SampleInfo const& sinfo, cbdg::Label::Tag kind) -> bool {
    auto const* support = evidence.Find(sinfo.SampleName());
    return sinfo.TagKind() == kind && support != nullptr && support->TotalAltCov() > 0;
  };

  bool const in_normal = std::ranges::any_of(
      samps, [&](auto const& sinfo) { return has_alt(sinfo, cbdg::Label::NORMAL); });

  bool const in_tumor = std::ranges::any_of(
      samps, [&](auto const& sinfo) { return has_alt(sinfo, cbdg::Label::TUMOR); });

  static constexpr std::array<RawVariant::State, 4> STATE_MAP = {{
      RawVariant::State::NONE,    // 00: Neither
      RawVariant::State::NORMAL,  // 01: Normal only
      RawVariant::State::TUMOR,   // 10: Tumor only
      RawVariant::State::SHARED   // 11: Both
  }};

  // We compute a 2-bit mapping index directly from the boolean parameters:
  // Bit 1 (Left bit):  in_tumor
  // Bit 0 (Right bit): in_normal
  // Evaluates exactly into [0=NONE, 1=NORMAL, 2=TUMOR, 3=SHARED].
  usize const state_idx = (static_cast<usize>(in_tumor) << 1) | static_cast<usize>(in_normal);
  mState = STATE_MAP[state_idx];
}

// ============================================================================
// BuildInfoField: assemble the VCF INFO string.
//
// Structure:  [STATE;]TYPE=<type>;LENGTH=<len>[;GRAPH_CX=...][;SEQ_CX=...]
//
//   STATE     — SHARED/NORMAL/TUMOR (tumor-normal mode only; omitted otherwise)
//   TYPE      — SNV, INS, DEL, MNP (always present)
//   LENGTH    — variant length in bp (always present)
//   GRAPH_CX  — graph complexity (GEI, TipToPathCovRatio, MaxDegree)
//   SEQ_CX    — sequence complexity (11 coverage-invariant features)
//
// Note: SCA, FLD, and MQCD are per-sample FORMAT fields, not site-level INFO.
// ============================================================================
void VariantCall::BuildInfoField(bool const tumor_normal_mode) {
  using namespace std::string_view_literals;

  static constexpr std::array<std::string_view, 6> TYPE_MAP = {
      {"REF"sv, "SNV"sv, "INS"sv, "DEL"sv, "MNP"sv, "CPX"sv}};

  std::vector<std::string_view> vcategories;
  vcategories.reserve(mCategories.size());
  for (auto const cat : mCategories) {
    vcategories.push_back(TYPE_MAP[static_cast<i8>(cat) + 1]);
  }

  std::string info;
  info.reserve(1024);

  if (tumor_normal_mode) {
    // State prefix — tumor-normal mode only
    // SHARED/NORMAL/TUMOR state — only in tumor-normal (somatic) mode
    static constexpr std::array<std::string_view, 5> STATE_MAP = {
        {"NONE"sv, "SHARED"sv, "NORMAL"sv, "TUMOR"sv, "UNKNOWN"sv}};

    absl::StrAppend(&info, STATE_MAP[static_cast<i8>(mState) + 1], ";");
  }

  if (mIsMultiallelic) absl::StrAppend(&info, "MULTIALLELIC;");

  absl::StrAppend(&info, "TYPE=", absl::StrJoin(vcategories, ","),
                  ";LENGTH=", absl::StrJoin(mVariantLengths, ","));

  absl::StrAppend(&info, ";GRAPH_CX=", mGraphCx.FormatVcfValue());
  absl::StrAppend(&info, ";SEQ_CX=", mSeqCx.FormatVcfValue());

  mInfoField = std::move(info);
}

// ============================================================================
// AsVcfRecord: emit a VCF record with comma-separated ALTs for multi-allelic.
// ============================================================================
auto VariantCall::AsVcfRecord() const -> std::string {
  auto const alt_field = absl::StrJoin(mAltAlleles, ",");
  std::vector<std::string> format_strings;
  format_strings.reserve(mSampleGenotypes.size() + 1);
  // clang-format off
  format_strings.emplace_back("GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ");
  // clang-format on

  for (auto const& sample : mSampleGenotypes) {
    format_strings.push_back(sample.RenderVcfString());
  }

  // clang-format off
  return fmt::format("{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL:.2f}\t.\t{INFO}\t{FORMAT}",
                     fmt::arg("CHROM", mChromName), fmt::arg("POS", mStartPos1),
                     fmt::arg("REF", mRefAllele), fmt::arg("ALT", alt_field),
                     fmt::arg("QUAL", mSiteQuality), fmt::arg("INFO", mInfoField),
                     fmt::arg("FORMAT", absl::StrJoin(format_strings, "\t")));
  // clang-format on
}

auto VariantCall::SampleGenotypeData::RenderVcfString() const -> std::string {
  if (mIsMissingSupport) {
    return "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.";
  }

  auto const gt_str = (mGenotypeIndices.first == -1 && mGenotypeIndices.second == -1)
                          ? std::string("./.")
                          : fmt::format("{}/{}", mGenotypeIndices.first, mGenotypeIndices.second);

  auto const pl_str = mPhredLikelihoods.empty() ? "." : absl::StrJoin(mPhredLikelihoods, ",");

  auto const format_f32 = [](std::string* out, f32 const val) {
    absl::StrAppendFormat(out, "%.1F", val);
  };

  auto const format_f64 = [](std::string* out, f64 const val) {
    absl::StrAppendFormat(out, "%.4F", val);
  };

  auto const ad_str = absl::StrJoin(mAlleleDepths, ",");
  auto const adf_str = absl::StrJoin(mFwdAlleleDepths, ",");
  auto const adr_str = absl::StrJoin(mRevAlleleDepths, ",");
  auto const rmq_str = absl::StrJoin(mRmsMappingQualities, ",", format_f32);
  auto const npbq_str = absl::StrJoin(mNormPosteriorBQs, ",", format_f32);

  // CMLOD: Number=A — strip index 0 (REF LOD = 0.0 by definition), emit only ALT LODs
  auto const cmlod_str = mContinuousMixtureLods.empty()
                             ? std::string(".")
                             : absl::StrJoin(mContinuousMixtureLods.begin() + 1,
                                             mContinuousMixtureLods.end(), ",", format_f64);

  // Optional-safe formatting: 9 metrics can be std::nullopt (untestable) → "." in VCF.
  // All other fields are always populated (not optional-capable).
  auto const fld_str = FormatOptional(mFragLenDelta, "{:.1f}");
  auto const rpcd_str = FormatOptional(mReadPosCohenD, "{:.4f}");
  auto const bqcd_str = FormatOptional(mBaseQualCohenD, "{:.4f}");
  auto const mqcd_str = FormatOptional(mMapQualCohenD, "{:.4f}");
  auto const asmd_str = FormatOptional(mAlleleMismatchDelta, "{:.3f}");
  auto const fsse_str = FormatOptional(mFragStartEntropy, "{:.4f}");
  auto const ahdd_str = FormatOptional(mAltHapDiscordDelta, "{:.3f}");
  auto const hse_str = FormatOptional(mHaplotypeSegEntropy, "{:.4f}");
  auto const pdcv_str = FormatOptional(mPathDepthCv, "{:.4f}");

  // clang-format off
  return fmt::format(
      "{GT}:{AD}:{ADF}:{ADR}:{DP}:{RMQ}:{NPBQ}:{SB:.3f}:{SCA:.4f}:{FLD}:{RPCD}:"
      "{BQCD}:{MQCD}:{ASMD}:{SDFC:.2f}:{PRAD:.4f}:{PANG:.4f}:{CMLOD}:"
      "{FSSE}:{AHDD}:{HSE}:{PDCV}:{PL}:{GQ}",
      fmt::arg("GT",    gt_str),           fmt::arg("AD",    ad_str),
      fmt::arg("ADF",   adf_str),          fmt::arg("ADR",   adr_str),
      fmt::arg("DP",    mTotalDepth),      fmt::arg("RMQ",   rmq_str),
      fmt::arg("NPBQ",  npbq_str),         fmt::arg("SB",    mStrandBias),
      fmt::arg("SCA",   mSoftClipAsym),    fmt::arg("FLD",   fld_str),
      fmt::arg("RPCD",  rpcd_str),         fmt::arg("BQCD",  bqcd_str),
      fmt::arg("MQCD",  mqcd_str),         fmt::arg("ASMD",  asmd_str),
      fmt::arg("SDFC",  mSiteDepthFoldChange), fmt::arg("PRAD",  mPolarRadius),
      fmt::arg("PANG",  mPolarAngle),      fmt::arg("CMLOD", cmlod_str),
      fmt::arg("FSSE",  fsse_str),         fmt::arg("AHDD",  ahdd_str),
      fmt::arg("HSE",   hse_str),          fmt::arg("PDCV",  pdcv_str),
      fmt::arg("PL",    pl_str),           fmt::arg("GQ",    mGenotypeQuality));
  // clang-format on
}

}  // namespace lancet::caller
