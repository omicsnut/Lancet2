#include "lancet/caller/variant_call.h"

#include "lancet/base/assert.h"
#include "lancet/base/polar_coords.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_support.h"
#include "lancet/caller/vcf_formatter.h"
#include "lancet/cbdg/label.h"

#include "absl/hash/hash.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"

#include <absl/container/flat_hash_map.h>
#include <absl/container/inlined_vector.h>
#include <algorithm>
#include <array>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>

namespace {

// Variant identity is defined by CHROM+POS+REF (locus-level, no ALTs).
// Uses mChromIndex (integer, always correctly set) instead of mChromName for robustness.
// See the identity design note on VariantCall::operator< for full rationale.
[[nodiscard]] inline auto HashRawVariant(lancet::caller::RawVariant const* var) -> u64 {
  return static_cast<u64>(absl::HashOf(var->mChromIndex, var->mGenomeChromPos1, var->mRefAllele));
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
      mIsMultiallelic(var->mAlts.size() > 1) {
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
// Determines case-control vs control-only mode, then delegates to
// three focused methods:
//   1. BuildFormatFields  — per-sample FORMAT strings + site quality
//   2. ComputeState       — SHARED/CTRL/CASE classification
//   3. BuildInfoField     — INFO string assembly
// ============================================================================
void VariantCall::Finalize(SupportArray const& evidence, Samples samps) {
  // Detect case-control mode by checking for at least one sample of each role.
  // TagKind serves as the role identifier: Label::CASE = focal, Label::CTRL = baseline.
  static auto const IS_CASE = [](auto const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::CASE;
  };
  static auto const IS_CTRL = [](auto const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::CTRL;
  };
  auto const case_ctrl_mode =
      std::ranges::any_of(samps, IS_CASE) && std::ranges::any_of(samps, IS_CTRL);

  BuildFormatFields(evidence, samps, case_ctrl_mode);
  ComputeState(evidence, samps, case_ctrl_mode);
  BuildInfoField(case_ctrl_mode);
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
                                    bool const case_ctrl_mode) {
  auto const num_alleles = mAltAlleles.size() + 1;  // +1 for REF

  mSampleGenotypes.reserve(samps.size());
  for (auto const& sinfo : samps) {
    auto const* support = evidence.Find(sinfo.SampleName());
    if (support == nullptr) {
      SampleFormatData missing;
      missing.SetMissingSupport(true);
      mSampleGenotypes.push_back(std::move(missing));
      continue;
    }

    mTotalSampleCov += support->TotalSampleCov();
    auto const pls = support->ComputePLs();
    SampleFormatData sample;
    sample.SetMissingSupport(false);

    AssignGenotype(sample, pls, num_alleles);
    AssignPerAlleleMetrics(sample, support, num_alleles);

    sample.SetTotalDepth(support->TotalSampleCov());
    sample.SetGenotypeQuality(VariantSupport::ComputeGQ(pls));

    UpdateSiteQuality(sinfo, support, evidence, samps, case_ctrl_mode, pls);
    mHasAltSupport = mHasAltSupport || (support->TotalAltCov() > 0);

    // Strand bias log odds ratio (Number=1, per-sample)
    sample.SetStrandBias(static_cast<f32>(support->StrandBiasLogOR()));

    // Alignment-derived per-sample annotations (coverage-normalized effect sizes).
    // SetField bridges std::optional<f64> (VariantSupport) → compact f32 storage directly.
    sample.SetSoftClipAsym(static_cast<f32>(support->SoftClipAsymmetry()));
    sample.SetField(SampleFormatData::FRAG_LEN_DELTA, support->FragLengthDelta());
    sample.SetField(SampleFormatData::READ_POS_COHEN_D, support->ReadPosCohenD());
    sample.SetField(SampleFormatData::BASE_QUAL_COHEN_D, support->BaseQualCohenD());
    sample.SetField(SampleFormatData::MAP_QUAL_COHEN_D, support->MappingQualCohenD());

    // ASMD: subtract max variant length so the variant's own edit distance
    // against REF doesn't inflate the mismatch delta (e.g., 50bp del = +50 NM).
    auto const max_var_len = std::transform_reduce(
        mVariantLengths.cbegin(), mVariantLengths.cend(), usize{0},
        [](usize acc, usize cur) { return std::max(acc, cur); },
        [](i64 len) { return static_cast<usize>(std::abs(len)); });
    sample.SetField(SampleFormatData::ALLELE_MISMATCH_DELTA,
                    support->AlleleMismatchDelta(max_var_len));

    sample.SetSiteDepthFoldChange(
        static_cast<f32>(SiteDepthFoldChange(sinfo.SampleName(), support->TotalSampleCov())));

    // Polar coordinate features for ML variant classification
    // PRAD/PANG separate allele identity from depth (see polar_coords.h)
    auto const ad_ref = static_cast<f64>(support->TotalRefCov());
    auto const ad_alt = static_cast<f64>(support->TotalAltCov());
    sample.SetPolarRadius(static_cast<f32>(base::PolarRadius(ad_ref, ad_alt)));
    sample.SetPolarAngle(static_cast<f32>(base::PolarAngle(ad_alt, ad_ref)));

    // Continuous Mixture LOD scores (CMLOD FORMAT field)
    auto const cmlod_scores = support->ComputeContinuousMixtureLods();
    absl::InlinedVector<f64, 4> cmlod_vec(cmlod_scores.cbegin(), cmlod_scores.cend());
    sample.SetContinuousMixtureLods(std::move(cmlod_vec));

    // ── Artifact detection metrics ─────────────────────────────────────
    sample.SetField(SampleFormatData::FRAG_START_ENTROPY, support->ComputeFSSE());
    sample.SetField(SampleFormatData::ALT_HAP_DISCORD_DELTA, support->ComputeAHDD());

    if (mRawVariant != nullptr) {
      sample.SetField(SampleFormatData::HAPLOTYPE_SEG_ENTROPY,
                      support->ComputeHSE(mRawVariant->mNumTotalHaps));
      // Bridge optional<f64> → SetField for PDCV
      sample.SetField(SampleFormatData::PATH_DEPTH_CV, mRawVariant->mMaxPathCv);
    }

    sample.SetPhredLikelihoods(pls);
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
void VariantCall::AssignGenotype(SampleFormatData& sample, absl::Span<u32 const> pls,
                                 [[maybe_unused]] usize num_alleles) {
  if (pls.empty()) {
    sample.SetGenotypeIndices(-1, -1);
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
  sample.SetGenotypeIndices(static_cast<i16>(aidx), static_cast<i16>(jdx));
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
//   Coverage-invariant case-vs-control enrichment metric.
// ============================================================================
void VariantCall::UpdateSiteQuality(core::SampleInfo const& sinfo,
                                    [[maybe_unused]] VariantSupport const* support,
                                    SupportArray const& evidence, Samples samps,
                                    bool case_ctrl_mode, absl::Span<u32 const> pls) {
  if (case_ctrl_mode) {
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
// log odds ratio. Compares ALT/REF allele counts between the current case
// sample and the average across control samples.
//
//   SOLOR = ln( ((case_alt+1)(ctrl_ref+1)) / ((case_ref+1)(ctrl_alt+1)) )
//
// Haldane correction (+1) handles zero-count edge cases without sentinels.
// A clean somatic produces SOLOR ≈ 5; germline (shared) produces SOLOR ≈ 0.
// Coverage-stable: once control VAF stabilizes (≥60×), SOLOR varies < 3%.
// ============================================================================
auto VariantCall::SomaticLogOddsRatio(core::SampleInfo const& curr, SupportArray const& supports,
                                      Samples samps) -> f64 {
  // SOLOR is only meaningful for focal (case) samples compared against baseline (control).
  if (curr.TagKind() != cbdg::Label::CASE) return 0.0;

  auto const* case_evidence = supports.Find(curr.SampleName());
  // Haldane correction (+1) mitigates undefined zero-division smoothly
  f64 const case_alt = case_evidence ? static_cast<f64>(case_evidence->TotalAltCov()) + 1.0 : 1.0;
  f64 const case_ref = case_evidence ? static_cast<f64>(case_evidence->TotalRefCov()) + 1.0 : 1.0;

  f64 sum_ctrl_alt = 0.0;
  f64 sum_ctrl_ref = 0.0;
  f64 count_ctrl = 0.0;

  for (auto const& sinfo : samps) {
    auto const* evidence = supports.Find(sinfo.SampleName());
    // Accumulate only baseline (control) samples for the denominator.
    if (sinfo.TagKind() != cbdg::Label::CTRL || evidence == nullptr) continue;

    sum_ctrl_alt += static_cast<f64>(evidence->TotalAltCov());
    sum_ctrl_ref += static_cast<f64>(evidence->TotalRefCov());
    count_ctrl += 1.0;
  }

  f64 const ctrl_count = std::max(count_ctrl, 1.0);
  f64 const ctrl_alt = (sum_ctrl_alt / ctrl_count) + 1.0;
  f64 const ctrl_ref = (sum_ctrl_ref / ctrl_count) + 1.0;

  return std::log((case_alt * ctrl_ref) / (case_ref * ctrl_alt));
}

void VariantCall::AssignPerAlleleMetrics(SampleFormatData& sample, VariantSupport const* support,
                                         usize num_alleles) {
  absl::InlinedVector<u16, 4> allele_depths;
  absl::InlinedVector<u16, 4> fwd_depths;
  absl::InlinedVector<u16, 4> rev_depths;
  absl::InlinedVector<f32, 4> rms_mapq;
  absl::InlinedVector<f32, 4> norm_pbq;

  allele_depths.reserve(num_alleles);
  fwd_depths.reserve(num_alleles);
  rev_depths.reserve(num_alleles);
  rms_mapq.reserve(num_alleles);
  norm_pbq.reserve(num_alleles);

  for (usize allele = 0; allele < num_alleles; ++allele) {
    auto const idx = static_cast<AlleleIndex>(allele);
    allele_depths.push_back(support->TotalAlleleCov(idx));
    fwd_depths.push_back(support->FwdCount(idx));
    rev_depths.push_back(support->RevCount(idx));
    rms_mapq.push_back(static_cast<f32>(support->RmsMappingQual(idx)));

    // NPBQ: raw posterior base quality divided by allele depth
    // Recovers the effective per-read quality (~30 for Q30 reads at any depth)
    auto const raw_pbq = support->RawPosteriorBaseQual(idx);
    auto const allele_cov = support->TotalAlleleCov(idx);
    auto const npbq = allele_cov > 0 ? raw_pbq / static_cast<f64>(allele_cov) : 0.0;
    norm_pbq.push_back(static_cast<f32>(npbq));
  }

  sample.SetAlleleDepths(std::move(allele_depths));
  sample.SetFwdAlleleDepths(std::move(fwd_depths));
  sample.SetRevAlleleDepths(std::move(rev_depths));
  sample.SetRmsMappingQualities(std::move(rms_mapq));
  sample.SetNormPosteriorBQs(std::move(norm_pbq));
}

// ============================================================================
// ComputeState: classify variant as SHARED, CTRL, CASE, UNKNOWN, or NONE.
//
// In case-control mode: SHARED/CTRL/CASE based on ALT presence, NONE if no support.
// In control-only mode: always UNKNOWN (not enough info to classify).
// ============================================================================
void VariantCall::ComputeState(SupportArray const& evidence, Samples samps,
                               bool const case_ctrl_mode) {
  if (!case_ctrl_mode) {
    mState = AlleleState::UNKNOWN;
    return;
  }

  auto const has_alt = [&evidence](core::SampleInfo const& sinfo, cbdg::Label::Tag role) -> bool {
    auto const* support = evidence.Find(sinfo.SampleName());
    return sinfo.TagKind() == role && support != nullptr && support->TotalAltCov() > 0;
  };

  bool const in_ctrl = std::ranges::any_of(
      samps, [&](auto const& sinfo) { return has_alt(sinfo, cbdg::Label::CTRL); });

  bool const in_case = std::ranges::any_of(
      samps, [&](auto const& sinfo) { return has_alt(sinfo, cbdg::Label::CASE); });

  static constexpr std::array<AlleleState, 4> STATE_MAP = {{
      AlleleState::NONE,   // 00: Neither control nor case
      AlleleState::CTRL,   // 01: Control only
      AlleleState::CASE,   // 10: Case only
      AlleleState::SHARED  // 11: Both control and case
  }};

  // 2-bit mapping index from boolean parameters:
  // Bit 1 (left):  in_case
  // Bit 0 (right): in_ctrl
  // Evaluates into [0=NONE, 1=CTRL, 2=CASE, 3=SHARED].
  usize const state_idx = (static_cast<usize>(in_case) << 1) | static_cast<usize>(in_ctrl);
  mState = STATE_MAP[state_idx];
}

// ============================================================================
// BuildInfoField: assemble the VCF INFO string.
//
// Structure:  [STATE;]TYPE=<type>;LENGTH=<len>[;GRAPH_CX=...][;SEQ_CX=...]
//
//   STATE     — SHARED/CTRL/CASE (case-control mode only; omitted otherwise)
//   TYPE      — SNV, INS, DEL, MNP (always present)
//   LENGTH    — variant length in bp (always present)
//   GRAPH_CX  — graph complexity (GEI, TipToPathCovRatio, MaxDegree)
//   SEQ_CX    — sequence complexity (11 coverage-invariant features)
//
// Note: SCA, FLD, and MQCD are per-sample FORMAT fields, not site-level INFO.
// ============================================================================
void VariantCall::BuildInfoField(bool const case_ctrl_mode) {
  using namespace std::string_view_literals;

  static constexpr std::array<std::string_view, 6> TYPE_MAP = {
      {"REF"sv, "SNV"sv, "INS"sv, "DEL"sv, "MNP"sv, "CPX"sv}};

  std::vector<std::string_view> vcategories(mCategories.size());
  std::ranges::transform(mCategories, vcategories.begin(),
                         [](auto categ) { return TYPE_MAP[static_cast<i8>(categ) + 1]; });

  std::string info;
  info.reserve(1024);

  if (case_ctrl_mode) {
    // SHARED/CTRL/CASE state prefix — case-control mode only
    static constexpr std::array<std::string_view, 5> STATE_MAP = {
        {"NONE"sv, "SHARED"sv, "CTRL"sv, "CASE"sv, "UNKNOWN"sv}};

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
  format_strings.emplace_back(FORMAT_HEADER);

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
}  // namespace lancet::caller
