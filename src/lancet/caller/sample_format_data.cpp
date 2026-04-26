#include "lancet/caller/sample_format_data.h"

#include "lancet/base/types.h"

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"

#include <optional>
#include <string>

namespace {

// Format a std::optional<f32> for VCF output, emitting "." for missing values.
// VCF 4.5 spec: "." = missing value. This is the single point where the
// optional → dot conversion happens for all FORMAT fields.
// NOTE: Cannot use IEEE NaN here because -ffast-math (cmake/defaults.cmake)
// implies -ffinite-math-only, which optimizes away std::isnan().
[[nodiscard]] inline auto FormatOptional(std::optional<f32> const& value, char const* fmt_spec)
    -> std::string {
  if (!value.has_value()) return ".";
  return fmt::format(fmt::runtime(fmt_spec), value.value());
}

}  // namespace

namespace lancet::caller {

auto SampleFormatData::RenderVcfString() const -> std::string {
  if (IsMissingSupport()) {
    return "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.";
  }

  auto const [gt_first, gt_second] = GenotypeIndices();
  auto const gt_str = (gt_first == -1 && gt_second == -1)
                          ? std::string("./.")
                          : fmt::format("{}/{}", gt_first, gt_second);

  auto const pls = PhredLikelihoods();
  auto const pl_str = pls.empty() ? "." : absl::StrJoin(pls, ",");

  auto const format_f32 = [](std::string* out, f32 const val) {
    absl::StrAppendFormat(out, "%.1F", val);
  };

  auto const format_f64 = [](std::string* out, f64 const val) {
    absl::StrAppendFormat(out, "%.4F", val);
  };

  auto const ad_str = absl::StrJoin(AlleleDepths(), ",");
  auto const adf_str = absl::StrJoin(FwdAlleleDepths(), ",");
  auto const adr_str = absl::StrJoin(RevAlleleDepths(), ",");
  auto const rmq_str = absl::StrJoin(RmsMappingQualities(), ",", format_f32);
  auto const npbq_str = absl::StrJoin(NormPosteriorBQs(), ",", format_f32);

  // CMLOD: Number=A — strip index 0 (REF LOD = 0.0 by definition), emit only ALT LODs.
  // Guard: size() < 2 catches both empty and REF-only vectors (defensive, upstream
  // should always produce num_alleles entries after the allele count fix).
  auto const cmlods = ContinuousMixtureLods();
  auto const cmlod_str = [&cmlods, &format_f64]() -> std::string {
    if (cmlods.size() < 2) return ".";
    return absl::StrJoin(cmlods.begin() + 1, cmlods.end(), ",", format_f64);
  }();

  // Optional-safe formatting via GetField: 9 metrics can be absent → "." in VCF.
  auto const fld_str = FormatOptional(GetField(FRAG_LEN_DELTA), "{:.1f}");
  auto const rpcd_str = FormatOptional(GetField(READ_POS_COHEN_D), "{:.4f}");
  auto const bqcd_str = FormatOptional(GetField(BASE_QUAL_COHEN_D), "{:.4f}");
  auto const mqcd_str = FormatOptional(GetField(MAP_QUAL_COHEN_D), "{:.4f}");
  auto const asmd_str = FormatOptional(GetField(ALLELE_MISMATCH_DELTA), "{:.3f}");
  auto const fsse_str = FormatOptional(GetField(FRAG_START_ENTROPY), "{:.4f}");
  auto const ahdd_str = FormatOptional(GetField(ALT_HAP_DISCORD_DELTA), "{:.3f}");
  auto const hse_str = FormatOptional(GetField(HAPLOTYPE_SEG_ENTROPY), "{:.4f}");
  auto const pdcv_str = FormatOptional(GetField(PATH_DEPTH_CV), "{:.4f}");
  auto const sdfc_str = FormatOptional(GetField(SITE_DEPTH_FOLD_CHANGE), "{:.2f}");

  // clang-format off
  return fmt::format(
      "{GT}:{AD}:{ADF}:{ADR}:{DP}:{RMQ}:{NPBQ}:{SB:.3f}:{SCA:.4f}:{FLD}:{RPCD}:"
      "{BQCD}:{MQCD}:{ASMD}:{SDFC}:{PRAD:.4f}:{PANG:.4f}:{CMLOD}:"
      "{FSSE}:{AHDD}:{HSE}:{PDCV}:{PL}:{GQ}",
      fmt::arg("GT",    gt_str),           fmt::arg("AD",    ad_str),
      fmt::arg("ADF",   adf_str),          fmt::arg("ADR",   adr_str),
      fmt::arg("DP",    TotalDepth()),      fmt::arg("RMQ",   rmq_str),
      fmt::arg("NPBQ",  npbq_str),         fmt::arg("SB",    StrandBias()),
      fmt::arg("SCA",   SoftClipAsym()),    fmt::arg("FLD",   fld_str),
      fmt::arg("RPCD",  rpcd_str),         fmt::arg("BQCD",  bqcd_str),
      fmt::arg("MQCD",  mqcd_str),         fmt::arg("ASMD",  asmd_str),
      fmt::arg("SDFC",  sdfc_str),         fmt::arg("PRAD",  PolarRadius()),
      fmt::arg("PANG",  PolarAngle()),      fmt::arg("CMLOD", cmlod_str),
      fmt::arg("FSSE",  fsse_str),         fmt::arg("AHDD",  ahdd_str),
      fmt::arg("HSE",   hse_str),          fmt::arg("PDCV",  pdcv_str),
      fmt::arg("PL",    pl_str),           fmt::arg("GQ",    GenotypeQuality()));
  // clang-format on
}

}  // namespace lancet::caller
