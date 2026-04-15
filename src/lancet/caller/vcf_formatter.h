#ifndef SRC_LANCET_CALLER_VCF_FORMATTER_H_
#define SRC_LANCET_CALLER_VCF_FORMATTER_H_

#include <string_view>

namespace lancet::caller {

// ============================================================================
// FORMAT_HEADER: authoritative single-source-of-truth for VCF FORMAT fields.
//
// All FORMAT field rendering (AsVcfRecord, RenderVcfString, missing-support
// strings) must be kept in sync with this constant. Changes here must be
// reflected in sample_format_data.cpp's RenderVcfString format template.
//
// clang-format off
// Field ordering rationale (VCF convention: genotype fields first, then
// per-allele depths, then site-level metrics, then model outputs):
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
//   FLD  - Number=1: Fragment Length Delta (mean ALT isize − mean REF isize)
//   RPCD - Number=1: Read Position Cohen's D (folded position effect size)
//   BQCD - Number=1: Base Quality Cohen's D (base quality effect size)
//   MQCD  - Number=1: Mapping Quality Cohen's D (MAPQ effect size)
//   ASMD  - Number=1: Allele-Specific Mismatch Delta
//   SDFC  - Number=1: Site Depth Fold Change
//   PRAD  - Number=1: Polar Radius log10(1 + sqrt(AD_Ref² + AD_Alt²))
//   PANG  - Number=1: Polar Angle atan2(AD_Alt, AD_Ref) in radians
//   CMLOD - Number=A: Continuous Mixture LOD per ALT
//   FSSE  - Number=1: Fragment Start Shannon Entropy [0,1]
//   AHDD  - Number=1: ALT-Haplotype Discordance Delta
//   HSE   - Number=1: Haplotype Segregation Entropy [0,1]
//   PDCV  - Number=1: Path Depth Coefficient of Variation
//   PL    - Number=G: Phred-scaled genotype likelihoods (Dirichlet-Multinomial)
//   GQ    - Genotype quality (second-lowest DM PL, capped at 99)
// ============================================================================
inline constexpr std::string_view FORMAT_HEADER = 
    "GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ";
// clang-format on

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VCF_FORMATTER_H_
