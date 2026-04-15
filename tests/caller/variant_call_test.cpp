#include "lancet/caller/variant_call.h"

#include "lancet/caller/sample_format_data.h"

#include "catch_amalgamated.hpp"

namespace lancet::caller::tests {

// FORMAT field order:
//   GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ
TEST_CASE("SampleFormatData compactly serializes VCF genotype structure strings precisely",
          "[lancet][caller][VariantCall]") {
  SampleFormatData sample;
  sample.SetPhredLikelihoods({10, 0, 100});
  sample.SetAlleleDepths({30, 20});
  sample.SetFwdAlleleDepths({15, 10});
  sample.SetRevAlleleDepths({15, 10});
  sample.SetRmsMappingQualities({59.5F, 60.0F});
  sample.SetNormPosteriorBQs({30.1F, 35.2F});
  sample.SetContinuousMixtureLods({0.0, 12.5432});
  sample.SetStrandBias(2.45F);
  sample.SetSoftClipAsym(0.5F);
  sample.SetField(SampleFormatData::FRAG_LEN_DELTA, 5.2);
  sample.SetField(SampleFormatData::READ_POS_COHEN_D, 0.1);
  sample.SetField(SampleFormatData::BASE_QUAL_COHEN_D, 0.2);
  sample.SetField(SampleFormatData::MAP_QUAL_COHEN_D, -0.5);
  sample.SetField(SampleFormatData::ALLELE_MISMATCH_DELTA, 0.05);
  sample.SetSiteDepthFoldChange(1.25F);
  sample.SetPolarRadius(5.51F);
  sample.SetPolarAngle(0.588F);
  sample.SetTotalDepth(50);
  sample.SetGenotypeQuality(10);
  sample.SetGenotypeIndices(0, 1);
  sample.SetMissingSupport(false);

  auto const vcf_str = sample.RenderVcfString();
  REQUIRE(vcf_str == "0/1:30,20:15,10:15,10:50:59.5,60.0:30.1,35.2:2.450:0.5000:5.2:"
                     "0.1000:0.2000:-0.5000:0.050:1.25:5.5100:0.5880:12.5432:"
                     ".:.:.:.:10,0,100:10");
  sample.SetMissingSupport(true);
  REQUIRE(sample.RenderVcfString() == "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.");
}

TEST_CASE("SampleFormatData renders nullopt metrics as VCF missing value dot",
          "[lancet][caller][VariantCall]") {
  SampleFormatData sample;
  sample.SetPhredLikelihoods({0, 10, 50});
  sample.SetAlleleDepths({10, 0});
  sample.SetFwdAlleleDepths({5, 0});
  sample.SetRevAlleleDepths({5, 0});
  sample.SetRmsMappingQualities({60.0F, 0.0F});
  sample.SetNormPosteriorBQs({30.0F, 0.0F});
  sample.SetContinuousMixtureLods({0.0, 0.0});
  sample.SetStrandBias(0.0F);
  sample.SetSoftClipAsym(0.0F);
  // Leave 9 optional metrics at their default (absent) — untestable → "." in VCF
  sample.SetSiteDepthFoldChange(1.0F);
  sample.SetPolarRadius(1.0F);
  sample.SetPolarAngle(0.0F);
  sample.SetTotalDepth(10);
  sample.SetGenotypeQuality(0);
  sample.SetGenotypeIndices(0, 0);
  sample.SetMissingSupport(false);

  auto const vcf_str = sample.RenderVcfString();
  // FLD, RPCD, BQCD, MQCD, ASMD, FSSE, AHDD, HSE, PDCV should all render as "." (absent → dot)
  REQUIRE(vcf_str == "0/0:10,0:5,0:5,0:10:60.0,0.0:30.0,0.0:0.000:0.0000:.:"
                     ".:.:.:.:1.00:1.0000:0.0000:0.0000:"
                     ".:.:.:.:0,10,50:0");
}

}  // namespace lancet::caller::tests
