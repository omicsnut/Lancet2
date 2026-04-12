#include "lancet/caller/variant_call.h"

#include "catch_amalgamated.hpp"

#include <optional>
#include <vector>

namespace lancet::caller::tests {

// FORMAT field order:
//   GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:PL:GQ
TEST_CASE("SampleGenotypeData compactly serializes VCF genotype structure strings precisely",
          "[lancet][caller][VariantCall]") {
  VariantCall::SampleGenotypeData sample;
  sample.mPhredLikelihoods = {10, 0, 100};
  sample.mAlleleDepths = {30, 20};
  sample.mFwdAlleleDepths = {15, 10};
  sample.mRevAlleleDepths = {15, 10};
  sample.mRmsMappingQualities = {59.5, 60.0};
  sample.mNormPosteriorBQs = {30.1, 35.2};
  sample.mContinuousMixtureLods = {0.0, 12.5432};
  sample.mStrandBias = 2.45F;
  sample.mSoftClipAsym = 0.5F;
  sample.mFragLenDelta = 5.2F;
  sample.mReadPosCohenD = 0.1F;
  sample.mBaseQualCohenD = 0.2F;
  sample.mMapQualCohenD = -0.5F;
  sample.mAlleleMismatchDelta = 0.05F;
  sample.mSiteDepthFoldChange = 1.25F;
  sample.mPolarRadius = 5.51F;
  sample.mPolarAngle = 0.588F;
  sample.mTotalDepth = 50;
  sample.mGenotypeQuality = 10;
  sample.mGenotypeIndices = {0, 1};
  sample.mIsMissingSupport = false;

  auto const vcf_str = sample.RenderVcfString();
  REQUIRE(vcf_str == "0/1:30,20:15,10:15,10:50:59.5,60.0:30.1,35.2:2.450:0.5000:5.2:"
                     "0.1000:0.2000:-0.5000:0.050:1.25:5.5100:0.5880:12.5432:"
                     "10,0,100:10");
  sample.mIsMissingSupport = true;
  REQUIRE(sample.RenderVcfString() == "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.");
}

TEST_CASE("SampleGenotypeData renders nullopt metrics as VCF missing value dot",
          "[lancet][caller][VariantCall]") {
  VariantCall::SampleGenotypeData sample;
  sample.mPhredLikelihoods = {0, 10, 50};
  sample.mAlleleDepths = {10, 0};
  sample.mFwdAlleleDepths = {5, 0};
  sample.mRevAlleleDepths = {5, 0};
  sample.mRmsMappingQualities = {60.0, 0.0};
  sample.mNormPosteriorBQs = {30.0, 0.0};
  sample.mContinuousMixtureLods = {0.0, 0.0};
  sample.mStrandBias = 0.0F;
  sample.mSoftClipAsym = 0.0F;
  // Leave 5 optional metrics at their default (std::nullopt) — untestable
  sample.mSiteDepthFoldChange = 1.0F;
  sample.mPolarRadius = 1.0F;
  sample.mPolarAngle = 0.0F;
  sample.mTotalDepth = 10;
  sample.mGenotypeQuality = 0;
  sample.mGenotypeIndices = {0, 0};
  sample.mIsMissingSupport = false;

  auto const vcf_str = sample.RenderVcfString();
  // FLD, RPCD, BQCD, MQCD, ASMD should all render as "." (nullopt -> dot)
  REQUIRE(vcf_str == "0/0:10,0:5,0:5,0:10:60.0,0.0:30.0,0.0:0.000:0.0000:.:"
                     ".:.:.:.:1.00:1.0000:0.0000:0.0000:"
                     "0,10,50:0");
}

}  // namespace lancet::caller::tests
