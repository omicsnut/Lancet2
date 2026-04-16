#include "lancet/cli/pipeline_runner.h"

#ifdef LANCET_PROFILE_MODE
#include "gperftools/profiler.h"
#endif

#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"
#include "lancet/cli/cli_params.h"
#include "lancet/cli/vcf_header_builder.h"
#include "lancet/core/active_region_detector.h"
#include "lancet/core/input_spec_parser.h"
#include "lancet/core/pipeline_executor.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/window_builder.h"
#include "lancet/hts/bgzf_ostream.h"
#include "lancet/hts/uri_utils.h"

#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/format.h"

#include <filesystem>
#include <memory>
#include <utility>

#include <cstdlib>

namespace lancet::cli {

PipelineRunner::PipelineRunner(std::shared_ptr<CliParams> params) : mParamsPtr(std::move(params)) {
#ifdef LANCET_PROFILE_MODE
  setenv("CPUPROFILE_PER_THREAD_TIMERS", "1", 1);
  setenv("CPUPROFILE_FREQUENCY", "10000", 1);
  auto const timestamp = absl::FormatTime("%Y%m%d%ET%H%M%S", absl::Now(), absl::LocalTimeZone());
  auto const fname = fmt::format("Lancet.cpu_profile.{}.bin", timestamp);
  ProfilerStart(fname.c_str());
#endif
}

// ---------------------------------------------------------------------------
// Run: primary pipeline entry point
//
// Validates params → sets up output → constructs VCF header → delegates
// the full batch-fed execution loop to PipelineExecutor.
// ---------------------------------------------------------------------------
void PipelineRunner::Run() {
  base::Timer timer;
  ValidateAndPopulateParams();
  SetupGraphOutputDir();

  hts::BgzfOstream output_vcf;
  OpenOutputVcf(output_vcf);
  output_vcf << BuildVcfHeader(*mParamsPtr);
  output_vcf.flush();

  // Initialize the window builder with sorted regions
  core::WindowBuilder window_builder(mParamsPtr->mVariantBuilder.mRdCollParams.mRefPath,
                                     mParamsPtr->mWindowBuilder);

  window_builder.AddBatchRegions(absl::MakeConstSpan(mParamsPtr->mInRegions));
  window_builder.AddBatchRegions(mParamsPtr->mBedFile);

  if (window_builder.IsEmpty()) {
    LOG_WARN("No input regions provided to build windows."
             " Using contigs in reference as input regions")
    window_builder.AddAllReferenceRegions();
  }

  // Sort input regions before batch emission to ensure deterministic genomic ordering
  window_builder.SortInputRegions();
  core::PipelineExecutor executor(
      std::move(window_builder),
      std::make_shared<core::VariantBuilder::Params const>(mParamsPtr->mVariantBuilder),
      mParamsPtr->mNumWorkerThreads, mParamsPtr->mWindowBuilder.mWindowLength);

  auto const stats = executor.Execute(output_vcf);
  output_vcf.Close();
  core::PipelineExecutor::LogWindowStats(stats);

  auto const total_rt = absl::FormatDuration(absl::Trunc(timer.Runtime(), absl::Milliseconds(1)));
  LOG_INFO("Successfully completed processing | Runtime={}", total_rt)
  std::exit(EXIT_SUCCESS);
}

// ---------------------------------------------------------------------------
// SetupGraphOutputDir — propagate and recreate graph debug output directory
// ---------------------------------------------------------------------------

void PipelineRunner::SetupGraphOutputDir() {
  if (!mParamsPtr->mVariantBuilder.mOutGraphsDir.empty()) {
    mParamsPtr->mVariantBuilder.mGraphParams.mOutGraphsDir =
        mParamsPtr->mVariantBuilder.mOutGraphsDir;
    std::filesystem::remove_all(mParamsPtr->mVariantBuilder.mOutGraphsDir);
    std::filesystem::create_directories(mParamsPtr->mVariantBuilder.mOutGraphsDir);
  }
}

// ---------------------------------------------------------------------------
// OpenOutputVcf — resolve path, validate cloud credentials, open BGZF stream
//
// Local paths are made absolute and parent directories are created.
// Cloud URIs (gs://, s3://) trigger an immediate zero-byte HTTP PUT via
// hopen("w") to validate credentials upfront rather than discovering auth
// failures after a 40-hour pipeline run.
// ---------------------------------------------------------------------------
void PipelineRunner::OpenOutputVcf(hts::BgzfOstream& output_vcf) {
  // Resolve local paths to absolute and ensure parent directories exist.
  // Cloud URIs (gs://, s3://) bypass local path resolution entirely.
  if (!hts::IsCloudUri(mParamsPtr->mOutVcfGz.string())) {
    mParamsPtr->mOutVcfGz = std::filesystem::absolute(mParamsPtr->mOutVcfGz);
    if (!std::filesystem::exists(mParamsPtr->mOutVcfGz.parent_path())) {
      std::filesystem::create_directories(mParamsPtr->mOutVcfGz.parent_path());
    }
  }

  // AWS and GCP use 5MB+ multipart chunk caching over libcurl. Small VCF
  // headers won't reach this threshold, so libcurl defers the HTTP handshake
  // until BgzfOstream::Close() flushes after the full pipeline run.
  // To avoid a silent authentication failure after 40 hours, force an
  // immediate zero-byte HTTP PUT via hopen("w") to validate cloud credentials
  // upfront before starting the pipeline.
  if (hts::IsCloudUri(mParamsPtr->mOutVcfGz.string())) {
    auto const err = hts::ValidateCloudAccess(mParamsPtr->mOutVcfGz.string(), "w");
    if (!err.empty()) {
      LOG_CRITICAL("Cloud authentication failed! Cannot write to remote bucket: {}",
                   mParamsPtr->mOutVcfGz.string())
      std::exit(EXIT_FAILURE);
    }
  }

  if (!output_vcf.Open(mParamsPtr->mOutVcfGz, hts::BgzfFormat::VCF)) {
    LOG_CRITICAL("Could not open output VCF file: {}", mParamsPtr->mOutVcfGz.string())
    std::exit(EXIT_FAILURE);
  }
}

// ---------------------------------------------------------------------------
// ValidateAndPopulateParams — pre-flight checks before pipeline execution
//
// 1. Determines case-control mode from parsed sample specs. This must happen
//    BEFORE the skip-active-region early-return so VCF SHARED/CTRL/CASE INFO
//    headers are emitted correctly regardless of --no-active-region.
// 2. Validates that all input BAM/CRAM files contain MD tags. If any file
//    lacks MD tags, active region detection is disabled with a warning.
// ---------------------------------------------------------------------------
void PipelineRunner::ValidateAndPopulateParams() {
  auto const& rdcoll = mParamsPtr->mVariantBuilder.mRdCollParams;
  auto const all_specs =
      core::ParseAllInputSpecs(rdcoll.mCtrlPaths, rdcoll.mCasePaths, rdcoll.mSampleSpecs);

  // Case-control mode: true when both control and case samples exist.
  // Computed before skip-active-region so VCF headers are always correct.
  auto const has_label = [&all_specs](cbdg::Label::Tag label_tag) -> bool {
    return std::ranges::any_of(all_specs,
                               [label_tag](auto const& spec) { return spec.mTag == label_tag; });
  };

  mParamsPtr->mIsCaseCtrlMode = has_label(cbdg::Label::CASE) && has_label(cbdg::Label::CTRL);
  if (mParamsPtr->mVariantBuilder.mSkipActiveRegion) return;

  auto const missing_md = std::ranges::find_if(all_specs, [&rdcoll](auto const& spec) {
    return core::IsMissingMdTag(spec.mFilePath, rdcoll.mRefPath);
  });

  if (missing_md != all_specs.end()) {
    LOG_WARN("MD tag missing in BAM/CRAM: {}. Turning off active region detection",
             missing_md->mFilePath.filename().string())
    mParamsPtr->mVariantBuilder.mSkipActiveRegion = true;
  }
}

}  // namespace lancet::cli
