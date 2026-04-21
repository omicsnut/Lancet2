#ifndef SRC_LANCET_CLI_PIPELINE_RUNNER_H_
#define SRC_LANCET_CLI_PIPELINE_RUNNER_H_

#include "lancet/cli/cli_params.h"
#include "lancet/hts/bgzf_ostream.h"

#include <memory>

namespace lancet::cli {

class PipelineRunner {
 public:
  explicit PipelineRunner(std::shared_ptr<CliParams> params);

  [[noreturn]] void Run();

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::shared_ptr<CliParams> mParamsPtr;

  // ---------------------------------------------------------------------------
  // Modularized helpers extracted from the monolithic Run() method
  // ---------------------------------------------------------------------------

  /// Validates BAM/CRAM inputs and populates derived params (e.g. MD tag check).
  void ValidateAndPopulateParams();

  /// Propagates graph output directory from CLI params and recreates it on disk.
  void SetupGraphOutputDir();

  /// Resolves the output VCF path (local or cloud), validates credentials,
  /// and opens the BGZF output stream. Exits on failure.
  void OpenOutputVcf(hts::BgzfOstream& output_vcf);
};

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_PIPELINE_RUNNER_H_
