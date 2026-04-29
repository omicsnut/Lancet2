#ifndef SRC_LANCET_CORE_VARIANT_BUILDER_H_
#define SRC_LANCET_CORE_VARIANT_BUILDER_H_

#include "lancet/base/tar_gz_writer.h"
#include "lancet/base/types.h"
#include "lancet/caller/genotyper.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_call.h"
#include "lancet/caller/variant_set.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/probe_results_writer.h"
#include "lancet/core/probe_diagnostics.h"
#include "lancet/core/read_collector.h"
#include "lancet/core/sample_info.h"
#include "lancet/core/variant_annotator.h"
#include "lancet/core/window.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace lancet::core {

class VariantBuilder {
 public:
  static constexpr u32 MIN_PHRED_SCORE = 0;
  static constexpr u32 MAX_PHRED_SCORE = 255;

  struct Params {
    // ── 8B Align ────────────────────────────────────────────────────────────
    /// Final user-facing archive path; ends in `.tar.gz`. Empty when
    /// `--out-graphs-tgz` is unset (no graph outputs are written).
    std::filesystem::path mOutGraphsTgz;

    /// Hidden subdirectory next to the final archive where each worker
    /// writes its `worker_<idx>.tar.gz` shard during the run. Populated
    /// by PipelineRunner only when `mOutGraphsTgz` is set; the merge step
    /// in PipelineRunner removes this directory after a successful merge.
    std::filesystem::path mShardsDir;

    std::filesystem::path mProbeVariantsPath;  // input missed_variants.txt for k-mer probing
    std::filesystem::path mProbeResultsPath;   // output probe_results.tsv (CLI parsing only)
    std::shared_ptr<cbdg::ProbeIndex const> mProbeIndex;  // precomputed global k-mer index
    std::shared_ptr<cbdg::ProbeResultsWriter> mProbeResultsWriter;  // thread-safe TSV writer

    /// Global genome GC fraction for LongdustQ bias correction.
    /// Default: 0.41 (human genome-wide average, Lander et al. 2001,
    /// Piovesan et al. 2019, Nurk et al. 2022 T2T-CHM13).
    /// Set to 0.5 for uniform distribution (no GC correction).
    /// See --genome-gc-bias CLI parameter.
    f64 mGcFraction = 0.41;

    cbdg::GraphParams mGraphParams;
    ReadCollector::Params mRdCollParams;

    /// Pre-built sample list — computed once at pipeline startup, immutable.
    /// Shared across all threads via shared_ptr<Params const>.
    std::vector<SampleInfo> mSampleList;

    // ── 1B Align ────────────────────────────────────────────────────────────
    bool mSkipActiveRegion = false;
  };

  /// `worker_index` is assigned by PipelineExecutor when constructing
  /// the worker pool (range `[0, num_threads)`). It is used to derive the
  /// per-worker shard filename `worker_<idx>.tar.gz` under
  /// `params->mShardsDir`. When `params->mShardsDir` is empty (i.e.
  /// `--out-graphs-tgz` unset), no shard is opened and snapshot emission
  /// is disabled entirely.
  VariantBuilder(std::shared_ptr<Params const> params, u32 window_length, u32 worker_index);

  enum class StatusCode : u8 {
    UNKNOWN = 0,

    SKIPPED_NONLY_REF_BASES = 1,
    SKIPPED_REF_REPEAT_SEEN = 2,
    SKIPPED_INACTIVE_REGION = 3,
    SKIPPED_ANCHOR_COVERAGE = 4,
    SKIPPED_NOASM_HAPLOTYPE = 5,
    MISSING_NO_MSA_VARIANTS = 6,
    FOUND_GENOTYPED_VARIANT = 7
  };

  [[nodiscard]] auto CurrentStatus() const noexcept -> StatusCode { return mCurrentCode; }

  using WindowResults = std::vector<std::unique_ptr<caller::VariantCall>>;
  [[nodiscard]] auto ProcessWindow(std::shared_ptr<Window const> const& window) -> WindowResults;

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  cbdg::Graph mDebruijnGraph;
  ReadCollector mReadCollector;
  caller::Genotyper mGenotyper;
  std::shared_ptr<Params const> mParamsPtr;
  caller::MsaBuilder mSpoaState;

  /// Variant annotator — produces coverage-invariant complexity features per variant.
  VariantAnnotator mAnnotator;

  /// Post-graph probe analysis: MSA extraction and genotyper read assignment.
  ProbeDiagnostics mProbeDiagnostics;

  /// Per-worker gzipped TAR shard. All graph outputs (DOTs from the
  /// snapshot buffer, GFA + FASTA from MSA extraction) for windows
  /// processed by this worker thread are appended as TAR entries here.
  /// Null when `--out-graphs-tgz` is unset (zero overhead). PipelineRunner
  /// merges the shards into the final user-visible archive after all
  /// workers have joined.
  std::unique_ptr<base::TarGzWriter> mGraphShardWriter;

  // ── 1B Align ────────────────────────────────────────────────────────────
  StatusCode mCurrentCode = StatusCode::UNKNOWN;

  // ── ProcessWindow helpers ───────────────────────────────────────────────
  [[nodiscard]] auto ShouldSkipWindow(Window const& window) -> bool;

  [[nodiscard]] auto ExtractVariants(cbdg::ComponentResult const& component, usize component_id,
                                     Window const& window) -> caller::VariantSet;

  static void CollectSupportedCalls(caller::VariantSet const& extracted,
                                    caller::Genotyper::Result& geno_result,
                                    absl::Span<SampleInfo const> samples, usize window_length,
                                    WindowResults& output_variant_calls);
};

[[nodiscard]] auto ToString(VariantBuilder::StatusCode status_code) -> std::string;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_VARIANT_BUILDER_H_
