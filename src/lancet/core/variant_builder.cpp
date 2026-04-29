#include "lancet/core/variant_builder.h"

#include "lancet/base/assert.h"
#include "lancet/base/logging.h"
#include "lancet/base/repeat.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/caller/msa_builder.h"
#include "lancet/caller/variant_call.h"
#include "lancet/caller/variant_set.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/active_region_detector.h"
#include "lancet/core/sample_info.h"
#include "lancet/core/window.h"

#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"
#include "spoa/alignment_engine.hpp"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace lancet::core {

namespace {

// ============================================================================
// HasAltSupport: true if any sample has > 0 ALT-supporting reads.
// ============================================================================
auto HasAltSupport(caller::SupportArray const& evidence) -> bool {
  return std::ranges::any_of(evidence, [](auto const& item) -> bool {
    return item.mData && item.mData->TotalAltCov() > 0;
  });
}

// ============================================================================
// SerializeSpoaState: Appends this component's GFA + FASTA to the per-worker
// tar.gz shard when `--out-graphs-tgz` is set (shard_writer is non-null).
// Both files for the component become regular-file TAR entries under the
// `poa_graph/<window>/` archive subdir. When shard_writer is null,
// the entire MSA-rendering path is skipped — zero overhead in production.
// ============================================================================
void SerializeSpoaState(caller::MsaBuilder& state, Window const& window, u32 component_id,
                        base::TarGzWriter& shard_writer) {
  auto const start1 = window.StartPos1();
  auto const end1 = window.EndPos1();
  auto const region_prefix = fmt::format("{}_{}_{}", window.ChromName(), start1, end1);
  auto const window_subdir = std::filesystem::path{"poa_graph"} / region_prefix;

  auto const gfa_contents = state.BuildGfaString();
  auto const msa = state.mGraph.GenerateMultipleSequenceAlignment(false);
  auto const fasta_contents = caller::MsaBuilder::BuildFastaString(absl::MakeConstSpan(msa));

  auto const gfa_filename = fmt::format("msa__{}__c{}.gfa", region_prefix, component_id);
  auto const fasta_filename = fmt::format("msa__{}__c{}.fasta", region_prefix, component_id);
  auto const gfa_entry_path = window_subdir / gfa_filename;
  auto const fasta_entry_path = window_subdir / fasta_filename;
  shard_writer.AddRegularFileEntry(gfa_entry_path.string(), gfa_contents);
  shard_writer.AddRegularFileEntry(fasta_entry_path.string(), fasta_contents);
}

}  // namespace

VariantBuilder::VariantBuilder(ParamsPtr params, u32 window_len, u32 worker_id)
    : mDebruijnGraph(params->mGraphParams),
      mReadCollector(params->mRdCollParams, absl::MakeConstSpan(params->mSampleList)),
      mParamsPtr(std::move(params)),
      mSpoaState(lancet::caller::MsaBuilder()),
      mAnnotator(mParamsPtr->mGcFraction) {
  static constexpr u8 DNA_ALPHABET_SIZE = 4;
  static constexpr u32 PREALLOC_MULTIPLIER = 3;
  mSpoaState.mEngine->Prealloc(window_len * PREALLOC_MULTIPLIER, DNA_ALPHABET_SIZE);

  // Initialize probe diagnostics and wire ProbeTracker into the graph.
  mProbeDiagnostics.Initialize(mParamsPtr->mProbeVariantsPath, mParamsPtr->mProbeResultsWriter,
                               mParamsPtr->mProbeIndex);
  mDebruijnGraph.SetProbeTracker(mProbeDiagnostics.Tracker());

  // Open this worker's per-thread gzipped TAR shard if `--out-graphs-tgz`
  // is set (PipelineRunner populates `mShardsDir` only in that case). The
  // shard is opened with EndOfArchive::OMIT — the shard merger appends a
  // single end-of-archive marker once, after all shards are concatenated.
  // When `mShardsDir` is empty, the writer stays null and Graph's
  // snapshot-emission paths short-circuit on the null shard writer.
  if (!mParamsPtr->mShardsDir.empty()) {
    auto const shard_filename = fmt::format("worker_{}.tar.gz", worker_id);
    auto const shard_path = mParamsPtr->mShardsDir / shard_filename;
    mGraphShardWriter =
        std::make_unique<base::TarGzWriter>(shard_path, base::TarGzWriter::EndOfArchive::OMIT);
    mDebruijnGraph.SetGraphShardWriter(mGraphShardWriter.get());
  }
}

// ============================================================================
// ShouldSkipWindow: pre-read guard checks for N-only, repeat, and inactive.
// ============================================================================
auto VariantBuilder::ShouldSkipWindow(Window const& window) -> bool {
  auto const region_string = window.AsRegionPtr()->ToSamtoolsRegion();

  if (std::ranges::all_of(window.SeqView(), [](char base) { return base == 'N'; })) {
    LOG_DEBUG("Skipping window {} as reference contains only N bases", region_string)
    mCurrentCode = StatusCode::SKIPPED_NONLY_REF_BASES;
    return true;
  }

  auto const max_k = mParamsPtr->mGraphParams.mMaxKmerLen;
  if (lancet::base::HasExactRepeat(lancet::base::SlidingView(window.SeqView(), max_k))) {
    LOG_DEBUG("Skipping window {} as reference contains {}-mer repeats", region_string, max_k)
    mCurrentCode = StatusCode::SKIPPED_REF_REPEAT_SEEN;
    return true;
  }

  if (!mParamsPtr->mSkipActiveRegion &&
      !core::IsActiveRegion(mReadCollector.SampleList(), mReadCollector.Extractors(),
                            *window.AsRegionPtr())) {
    LOG_DEBUG("Skipping window {} as it has no evidence of mutation in any sample", region_string)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return true;
  }

  return false;
}

// ============================================================================
// ExtractVariants: build MSA from component haplotypes, extract variants
// from the alignment graph, and annotate complexity + path metrics.
//
// Uses ComponentResult's zero-copy HaplotypeSequenceViews() for SPOA/scorer
// and pre-computed MaxAltPathCv()/NumPaths() for path metric annotation.
// ============================================================================
auto VariantBuilder::ExtractVariants(cbdg::ComponentResult const& component,
                                     usize const component_id, Window const& window)
    -> caller::VariantSet {
  auto const region_string = window.AsRegionPtr()->ToSamtoolsRegion();
  auto const ref_anchor_pos1 = window.StartPos1() + component.AnchorStartOffset();
  auto const hap_views = component.HaplotypeSequenceViews();
  auto const weights = component.HaplotypeWeights();

  LOG_DEBUG("Building MSA for graph component {} from window {} with {} assembled haplotypes",
            component_id, region_string, component.NumPaths())

  mSpoaState.UpdateSpoaState(absl::MakeConstSpan(hap_views), absl::MakeConstSpan(weights));
  if (mGraphShardWriter) {
    SerializeSpoaState(mSpoaState, window, component_id, *mGraphShardWriter);
  }

  caller::VariantSet vset(mSpoaState.mGraph, window, ref_anchor_pos1);

  mAnnotator.AnnotateSequenceComplexity(vset, absl::MakeConstSpan(hap_views));
  VariantAnnotator::AnnotateGraphComplexity(vset, component.Metrics());

  // Path metrics annotation: nhaps + max ALT path depth CV.
  // Deterministic — depends only on component paths, not on genotyping.
  auto const max_alt_cv = component.MaxAltPathCv();
  for (auto const& var : vset) {
    var.mNumTotalHaps = component.NumPaths();
    var.mMaxPathCv = max_alt_cv;
  }

  if (vset.IsEmpty()) {
    LOG_DEBUG("No variants found in graph component {} for window {} with {} haplotypes",
              component_id, region_string, component.NumPaths())
  }

  return vset;
}

// ============================================================================
// CollectSupportedCalls: filter genotyped variants for ALT support and
// construct VariantCalls. Appends to output_variant_calls container.
// ============================================================================
void VariantBuilder::CollectSupportedCalls(caller::VariantSet const& extracted,
                                           caller::Genotyper::Result& geno_result,
                                           absl::Span<SampleInfo const> samples,
                                           usize const window_length,
                                           WindowResults& output_variant_calls) {
  for (auto const& var : extracted) {
    auto iter = geno_result.find(&var);
    caller::VariantCall::SupportsByVariant var_supports;
    if (iter != geno_result.end() && HasAltSupport(iter->second)) {
      var_supports.emplace(&var, std::move(iter->second));
    }

    if (var_supports.empty()) continue;

    output_variant_calls.emplace_back(std::make_unique<caller::VariantCall>(
        &var, std::move(var_supports), samples, window_length));
  }
}

auto VariantBuilder::ProcessWindow(std::shared_ptr<Window const> const& window) -> WindowResults {
  auto const region = window->AsRegionPtr();
  auto const region_string = region->ToSamtoolsRegion();

  static thread_local auto const CURRENT_TID = std::this_thread::get_id();
  static thread_local auto const THREAD_ID = absl::Hash<std::thread::id>()(CURRENT_TID);
  LOG_DEBUG("Processing window {} in thread {:#x}", region_string, THREAD_ID)

  // Phase 1: Pre-read qualification — N-only, repeat k-mers, active region
  if (ShouldSkipWindow(*window)) return {};

  // Phase 2: Read collection and depth qualification
  LOG_DEBUG("Collecting all available sample reads for window {}", region_string)
  auto const rc_result = mReadCollector.CollectRegionResult(*region);
  auto const reads = absl::MakeConstSpan(rc_result.mSampleReads);
  auto const samples = absl::MakeConstSpan(rc_result.mSampleList);

  auto const cross_sample_cov = SampleInfo::CrossSampleMeanCoverage(samples, window->Length());
  if (cross_sample_cov < static_cast<f64>(mParamsPtr->mGraphParams.mMinAnchorCov)) {
    LOG_DEBUG("Skipping window {} with {:.2f}x total coverage as min. anchor coverage is {}x",
              region_string, cross_sample_cov, mParamsPtr->mGraphParams.mMinAnchorCov)
    mCurrentCode = StatusCode::SKIPPED_ANCHOR_COVERAGE;
    return {};
  }

  // Phase 3: de Bruijn graph assembly and haplotype enumeration
  LOG_DEBUG("Building graph for {} with {} extracted sample reads and {:.2f}x total coverage",
            region_string, reads.size(), cross_sample_cov)
  auto const components = mDebruijnGraph.BuildComponentResults(window->AsRegionPtr(), reads);

  auto const num_assembled_haps = std::accumulate(
      components.cbegin(), components.cend(), u64{0},
      [](u64 const sum, auto const& comp) -> u64 { return sum + comp.NumAltHaplotypes(); });

  if (num_assembled_haps == 0) {
    LOG_DEBUG("Could not successfully assemble any haplotypes for window {} with k={}",
              region_string, mDebruijnGraph.CurrentK())
    mCurrentCode = StatusCode::SKIPPED_NOASM_HAPLOTYPE;
    return {};
  }

  // Phase 4: Per-component MSA, genotyping, and variant collection
  WindowResults variant_calls;
  for (usize component_idx = 0; component_idx < components.size(); ++component_idx) {
    auto const& component = components[component_idx];

    auto extracted = ExtractVariants(component, component_idx, *window);
    if (extracted.IsEmpty()) continue;

    // Probes with paths in skipped (empty) components keep all MSA flags false.
    // The Python attribution engine classifies these as msa_not_extracted.
    mProbeDiagnostics.CheckMsaExtraction(extracted, *window);
    LOG_DEBUG("Found variant(s) in graph component {} for window {} with {} haplotypes",
              component_idx, region_string, component.NumPaths())

    // HaplotypeSequences() allocates — only called when variants exist.
    // Genotyper's minimap2 requires null-terminated c_str() pointers.
    auto const hap_seqs = component.HaplotypeSequences();
    auto geno_result = mGenotyper.Genotype(hap_seqs, reads, extracted);
    mProbeDiagnostics.CheckGenotyperResult(geno_result, extracted);
    CollectSupportedCalls(extracted, geno_result, samples, window->Length(), variant_calls);
  }

  mProbeDiagnostics.SubmitCompleted();
  if (variant_calls.empty()) {
    LOG_DEBUG("No variants found for window {} from {} assembled graph contigs/haplotypes",
              region_string, num_assembled_haps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped {} variant(s) for window {} by re-aligning sample reads",
            variant_calls.size(), region_string)
  return variant_calls;
}

auto ToString(VariantBuilder::StatusCode const status_code) -> std::string {
  using enum VariantBuilder::StatusCode;

  switch (status_code) {
    case SKIPPED_NONLY_REF_BASES:
      return "SKIPPED_NONLY_REF_BASES";
    case SKIPPED_REF_REPEAT_SEEN:
      return "SKIPPED_REF_REPEAT_SEEN";
    case SKIPPED_INACTIVE_REGION:
      return "SKIPPED_INACTIVE_REGION";
    case SKIPPED_ANCHOR_COVERAGE:
      return "SKIPPED_ANCHOR_COVERAGE";
    case SKIPPED_NOASM_HAPLOTYPE:
      return "SKIPPED_NOASM_HAPLOTYPE";
    case MISSING_NO_MSA_VARIANTS:
      return "MISSING_NO_MSA_VARIANTS";
    case FOUND_GENOTYPED_VARIANT:
      return "FOUND_GENOTYPED_VARIANT";
    default:
      break;
  }

  return "UNKNOWN";
}

}  // namespace lancet::core
