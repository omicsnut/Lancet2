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

#include "absl/base/call_once.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"
#include "spoa/alignment_engine.hpp"

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace lancet::core {

namespace {

/*
 * ============================================================================
 * SPOA MSA Parameter Rationale for Lancet2 Variant Extraction
 * ============================================================================
 * Values: Match: 0, Mismatch: -6, Gap1: -6,-2, Gap2: -26,-1
 *
 * **SIMD Lane Width Note**: All classical parameters (typically +2/-4) have
 * been shifted downwards by 2.  Setting the Match score to 0 keeps all
 * runtime scores non-positive, which keeps SPOA's WorstCaseAlignmentScore()
 * above the int16 threshold — so the engine selects the faster int16 SIMD
 * path (16 lanes per AVX2 register) instead of falling back to int32
 * (8 lanes, half throughput).  SPOA 4.1.5 dynamically dispatches between
 * int16 and int32 via WorstCaseAlignmentScore().
 * With Match=0 and our gap parameters, the worst-case score for typical
 * alignments (~-2300) is safely above int16 minimum (-32768).
 *
 * Unlike minimap2's `asm5` preset (which aggressively splits contigs at major
 * divergences for whole-genome synteny filtering), these parameters are tuned
 * to force end-to-end global alignment within a specific micro-assembly window
 * to capture dense somatic mutations and large Insertions/Deletions.
 *
 * 1. Why Convex (Dual-Affine) vs. Affine or Linear Scoring:
 *    - Linear Scoring applies a flat penalty per gap base, which is biologically
 *      inaccurate (one 50bp deletion is one biological event, not fifty 1bp
 *      independent events).
 *    - Single Affine Scoring forces a compromise: tune for small variants (strict
 *      extension) and you penalize/clip large insertions/deletions; tune for large
 *      insertions/deletions (loose extension) and sequencer noise creates messy,
 *      spurious small gaps.
 *    - Convex (Dual-Affine) Scoring solves this by taking the minimum of two
 *      intersecting models. It is strict for short gaps to suppress sequencer
 *      noise, but switches to a cheap extension penalty for large
 *      biological gaps.
 *
 * 2. Mismatch Tolerance (Multi-Nucleotide Variants / MNVs):
 *    asm5 uses a +1 match / -19 mismatch, which shatters alignments at dense
 *    mutation clusters. We use 0 / -6. The geometrical difference (6) keeps
 *    the MSA robustly intact while globally forcing alignments through complex variants.
 *
 * 3. Micro-Indel Sensitivity (Convex Model 1: -6, -2):
 *    asm5's -39 gap open penalty prevents small indels, forcing them to misalign
 *    as false-positive SNPs. Our -6 open / -2 extend penalty allows true small
 *    biological indels to open naturally while still applying enough friction
 *    to prevent 1bp sequencing errors (e.g., homopolymer stutters) from opening gaps.
 *
 * 4. Large Insertion/Deletion Continuity (Convex Model 2: -26, -1):
 *    asm5's -81 penalty for large gaps will soft-clip contigs right at an insertion/deletion
 *    breakpoint. Our parameters mathematically intersect at exactly 20bp (6 + 2L = 26 + 1L).
 *    For gaps > 20bp, the algorithm switches to Model 2 where the extension cost
 *    drops to -1. This "cheap extension" forces the DP matrix into mapping massive
 *    insertions/deletions as single, contiguous blocks in the MSA rather than
 *    dropping the alignment entirely.
 *
 *    https://curiouscoding.nl/posts/pairwise-alignment ->
 *    – Convex dual affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
 */
constexpr i8 MSA_MATCH_SCORE = 0;
constexpr i8 MSA_MISMATCH_SCORE = -6;
constexpr i8 MSA_OPEN1_SCORE = -6;
constexpr i8 MSA_EXTEND1_SCORE = -2;
constexpr i8 MSA_OPEN2_SCORE = -26;
constexpr i8 MSA_EXTEND2_SCORE = -1;
constexpr u8 DNA_ALPHABET_SIZE = 4;
constexpr u32 PREALLOC_WINDOW_LENGTH_MULTIPLIER = 3;

// ============================================================================
// AnnotatePathMetrics: populate haplotype count + path depth CV on a RawVariant.
//
// Sets two mutable fields on the variant (same pattern as mGraphMetrics):
//   mNumTotalHaps — total SPOA paths in this component (for HSE normalization)
//   mMaxPathCv    — max path depth CV across ALT paths (for PDCV FORMAT field)
// ============================================================================
void AnnotatePathMetrics(caller::RawVariant const& var, absl::Span<cbdg::Path const> comp_paths) {
  var.mNumTotalHaps = comp_paths.size();

  // Max path depth CV across ALT paths — O(nhaps) single pass.
  std::optional<f64> max_cv;
  for (usize hap = 1; hap < comp_paths.size(); ++hap) {
    auto const cv_val = comp_paths[hap].CoefficientOfVariationCoverage();
    max_cv = max_cv.has_value() ? std::max(*max_cv, cv_val) : cv_val;
  }
  var.mMaxPathCv = max_cv;
}

// ============================================================================
// CollectPathSequences: extract haplotype sequences from graph paths.
// ============================================================================
auto CollectPathSequences(absl::Span<cbdg::Path const> paths) -> std::vector<std::string> {
  std::vector<std::string> seqs;
  seqs.reserve(paths.size());
  std::ranges::transform(paths, std::back_inserter(seqs), [](auto const& path) -> std::string {
    return std::string(path.Sequence());
  });
  return seqs;
}

// ============================================================================
// CollectHapWeights: extract per-base SPOA weights from graph paths.
// ============================================================================
auto CollectHapWeights(absl::Span<cbdg::Path const> paths) -> cbdg::Path::HapWeights {
  cbdg::Path::HapWeights path_weights;
  path_weights.reserve(paths.size());
  std::ranges::transform(
      paths, std::back_inserter(path_weights),
      [](auto const& path) -> cbdg::Path::BaseWeights { return path.PerBaseWeights(); });
  return path_weights;
}

// ============================================================================
// BuildPerSampleCov: compute per-sample window coverage for SDFC normalization.
// ============================================================================
auto BuildPerSampleCov(absl::Span<SampleInfo const> samples, usize win_len)
    -> caller::VariantCall::PerSampleCov {
  caller::VariantCall::PerSampleCov result;
  for (auto const& sinfo : samples) {
    result[sinfo.SampleName()] = sinfo.SampledCov(win_len);
  }
  return result;
}

// ============================================================================
// CountAssembledHaplotypes: total non-REF haplotypes across all graph components.
// ============================================================================
auto CountAssembledHaplotypes(absl::Span<std::vector<cbdg::Path> const> components) -> u64 {
  return std::accumulate(components.begin(), components.end(), u64{0},
                         [](u64 sum, auto const& comp) -> u64 { return sum + comp.size() - 1; });
}

// ============================================================================
// HasAltSupport: true if any sample has > 0 ALT-supporting reads.
// ============================================================================
auto HasAltSupport(caller::SupportArray const& evidence) -> bool {
  return std::ranges::any_of(evidence, [](auto const& item) -> bool {
    return item.mData && item.mData->TotalAltCov() > 0;
  });
}

}  // namespace

VariantBuilder::VariantBuilder(std::shared_ptr<Params const> params, u32 const window_length)
    : mDebruijnGraph(params->mGraphParams),
      mReadCollector(params->mRdCollParams, absl::MakeConstSpan(params->mSampleList)),
      mParamsPtr(std::move(params)),
      mSpoaState{.mEngine = spoa::AlignmentEngine::Create(
                     spoa::AlignmentType::kNW, MSA_MATCH_SCORE, MSA_MISMATCH_SCORE, MSA_OPEN1_SCORE,
                     MSA_EXTEND1_SCORE, MSA_OPEN2_SCORE, MSA_EXTEND2_SCORE),
                 .mGraph = spoa::Graph()},
      mAnnotator(mParamsPtr->mGcFraction) {
  mSpoaState.mEngine->Prealloc(window_length * PREALLOC_WINDOW_LENGTH_MULTIPLIER,
                               DNA_ALPHABET_SIZE);

  // Initialize probe diagnostics and wire ProbeTracker into the graph.
  mProbeDiagnostics.Initialize(mParamsPtr->mProbeVariantsPath, mParamsPtr->mProbeResultsWriter,
                               mParamsPtr->mProbeIndex);
  mDebruijnGraph.SetProbeTracker(mProbeDiagnostics.Tracker());
}

auto VariantBuilder::ProcessWindow(std::shared_ptr<Window const> const& window) -> WindowResults {
  auto const region = window->AsRegionPtr();
  auto const reg_str = region->ToSamtoolsRegion();
  static thread_local auto const THREAD_ID =
      absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Processing window {} in thread {:#x}", reg_str, THREAD_ID)

  if (std::ranges::all_of(window->SeqView(), [](char base) { return base == 'N'; })) {
    LOG_DEBUG("Skipping window {} since it has only N bases in reference", reg_str)
    mCurrentCode = StatusCode::SKIPPED_NONLY_REF_BASES;
    return {};
  }

  if (lancet::base::HasExactRepeat(
          lancet::base::SlidingView(window->SeqView(), mParamsPtr->mGraphParams.mMaxKmerLen))) {
    LOG_DEBUG("Skipping window {} since reference has repeat {}-mers", reg_str,
              mParamsPtr->mGraphParams.mMaxKmerLen)
    mCurrentCode = StatusCode::SKIPPED_REF_REPEAT_SEEN;
    return {};
  }

  if (!mParamsPtr->mSkipActiveRegion &&
      !core::IsActiveRegion(mReadCollector.SampleList(), mReadCollector.Extractors(), *region)) {
    LOG_DEBUG("Skipping window {} since it has no evidence of mutation in any sample", reg_str)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Collecting all available sample reads for window {}", reg_str)
  auto const rc_result = mReadCollector.CollectRegionResult(*region);
  absl::Span<cbdg::Read const> const reads = absl::MakeConstSpan(rc_result.mSampleReads);
  absl::Span<SampleInfo const> const samples = absl::MakeConstSpan(rc_result.mSampleList);

  auto const total_cov = SampleInfo::CombinedSampledCov(samples, window->Length());
  if (total_cov < static_cast<f64>(mParamsPtr->mGraphParams.mMinAnchorCov)) {
    LOG_DEBUG("Skipping window {} since it has only {:.2f}x total sample coverage", reg_str,
              total_cov)
    mCurrentCode = StatusCode::SKIPPED_INACTIVE_REGION;
    return {};
  }

  LOG_DEBUG("Building graph for {} with {} sample reads and {:.2f}x total coverage", reg_str,
            reads.size(), total_cov)
  // First haplotype from each component will always be the reference haplotype sequence for the
  // graph
  auto const dbg_rslt = mDebruijnGraph.BuildComponentHaplotypes(window->AsRegionPtr(), reads);
  auto const& component_haplotypes = dbg_rslt.mGraphHaplotypes;

  auto const num_asm_haps = CountAssembledHaplotypes(absl::MakeConstSpan(component_haplotypes));
  if (num_asm_haps == 0) {
    LOG_DEBUG("Could not assemble any haplotypes for window {} with k={}", reg_str,
              mDebruijnGraph.CurrentK())
    mCurrentCode = StatusCode::SKIPPED_NOASM_HAPLOTYPE;
    return {};
  }

  WindowResults variants;
  auto const per_sample_cov = BuildPerSampleCov(samples, window->Length());
  for (usize idx = 0; idx < component_haplotypes.size(); ++idx) {
    auto const nhaps = component_haplotypes[idx].size();
    auto const anchor_start = window->StartPos1() + dbg_rslt.mAnchorStartIdxs[idx];
    std::vector<cbdg::Path> const& comp_paths = component_haplotypes[idx];

    auto const comp_haps = CollectPathSequences(comp_paths);
    auto const comp_weights = CollectHapWeights(comp_paths);

    LOG_DEBUG("Building MSA for graph component {} from window {} with {} haplotypes", idx, reg_str,
              nhaps)

    mSpoaState.UpdateSpoaState(absl::MakeConstSpan(comp_haps), absl::MakeConstSpan(comp_weights));
    // SerializeGraph is a no-op when MakeGfaPath returns an empty path
    // (i.e., --out-graphs-dir was not specified on the CLI).
    mSpoaState.SerializeGraph(MakeGfaPath(*window, idx));

    caller::VariantSet const vset(mSpoaState.mGraph, *window, anchor_start);

    // Annotate complexity features on every variant
    mAnnotator.AnnotateSequenceComplexity(vset, absl::MakeConstSpan(comp_haps));
    LANCET_ASSERT(idx < dbg_rslt.mComponentMetrics.size())
    VariantAnnotator::AnnotateGraphComplexity(vset, dbg_rslt.mComponentMetrics[idx]);

    if (vset.IsEmpty()) {
      LOG_DEBUG("No variants found in graph component {} for window {} with {} haplotypes", idx,
                reg_str, nhaps)
      continue;
    }

    // Check whether probe variants appear in the MSA-extracted variant set.
    mProbeDiagnostics.CheckMsaExtraction(vset, *window, idx);

    LOG_DEBUG("Found variant(s) in graph component {} for window {} with {} haplotypes", idx,
              reg_str, nhaps)
    auto genotyped = mGenotyper.Genotype(absl::MakeConstSpan(comp_haps), reads, vset);

    // Check genotyper read assignment for probe variants.
    mProbeDiagnostics.CheckGenotyperResult(genotyped, vset, idx);
    for (auto const& var : vset) {
      auto iter = genotyped.find(&var);
      caller::VariantCall::SupportsByVariant var_supports;
      if (iter != genotyped.end() && HasAltSupport(iter->second)) {
        var_supports.emplace(&var, std::move(iter->second));
      }

      if (var_supports.empty()) continue;

      AnnotatePathMetrics(var, comp_paths);
      variants.emplace_back(std::make_unique<caller::VariantCall>(&var, std::move(var_supports),
                                                                  samples, per_sample_cov));
    }
  }

  // Flush this window's probe records to the shared writer (no-op if inactive).
  mProbeDiagnostics.SubmitCompleted();

  if (variants.empty()) {
    LOG_DEBUG("No variants found for window {} from {} assembled graph paths", reg_str,
              num_asm_haps)
    mCurrentCode = StatusCode::MISSING_NO_MSA_VARIANTS;
    return {};
  }

  mCurrentCode = StatusCode::FOUND_GENOTYPED_VARIANT;
  LOG_DEBUG("Genotyped {} variant(s) for window {} by re-aligning sample reads", variants.size(),
            reg_str)
  return variants;
}

auto VariantBuilder::MakeGfaPath(Window const& win, usize const comp_id) const
    -> std::filesystem::path {
  if (mParamsPtr->mOutGraphsDir.empty()) return {};

  auto const graph_dir = mParamsPtr->mOutGraphsDir / "poa_graph";
  static absl::once_flag dir_created;
  absl::call_once(dir_created, [&] { std::filesystem::create_directories(graph_dir); });

  return graph_dir / fmt::format("msa__{}_{}_{}__c{}.gfa", win.ChromName(), win.StartPos1(),
                                 win.EndPos1(), comp_id);
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
