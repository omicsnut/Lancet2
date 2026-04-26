#include "lancet/cbdg/probe_results_writer.h"

#include "lancet/base/logging.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/probe_tracker.h"

#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/hash/hash.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ranges.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

ProbeResultsWriter::ProbeResultsWriter(std::filesystem::path results_path,
                                       std::vector<ProbeVariant> variants)
    : mResultsPath(std::move(results_path)), mVariants(std::move(variants)) {}

// ============================================================================
// Append: thread-safe batch write of probe records to the TSV file.
//
// Acquires the mutex, writes the header on the first call, then appends
// formatted TSV lines. The file is opened in append mode per call rather
// than held open, so partial results survive crashes.
// ============================================================================
void ProbeResultsWriter::Append(absl::Span<ProbeKRecord const> records) {
  if (records.empty()) return;

  absl::MutexLock const lock(mMutex);
  std::ofstream out(mResultsPath, std::ios::app);
  if (!out.is_open()) {
    LOG_WARN("KMER_PROBE could not open results file for writing: {}", mResultsPath.string())
    return;
  }

  // clang-format off
  if (!mHeaderWritten) {
    static constexpr auto HEADER_LINE =
        "probe_id\tchrom\tpos\tref\talt\tn_tier1_reads\tkmer_size\tn_expected_alt_kmers\tn_alt_kmers_in_reads\t"
        "comp_id\tn_comp_nodes\tn_split_across_comps\tn_surviving_build\tn_surviving_lowcov1\tn_surviving_compress1\t"
        "n_surviving_lowcov2\tn_surviving_compress2\tn_surviving_tips\thap_indices\tis_traversal_limited\t"
        "is_cycle_retry\tis_complex_retry\tis_no_anchor\tis_variant_in_anchor\t"
        "n_last_edge_kmers\tn_last_interior_kmers\tn_last_chain_gaps\t"
        "msa_shift_bp\tis_msa_exact_match\tis_msa_shifted\tis_msa_representation\t"
        "n_geno_true_alt_reads\tn_geno_total_ref_reads\tn_geno_stolen_to_ref\t"
        "n_geno_stolen_to_wrong_alt\tn_geno_non_overlapping\tis_geno_has_alt_support\tis_geno_no_result\t"
        "lost_at_stage\n";
    out << HEADER_LINE;
    mHeaderWritten = true;
  }

  for (auto const& record : records) {
    auto const& probe = mVariants[record.mProbeId];
    auto const lost_at = DeriveLostAt(record);
    auto const paths_str = fmt::format("{}", fmt::join(record.mHapIndices, ","));
    mWrittenProbeIds.insert(record.mProbeId);

    out << fmt::format(
        "{PROBE_ID}\t{CHROM}\t{POS}\t{REF}\t{ALT}\t{N_TIER1_READS}\t{KMER_SIZE}\t"
        "{N_EXPECTED_ALT_KMERS}\t{N_ALT_KMERS_IN_READS}\t{COMP_ID}\t{N_COMP_NODES}\t"
        "{N_SPLIT_ACROSS_COMPS}\t{N_SURVIVING_BUILD}\t{N_SURVIVING_LOWCOV1}\t"
        "{N_SURVIVING_COMPRESS1}\t{N_SURVIVING_LOWCOV2}\t{N_SURVIVING_COMPRESS2}\t"
        "{N_SURVIVING_TIPS}\t{HAP_INDICES}\t{IS_TRAVERSAL_LIMITED}\t{IS_CYCLE_RETRY}\t"
        "{IS_COMPLEX_RETRY}\t{IS_NO_ANCHOR}\t{IS_VARIANT_IN_ANCHOR}\t{N_LAST_EDGE}\t"
        "{N_LAST_INTERIOR}\t{N_LAST_GAPS}\t{MSA_SHIFT}\t{IS_MSA_EXACT}\t{IS_MSA_SHIFTED}\t"
        "{IS_MSA_REPR}\t{N_GENO_ALT}\t{N_GENO_REF}\t{N_GENO_STOLEN_REF}\t{N_GENO_STOLEN_ALT}\t"
        "{N_GENO_NON_OVL}\t{IS_GENO_ALT}\t{IS_GENO_NORES}\t{LOST_AT_STAGE}\n",
        fmt::arg("PROBE_ID", record.mProbeId),
        fmt::arg("CHROM", probe.mChrom),
        fmt::arg("POS", probe.mGenomeStart0),
        fmt::arg("REF", probe.mRef),
        fmt::arg("ALT", probe.mAlt),
        fmt::arg("N_TIER1_READS", probe.mTier1AltCount),
        fmt::arg("KMER_SIZE", record.mKmerSize),
        fmt::arg("N_EXPECTED_ALT_KMERS", record.mExpectedAltKmers),
        fmt::arg("N_ALT_KMERS_IN_READS", record.mAltKmersInReads),
        fmt::arg("COMP_ID", record.mCompId),
        fmt::arg("N_COMP_NODES", record.mCompNumNodes),
        fmt::arg("N_SPLIT_ACROSS_COMPS", record.mSplitAcrossComps),
        fmt::arg("N_SURVIVING_BUILD", record.mSurvivingBuild),
        fmt::arg("N_SURVIVING_LOWCOV1", record.mSurvivingLowcov1),
        fmt::arg("N_SURVIVING_COMPRESS1", record.mSurvivingCompress1),
        fmt::arg("N_SURVIVING_LOWCOV2", record.mSurvivingLowcov2),
        fmt::arg("N_SURVIVING_COMPRESS2", record.mSurvivingCompress2),
        fmt::arg("N_SURVIVING_TIPS", record.mSurvivingTips),
        fmt::arg("HAP_INDICES", paths_str.empty() ? "." : paths_str),
        fmt::arg("IS_TRAVERSAL_LIMITED", static_cast<int>(record.mIsTraversalLimited)),
        fmt::arg("IS_CYCLE_RETRY", static_cast<int>(record.mIsCycleRetry)),
        fmt::arg("IS_COMPLEX_RETRY", static_cast<int>(record.mIsComplexRetry)),
        fmt::arg("IS_NO_ANCHOR", static_cast<int>(record.mIsNoAnchor)),
        fmt::arg("IS_VARIANT_IN_ANCHOR", static_cast<int>(record.mIsVariantInAnchor)),
        fmt::arg("N_LAST_EDGE", record.mLastEdgeCount),
        fmt::arg("N_LAST_INTERIOR", record.mLastInteriorCount),
        fmt::arg("N_LAST_GAPS", record.mLastGapCount),
        fmt::arg("MSA_SHIFT", record.mMsaShiftBp),
        fmt::arg("IS_MSA_EXACT", static_cast<int>(record.mIsMsaExactMatch)),
        fmt::arg("IS_MSA_SHIFTED", static_cast<int>(record.mIsMsaShifted)),
        fmt::arg("IS_MSA_REPR", static_cast<int>(record.mIsMsaRepresentation)),
        fmt::arg("N_GENO_ALT", record.mGenoTrueAltReads),
        fmt::arg("N_GENO_REF", record.mGenoTotalRefReads),
        fmt::arg("N_GENO_STOLEN_REF", record.mGenoStolenToRef),
        fmt::arg("N_GENO_STOLEN_ALT", record.mGenoStolenToWrongAlt),
        fmt::arg("N_GENO_NON_OVL", record.mGenoNonOverlapping),
        fmt::arg("IS_GENO_ALT", static_cast<int>(record.mIsGenoHasAltSupport)),
        fmt::arg("IS_GENO_NORES", static_cast<int>(record.mIsGenoNoResult)),
        fmt::arg("LOST_AT_STAGE", lost_at));
    // clang-format on
  }

  LOG_DEBUG("KMER_PROBE wrote {} probe records to {}", records.size(), mResultsPath.string())
}

// ============================================================================
// EmitUnprocessedProbes: emit default records for loaded variants never written.
//
// Probes that fall in skipped windows (N-only, repeat, inactive region) are
// never activated by GenerateAndTag, so no ProbeKRecord is ever created.
// This method ensures every input variant has exactly one output row.
// ============================================================================
void ProbeResultsWriter::EmitUnprocessedProbes() {
  std::vector<ProbeKRecord> unprocessed;
  {
    absl::MutexLock const lock(mMutex);
    for (u16 idx = 0; idx < static_cast<u16>(mVariants.size()); ++idx) {
      if (mWrittenProbeIds.contains(idx)) continue;
      unprocessed.push_back(ProbeKRecord{.mProbeId = idx, .mIsNotProcessed = true});
    }
  }

  if (!unprocessed.empty()) {
    LOG_INFO("KMER_PROBE emitting {} not-processed probe records", unprocessed.size())
    Append(absl::MakeConstSpan(unprocessed));
  }
}

// ============================================================================
// DeriveLostAt: priority-based attribution of where the variant was lost.
//
// All values use consistent snake_case. Categories (21-level cascade):
//   Not processed (highest priority):
//     0. not_processed     — probe never activated in any window (skipped region)
//   Structural (component-level impossibilities):
//     1. variant_in_anchor — variant overlaps source/sink k-mer range
//     2. no_anchor         — component had no valid source/sink anchors
//     3. short_anchor      — anchor region too short (<150bp)
//     4. cycle_retry       — cycle detected, k was increased
//     5. complex_retry     — graph too complex, k was increased
//   Pruning (first stage where all ALT k-mers drop to zero):
//     6-11. pruned_at_{build,lowcov1,compress1,lowcov2,compress2,tips}
//   Path enumeration:
//     12. traversal_limited — k-mers survived but BFS budget exhausted
//     13. no_path           — k-mers survived all pruning but no walk carries ALT
//   MSA extraction:
//     14. msa_no_variant    — variant survived path but not extracted by MSA
//     15. msa_shifted       — extracted at wrong position (shifted)
//     16. msa_representation — subsumed by larger MNV representation
//   Genotyper:
//     17. geno_no_result    — variant extracted but genotyper produced no result
//     18. geno_no_support   — genotyper ran but no reads assigned to truth ALT
//     19. geno_stolen       — ALT-carrying reads misassigned to REF or wrong ALT
//   Success:
//     20. survived          — variant found and genotyped with ALT support
// ============================================================================
auto ProbeResultsWriter::DeriveLostAt(ProbeKRecord const& record) -> std::string_view {
  if (record.mIsNotProcessed) return "not_processed";
  if (record.mIsVariantInAnchor) return "variant_in_anchor";
  if (record.mIsNoAnchor) return "no_anchor";
  if (record.mIsShortAnchor) return "short_anchor";
  if (record.mIsCycleRetry) return "cycle_retry";
  if (record.mIsComplexRetry) return "complex_retry";

  // Check each pruning stage in pipeline order for first drop to zero.
  // Strings match PRUNE_STAGE_NAMES so lost_at_stage values are consistent.
  if (record.mExpectedAltKmers > 0 && record.mSurvivingBuild == 0) {
    return "pruned_at_build";
  }

  if (record.mSurvivingBuild > 0 && record.mSurvivingLowcov1 == 0) {
    return "pruned_at_lowcov1";
  }

  if (record.mSurvivingLowcov1 > 0 && record.mSurvivingCompress1 == 0) {
    return "pruned_at_compress1";
  }

  if (record.mSurvivingCompress1 > 0 && record.mSurvivingLowcov2 == 0) {
    return "pruned_at_lowcov2";
  }

  if (record.mSurvivingLowcov2 > 0 && record.mSurvivingCompress2 == 0) {
    return "pruned_at_compress2";
  }

  if (record.mSurvivingCompress2 > 0 && record.mSurvivingTips == 0) {
    return "pruned_at_tips";
  }

  if (record.mIsTraversalLimited && record.mHapIndices.empty()) return "traversal_limited";
  if (record.mHapIndices.empty()) return "no_path";

  // MSA extraction attribution
  if (!record.mIsMsaExactMatch && !record.mIsMsaShifted && !record.mIsMsaRepresentation) {
    return "msa_no_variant";
  }

  if (record.mIsMsaShifted && !record.mIsMsaExactMatch) return "msa_shifted";
  if (record.mIsMsaRepresentation && !record.mIsMsaExactMatch) return "msa_representation";

  // Genotyper attribution
  if (record.mIsGenoNoResult) return "geno_no_result";
  if (!record.mIsGenoHasAltSupport) return "geno_no_support";
  if (record.mGenoStolenToRef > 0 || record.mGenoStolenToWrongAlt > 0) return "geno_stolen";

  return "survived";
}

}  // namespace lancet::cbdg
