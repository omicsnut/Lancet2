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
#include <utility>
#include <vector>

namespace lancet::cbdg {

ProbeResultsWriter::ProbeResultsWriter(std::filesystem::path results_path,
                                       std::vector<ProbeVariant> variants)
    : mResultsPath(std::move(results_path)), mVariants(std::move(variants)) {}

// ============================================================================
// Append: thread-safe batch write of probe records to the TSV file.
//
// Pure data emitter — writes every observed fact per (probe_id, window,
// comp_id, k) attempt. No attribution logic. The lost_at_stage derivation
// is performed downstream in the Python analysis script.
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
        "probe_id\tchrom\tpos\tref\talt\twindow\tn_raw_alt_reads\tkmer_size\tn_expected_alt_kmers\tn_alt_kmers_in_reads\t"
        "comp_id\tn_comp_nodes\tn_split_across_comps\tn_surviving_build\tn_surviving_lowcov1\tn_surviving_compress1\t"
        "n_surviving_lowcov2\tn_surviving_compress2\tn_surviving_tips\thap_indices\tis_traversal_limited\t"
        "is_graph_cycle\tis_graph_complex\tis_no_anchor\tis_short_anchor\tis_variant_in_anchor\t"
        "n_last_edge_kmers\tn_last_interior_kmers\tn_last_chain_gaps\t"
        "msa_shift_bp\tis_msa_exact_match\tis_msa_shifted\tis_msa_subsumed\t"
        "n_geno_true_alt_reads\tn_geno_total_ref_reads\tn_geno_reassigned_to_ref\t"
        "n_geno_reassigned_to_wrong_alt\tn_geno_non_overlapping\tis_geno_has_alt_support\tis_geno_no_overlap\n";
    out << HEADER_LINE;
    mHeaderWritten = true;
  }

  for (auto const& record : records) {
    auto const& probe = mVariants[record.mProbeId];
    auto const paths_str = fmt::format("{}", fmt::join(record.mHapIndices, ","));
    mWrittenProbeIds.insert(record.mProbeId);

    out << fmt::format(
        "{PROBE_ID}\t{CHROM}\t{POS}\t{REF}\t{ALT}\t{WINDOW}\t{N_RAW_ALT_READS}\t{KMER_SIZE}\t"
        "{N_EXPECTED_ALT_KMERS}\t{N_ALT_KMERS_IN_READS}\t{COMP_ID}\t{N_COMP_NODES}\t"
        "{N_SPLIT_ACROSS_COMPS}\t{N_SURVIVING_BUILD}\t{N_SURVIVING_LOWCOV1}\t"
        "{N_SURVIVING_COMPRESS1}\t{N_SURVIVING_LOWCOV2}\t{N_SURVIVING_COMPRESS2}\t"
        "{N_SURVIVING_TIPS}\t{HAP_INDICES}\t{IS_TRAVERSAL_LIMITED}\t{IS_GRAPH_CYCLE}\t"
        "{IS_GRAPH_COMPLEX}\t{IS_NO_ANCHOR}\t{IS_SHORT_ANCHOR}\t{IS_VARIANT_IN_ANCHOR}\t"
        "{N_LAST_EDGE}\t{N_LAST_INTERIOR}\t{N_LAST_GAPS}\t{MSA_SHIFT}\t{IS_MSA_EXACT}\t"
        "{IS_MSA_SHIFTED}\t{IS_MSA_SUBSUMED}\t{N_GENO_ALT}\t{N_GENO_REF}\t{N_GENO_REASSIGNED_REF}\t"
        "{N_GENO_REASSIGNED_ALT}\t{N_GENO_NON_OVL}\t{IS_GENO_ALT}\t{IS_GENO_NO_OVL}\n",
        fmt::arg("PROBE_ID", record.mProbeId),
        fmt::arg("CHROM", probe.mChrom),
        fmt::arg("POS", probe.mGenomeStart0),
        fmt::arg("REF", probe.mRef),
        fmt::arg("ALT", probe.mAlt),
        fmt::arg("WINDOW", record.mRegion.empty() ? "." : record.mRegion),
        fmt::arg("N_RAW_ALT_READS", probe.mRawAltCount),
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
        fmt::arg("IS_GRAPH_CYCLE", static_cast<int>(record.mIsGraphCycle)),
        fmt::arg("IS_GRAPH_COMPLEX", static_cast<int>(record.mIsGraphComplex)),
        fmt::arg("IS_NO_ANCHOR", static_cast<int>(record.mIsNoAnchor)),
        fmt::arg("IS_SHORT_ANCHOR", static_cast<int>(record.mIsShortAnchor)),
        fmt::arg("IS_VARIANT_IN_ANCHOR", static_cast<int>(record.mIsVariantInAnchor)),
        fmt::arg("N_LAST_EDGE", record.mLastEdgeCount),
        fmt::arg("N_LAST_INTERIOR", record.mLastInteriorCount),
        fmt::arg("N_LAST_GAPS", record.mLastGapCount),
        fmt::arg("MSA_SHIFT", record.mMsaShiftBp),
        fmt::arg("IS_MSA_EXACT", static_cast<int>(record.mIsMsaExactMatch)),
        fmt::arg("IS_MSA_SHIFTED", static_cast<int>(record.mIsMsaShifted)),
        fmt::arg("IS_MSA_SUBSUMED", static_cast<int>(record.mIsMsaSubsumed)),
        fmt::arg("N_GENO_ALT", record.mGenoTrueAltReads),
        fmt::arg("N_GENO_REF", record.mGenoTotalRefReads),
        fmt::arg("N_GENO_REASSIGNED_REF", record.mGenoReassignedToRef),
        fmt::arg("N_GENO_REASSIGNED_ALT", record.mGenoReassignedToWrongAlt),
        fmt::arg("N_GENO_NON_OVL", record.mGenoNonOverlapping),
        fmt::arg("IS_GENO_ALT", static_cast<int>(record.mIsGenoHasAltSupport)),
        fmt::arg("IS_GENO_NO_OVL", static_cast<int>(record.mIsGenoNoOverlap)));
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

}  // namespace lancet::cbdg
