#ifndef SRC_LANCET_CORE_PROBE_DIAGNOSTICS_H_
#define SRC_LANCET_CORE_PROBE_DIAGNOSTICS_H_

#include "lancet/base/types.h"
#include "lancet/caller/genotyper.h"
#include "lancet/caller/variant_set.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/probe_results_writer.h"
#include "lancet/cbdg/probe_tracker.h"
#include "lancet/core/window.h"

#include <filesystem>
#include <memory>

namespace lancet::core {

// ============================================================================
// ProbeDiagnostics: post-graph analysis for missed variant probes.
//
// Composes cbdg::ProbeTracker (graph-stage k-mer tracking) with MSA and
// Genotyper analysis methods that require caller-layer types.
//
// Ownership:
//   VariantBuilder owns ProbeDiagnostics.
//   ProbeDiagnostics owns ProbeTracker (submits records to shared writer).
//   Graph holds a non-owning ProbeTracker* via Tracker().
// ============================================================================
class ProbeDiagnostics {
 public:
  ProbeDiagnostics() = default;

  /// Initialize the diagnostics subsystem: load truth variants, set the
  /// shared results writer, and attach the precomputed probe index.
  void Initialize(std::filesystem::path const& variants_path,
                  std::shared_ptr<cbdg::ProbeResultsWriter> results_writer,
                  std::shared_ptr<cbdg::ProbeIndex const> probe_index);

  /// Return a pointer to the internal ProbeTracker for Graph to use.
  /// Returns nullptr when probe tracking is not active (no --probe-variants).
  [[nodiscard]] auto Tracker() -> cbdg::ProbeTracker* {
    return mTracker.IsActive() ? &mTracker : nullptr;
  }

  /// Flush accumulated probe records to the shared writer. No-op if inactive.
  void SubmitCompleted() {
    if (mTracker.IsActive()) mTracker.SubmitCompleted();
  }

  /// Check whether each probe variant was extracted by the MSA.
  /// Three match tiers (highest priority first):
  ///   1. Exact match — same position and alleles
  ///   2. Shifted match — same alleles but different position
  ///   3. Representation match — variant subsumed by larger MNV
  void CheckMsaExtraction(caller::VariantSet const& variant_set, Window const& window,
                          usize component_idx);

  /// Check whether each probe variant received correct genotyper read support.
  /// Tracks stolen reads (ALT-carrying reads misassigned to REF or wrong ALT),
  /// non-overlapping reads, and total ALT/REF coverage.
  void CheckGenotyperResult(caller::Genotyper::Result const& genotyped,
                            caller::VariantSet const& variant_set, usize component_idx);

 private:
  cbdg::ProbeTracker mTracker;

  /// Check if a probe variant is represented as a sub-mutation within a larger
  /// same-length substitution (MNV representation match).
  static void CheckMsaRepresentationMatch(caller::VariantSet const& variant_set,
                                          cbdg::ProbeVariant const& probe,
                                          cbdg::ProbeKRecord& record);
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_PROBE_DIAGNOSTICS_H_
