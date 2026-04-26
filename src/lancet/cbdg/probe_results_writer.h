#ifndef SRC_LANCET_CBDG_PROBE_RESULTS_WRITER_H_
#define SRC_LANCET_CBDG_PROBE_RESULTS_WRITER_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/probe_index.h"

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_set.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"

#include <filesystem>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

struct ProbeKRecord;

// ============================================================================
// ProbeResultsWriter: thread-safe, append-only TSV writer for probe records.
//
// One row per (probe_id, kmer_size) pair. Records arrive from worker threads
// via Append() and are written to disk immediately under mutex protection.
//
// Semi-realtime chunked write architecture:
//
//   Thread 0 ──ProcessWindow──▶ SubmitCompleted() ─┐
//   Thread 1 ──ProcessWindow──▶ SubmitCompleted() ─┤
//   Thread 2 ──ProcessWindow──▶ SubmitCompleted() ─┤
//      ...                                         │
//   Thread N ──ProcessWindow──▶ SubmitCompleted() ─┤
//                                                  ▼
//                                    ┌──────────────────────┐
//                                    │ ProbeResultsWriter   │
//                                    │                      │
//                                    │ absl::MutexLock(mu)  │
//                                    │ write header (once)  │
//                                    │ append TSV lines     │
//                                    │ flush to disk        │
//                                    └──────────────────────┘
//                                              │
//                                              ▼
//                                    probe_results.tsv
//
// Ownership:
//   PipelineRunner creates one shared instance at startup.
//   Each ProbeTracker holds a shared_ptr and submits after each window.
//   The writer outlives all ProbeTrackers.
//
// Thread safety:
//   Append() acquires absl::Mutex, serializing all file writes and
//   guaranteeing exactly one header line.
//
// Output schema (probe_results.tsv — 39 columns):
//
// ┌────────────────────────────────────────────────────────────────────┐
// │ §1 Identity         probe_id chrom pos ref alt                     │
// │ §2 Read evidence    n_tier1_reads n_alt_kmers_in_reads             │
// │ §3 Graph assembly   kmer_size n_expected_alt_kmers comp_id         │
// │                     n_comp_nodes n_split_across_comps              │
// │ §4 Survival funnel  n_surviving_{build,lowcov1,compress1,          │
// │                                  lowcov2,compress2,tips}           │
// │ §5 Path & anchor    hap_indices is_traversal_limited               │
// │                     is_cycle_retry is_complex_retry                │
// │                     is_no_anchor is_variant_in_anchor              │
// │ §6 Offset analysis  n_last_edge_kmers n_last_interior_kmers        │
// │                     n_last_chain_gaps                              │
// │ §7 MSA extraction   msa_shift_bp is_msa_exact_match                │
// │                     is_msa_shifted is_msa_representation           │
// │ §8 Genotyper        n_geno_true_alt_reads n_geno_total_ref_reads   │
// │                     n_geno_stolen_to_ref n_geno_stolen_to_wrong_alt│
// │                     n_geno_non_overlapping is_geno_has_alt_support │
// │                     is_geno_no_result                              │
// │ §9 Attribution      lost_at_stage                                  │
// └────────────────────────────────────────────────────────────────────┘
//
// Column reference:
//
//   §1 Identity
//     probe_id              u16  [0, N)    index into the input variants file
//     chrom                 str            chromosome name
//     pos                   u64            0-based genome start position
//     ref, alt              str            REF and ALT allele sequences
//
//   §2 Read evidence — raw signal before graph construction.
//     n_tier1_reads         u16  ≥1        ALT reads reported by tier-1 caller
//     n_alt_kmers_in_reads  u16  ≥0        ALT-unique k-mers found in raw reads
//
//   §3 Graph assembly — graph topology at the current k.
//     kmer_size             usize          odd, from --min-kmer to --max-kmer
//                                          by --kmer-step (defaults: 13–127, step 6)
//     n_expected_alt_kmers  u16  ≥0        ALT-unique k-mers at this k
//     comp_id               usize          connected component containing majority
//                                          of this probe's tagged nodes
//     n_comp_nodes          usize          nodes in that component
//     n_split_across_comps  u16  ≥0        0 = single component; >0 = fragmented
//
//   §4 Survival funnel — k-mer attrition through pruning.
//     n_surviving_build      u16  [0, n_expected]  after graph construction
//     n_surviving_lowcov1    u16                    after 1st low-cov removal
//     n_surviving_compress1  u16                    after 1st unitig compression
//     n_surviving_lowcov2    u16                    after 2nd low-cov removal
//     n_surviving_compress2  u16                    after 2nd unitig compression
//     n_surviving_tips       u16                    after tip removal (final)
//
//     Monotonically non-increasing: build ≥ lowcov1 ≥ … ≥ tips. First drop
//     from >0 to 0 identifies the destructive stage. See lost_at_stage (§9).
//
//   §5 Path & anchor — traversal and structural failure modes.
//     hap_indices           str  "0,1,3"|"."  paths containing variant ("."=none)
//     is_traversal_limited  0|1               BFS budget exhausted
//     is_cycle_retry        0|1               cycle detected → incremented k, retried
//     is_complex_retry      0|1               graph too complex → incremented k, retried
//     is_no_anchor          0|1               no valid source/sink pair
//     is_variant_in_anchor  0|1               variant overlaps anchor range
//
//   §6 Offset analysis — positional profile of surviving k-mers.
//     n_last_edge_kmers      u16  [0, 2]          surviving at chain edges
//     n_last_interior_kmers  u16  [0, n_exp−2]    surviving in chain interior
//     n_last_chain_gaps      u16  ≥0              contiguous offset gaps
//
//   §7 MSA extraction — did the VariantSet contain this truth variant?
//     msa_shift_bp            i16               signed bp shift (0 = exact)
//     is_msa_exact_match      0|1               position + alleles match
//     is_msa_shifted          0|1               same alleles, shifted position
//     is_msa_representation   0|1               subsumed by larger MNV
//
//   §8 Genotyper — read assignment at the truth variant site.
//     n_geno_true_alt_reads       u16  ≥0  reads on truth ALT allele
//     n_geno_total_ref_reads      u16  ≥0  total reads on REF
//     n_geno_stolen_to_ref        u16  ≥0  ALT reads misassigned to REF
//     n_geno_stolen_to_wrong_alt  u16  ≥0  ALT reads misassigned to wrong ALT
//     n_geno_non_overlapping      u16  ≥0  ALT reads outside alignment span
//     is_geno_has_alt_support     0|1      ≥1 read supports truth ALT
//     is_geno_no_result           0|1      variant absent from result map
//
//   §9 Attribution — first pipeline stage where the variant was lost.
//     lost_at_stage  str  Priority-ordered cascade (first match wins):
//
//       Not processed (highest priority):
//         "not_processed"        probe never activated in any window
//
//       Structural:
//         "variant_in_anchor"    variant overlaps source/sink k-mer range
//         "no_anchor"            no valid source/sink pair
//         "short_anchor"         anchor region too short
//         "cycle_retry"          cycle detected, incremented k, retried
//         "complex_retry"        graph too complex, incremented k, retried
//
//       Pruning (checked in pipeline order):
//         "pruned_at_build"      ALT k-mers absent after construction
//         "pruned_at_lowcov1"    lost during 1st low-cov removal
//         "pruned_at_compress1"  lost during 1st unitig compression
//         "pruned_at_lowcov2"    lost during 2nd low-cov removal
//         "pruned_at_compress2"  lost during 2nd unitig compression
//         "pruned_at_tips"       lost during tip removal
//
//       Path enumeration:
//         "traversal_limited"    BFS exhausted, variant not in any path
//         "no_path"              variant absent from enumerated haplotypes
//
//       MSA:
//         "msa_no_variant"       not found in VariantSet
//         "msa_shifted"          found but position-shifted
//         "msa_representation"   subsumed by larger MNV
//
//       Genotyper:
//         "geno_no_result"       absent from genotyper result map
//         "geno_no_support"      0 reads on truth ALT
//         "geno_stolen"          ALT reads misassigned to REF/wrong ALT
//
//       Success:
//         "survived"             variant found and genotyped with ALT support
// ============================================================================
class ProbeResultsWriter {
 public:
  ProbeResultsWriter(std::filesystem::path results_path, std::vector<ProbeVariant> variants);

  /// Thread-safe: format and write a batch of records to disk.
  /// Called from ProbeTracker::SubmitCompleted() after each window.
  void Append(absl::Span<ProbeKRecord const> records);

  /// Emit default "not_processed" records for any loaded variant that never
  /// appeared in any thread's output. Call exactly once after all workers finish.
  void EmitUnprocessedProbes();

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::filesystem::path mResultsPath;
  std::vector<ProbeVariant> mVariants;
  absl::flat_hash_set<u16> mWrittenProbeIds ABSL_GUARDED_BY(mMutex);
  mutable absl::Mutex mMutex;

  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mHeaderWritten ABSL_GUARDED_BY(mMutex) = false;

  /// Derive the lost_at attribution string from a record's state.
  /// Priority-ordered cascade: structural → pruning → path → MSA → genotyper → survived.
  [[nodiscard]] static auto DeriveLostAt(ProbeKRecord const& record) -> std::string_view;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PROBE_RESULTS_WRITER_H_
