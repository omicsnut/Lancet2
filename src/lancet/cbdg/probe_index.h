#ifndef SRC_LANCET_CBDG_PROBE_INDEX_H_
#define SRC_LANCET_CBDG_PROBE_INDEX_H_

#include "lancet/base/types.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <filesystem>
#include <string>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// ProbeVariant: a truth variant loaded from the missed_variants.txt file.
//
// Each entry represents a variant the pipeline failed to call. The probe
// system generates ALT-unique k-mers for each variant and tracks their
// survival through graph construction, pruning, and path enumeration to
// identify exactly where the variant signal is lost.
// ============================================================================
struct ProbeVariant {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string mChrom;       // 32B
  std::string mRef;         // 32B
  std::string mAlt;         // 32B
  usize mGenomeStart0 = 0;  // 8B — 0-based genomic start position

  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mProbeId = 0;      // 2B — index into the variants vector
  u16 mRawAltCount = 0;  // 2B — raw (unfiltered) ALT read count from truth concordance
};

// ============================================================================
// ProbeIndex: precomputed global ALT k-mer lookup for all probe variants.
//
// Built once at pipeline startup. For each k in [minK..maxK, stepK] and
// each probe variant, fetches the flanking reference context from
// hts::Reference, builds an ALT context, hashes canonical ALT-unique
// k-mers, and stores them indexed by k-value and canonical hash.
//
// Thread safety: immutable after construction. Shared across all worker
// threads via std::shared_ptr<ProbeIndex const>.
// ============================================================================
class ProbeIndex {
 public:
  struct KmerEntry {
    // ── 2B Align ──────────────────────────────────────────────────────────
    u16 mProbeId = 0;           // 2B — which probe variant
    u16 mAltContextOffset = 0;  // 2B — position in ALT context k-mer chain
  };

  using KmerToProbes = absl::flat_hash_map<u64, absl::InlinedVector<KmerEntry, 2>>;

  /// Parse missed_variants.txt — loads all probe variants with their raw ALT read count.
  [[nodiscard]] static auto LoadVariantsFromFile(std::filesystem::path const& path)
      -> std::vector<ProbeVariant>;

  /// Retrieve the precomputed k-mer→probe map for a specific kmer size.
  /// Returns nullptr if no index was built for that k.
  [[nodiscard]] auto ForKmerSize(usize kmer_size) const -> KmerToProbes const*;

  /// Build the full index across all k-values in [min_kmer_len..max_kmer_len, step].
  /// Fetches flanking reference context from the reference FASTA for each probe.
  [[nodiscard]] static auto Build(absl::Span<ProbeVariant const> variants,
                                  hts::Reference const& reference, usize min_kmer_len,
                                  usize max_kmer_len, usize kmer_step_len) -> ProbeIndex;

  /// True when the index has no entries (no variants loaded or no k-values built).
  [[nodiscard]] auto IsEmpty() const noexcept -> bool;

  /// Total number of (kmer_hash, probe_entry) pairs across all k-values.
  [[nodiscard]] auto TotalEntries() const noexcept -> usize;

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  absl::flat_hash_map<usize, KmerToProbes> mPerKIndex;  // keyed by kmer_size
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PROBE_INDEX_H_
