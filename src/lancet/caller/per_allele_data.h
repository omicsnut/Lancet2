#ifndef SRC_LANCET_CALLER_PER_ALLELE_DATA_H_
#define SRC_LANCET_CALLER_PER_ALLELE_DATA_H_

#include "lancet/base/types.h"

#include "absl/container/flat_hash_map.h"

#include <vector>

namespace lancet::caller {

enum class Strand : bool { FWD, REV };

// ============================================================================
// PerAlleleData: per-allele evidence storage for a single variant site.
//
// Aggregates all read-level measurements for one specific allele (REF, ALT1,
// ALT2, ...) at a variant position. Used by VariantSupport to compute
// aggregate metrics for VCF FORMAT fields.
//
// Each vector stores one entry per read (not per-base). For indels, base
// quality is the MINIMUM PBQ across the variant region, ensuring that each
// read's evidence is represented by a single conservative value.
// ============================================================================
struct PerAlleleData {
  // ── 8B Align ────────────────────────────────────────────────────────────
  // All members are 8B-aligned (std::vector = {ptr, size, capacity} = 24B each,
  // flat_hash_map and usize are also 8B-aligned).

  // Deduplicate reads: a read can support an allele only once per strand.
  // Key = read name hash, value = strand seen.
  absl::flat_hash_map<u32, Strand> mNameHashes;  // 8B aligned (heap pointer)

  // Per-read representative base quality at the variant position, split by
  // strand. For indels, this is the MINIMUM PBQ across the variant region
  // (one entry per read, NOT per base position -- see AddEvidence comment).
  std::vector<u8> mFwdBaseQuals;  // 8B aligned (24B)
  std::vector<u8> mRevBaseQuals;  // 8B aligned (24B)

  // Per-read mapping quality (for RMS RMQ computation).
  std::vector<u8> mMapQuals;  // 8B aligned (24B)

  // Per-read normalized alignment score (for filtering / annotation).
  std::vector<f64> mAlnScores;  // 8B aligned (24B)

  // Insert sizes from properly-paired reads (for FLD FORMAT tag).
  // Only non-zero insert sizes from proper pairs are tracked.
  std::vector<f64> mProperPairIsizes;  // 8B aligned (24B)

  // Folded read positions: min(p, 1-p) for RPCD effect size test.
  std::vector<f64> mFoldedReadPositions;  // 8B aligned (24B)

  // Edit distances (NM) against REF haplotype for ASMD delta.
  // Stored as f64 for mean computation.
  std::vector<f64> mRefNmValues;  // 8B aligned (24B)

  // Edit distances (NM) against each read's assigned haplotype for AHDD.
  // Stored as f64 for mean computation. Populated for ALL reads — REF reads
  // are scored against haplotype 0 (the REF), ALT reads against their
  // winning ALT haplotype.
  std::vector<f64> mOwnHapNmValues;  // 8B aligned (24B)

  // Fragment alignment start positions for FSSE (Fragment Start Shannon Entropy).
  // Genomic coordinates — high repeat count at the same position signals PCR duplicates.
  std::vector<i64> mAlignmentStarts;  // 8B aligned (24B)

  // SPOA path IDs for HSE (Haplotype Segregation Entropy).
  // Tracks which assembled haplotype each ALT read was assigned to.
  std::vector<u32> mHaplotypeIds;  // 8B aligned (24B)

  // Count of soft-clipped reads supporting this allele (for SCA FORMAT tag).
  usize mSoftClipCount = 0;  // 8B
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_PER_ALLELE_DATA_H_
