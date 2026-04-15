#ifndef SRC_LANCET_CALLER_ALT_ALLELE_H_
#define SRC_LANCET_CALLER_ALT_ALLELE_H_

#include "lancet/base/types.h"

#include "absl/container/flat_hash_map.h"

#include <compare>
#include <string>
#include <utility>

namespace lancet::caller {

enum class AlleleType : i8 { REF = -1, SNV = 0, INS = 1, DEL = 2, MNP = 3, CPX = 4 };
enum class AlleleState : i8 { NONE = -1, SHARED = 0, CTRL = 1, CASE = 2, UNKNOWN = 3 };

// ===========================================================================
// DATA STRUCTURE CHOICE: `AltAllele` SUB-PAYLOAD
// ---------------------------------------------------------------------------
// A single genomic locus can mutate into multiple alternative forms across
// haplotypes or heterogeneous tumor populations (e.g., A→C on haplotype 1,
// A→T on haplotype 2).
//
// `AltAllele` packs each distinct ALT adjacent within a parent variant block.
// This maps 1:1 with VCF spec (ALT column holds "C,T") and prevents broken-
// apart biallelic records from misaligning or receiving conflicting parsimony
// trims later in the pipeline.
// ===========================================================================
struct AltAllele {
  std::string mSequence;

  // MULTI-ALLELIC LOCAL MATRIX MAP (ALTs):
  // A single multi-allelic locus has different string offsets depending on
  // which haplotype path a read traversed to reach it. A 100bp insertion
  // earlier in Haplotype 3 shifts this ALT's local index by +100 relative
  // to Haplotype 1.
  // Maps: Haplotype ID -> variant's exact local matrix start on that string.
  //
  // DATA STRUCTURE CHOICE: `absl::flat_hash_map`
  // Maps a supporting physical haplotype ID to its 0-indexed start position
  // within the variant. Abseil's FlatHashMap uses Google's SwissTable SIMD-
  // accelerated open-addressing implementation, storing keys/values
  // contiguously in memory. Outperforms std::map (O(logN) pointer chasing)
  // and std::unordered_map (cache-thrashing chunked linked lists).
  //
  // Moved up to satisfy 8B -> 4B -> 2B -> 1B alignment constraints.
  absl::flat_hash_map<usize, usize> mLocalHapStart0Idxs;

  i64 mLength = -1;
  AlleleType mType = AlleleType::REF;

  friend auto operator==(AltAllele const& lhs, AltAllele const& rhs) -> bool {
    return lhs.mSequence == rhs.mSequence && lhs.mLength == rhs.mLength && lhs.mType == rhs.mType;
  }

  // Lexicographic ordering by sequence ensures deterministic sort when the parent
  // vector of multiple ALTs is sorted within btree_set.
  friend auto operator<=>(AltAllele const& lhs, AltAllele const& rhs) {
    return lhs.mSequence <=> rhs.mSequence;
  }

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, AltAllele const& alt) -> HashState {
    // Use Abseil's hash combining
    return HashState::combine(std::move(hash_state), alt.mSequence, alt.mLength,
                              static_cast<i8>(alt.mType));
  }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_ALT_ALLELE_H_
