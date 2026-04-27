#ifndef SRC_LANCET_CALLER_VARIANT_BUBBLE_H_
#define SRC_LANCET_CALLER_VARIANT_BUBBLE_H_

#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"

#include "absl/container/flat_hash_map.h"

#include <string>
#include <string_view>
#include <vector>

namespace lancet::caller {

// ============================================================================
// CalculateVariantLength: biological mutation length after VCF padding removal.
//
// Companion to RawVariant::ClassifyVariant (which determines allele TYPE).
// For INS/DEL/CPX: returns alt_len − ref_len (net structural change).
// For MNP: trims shared prefix/suffix to extract the mutation core length.
// For SNV: always returns 1.
// ============================================================================
[[nodiscard]] auto CalculateVariantLength(std::string_view ref_allele, std::string_view alt_allele,
                                          AlleleType vtype) -> i64;

// ============================================================================
// VariantBubble: VCF string trimming and left-alignment container.
//
// Encapsulates per-bubble allele sequences extracted from the SPOA DAG,
// separate from graph traversal logic. Applies multi-allelic VCF parsimony
// trimming (right-trim then left-trim) to normalize the REF/ALT strings.
// https://genome.sph.umich.edu/wiki/Variant_Normalization
// ============================================================================
class VariantBubble {
 public:
  // ── 8B Align ────────────────────────────────────────────────────────────
  // Maps each distinct ALT allele string to the haplotype IDs that produced it.
  // O(1) lookup per allele via Abseil flat hash map.
  absl::flat_hash_map<std::string, std::vector<usize>> mAltAllelesToHaps;  // 8B (heap)
  std::string mRefAllele;            // 8B (heap pointer) + SSO
  std::vector<usize> mHapStarts;     // 8B (heap pointer)
  usize mGenomeStartPos = SIZE_MAX;  // 8B

  // Multi-allelic normalization trims shared prefix/suffix across ALL alleles simultaneously.
  // Trimming must stop if any ALT would become empty (losing anchor integrity for INDELs).
  void NormalizeVcfParsimony();
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_BUBBLE_H_
