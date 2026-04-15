#include "lancet/caller/raw_variant.h"

#include "lancet/caller/alt_allele.h"

namespace lancet::caller {

// =========================================================================================
// STRICT SEQUENCE CORE CLASSIFICATION ENGINE
// =========================================================================================
// WHY DO WE SQUEEZE THE STRINGS AGAIN IF THE `VariantExtractor` FSM ALREADY APPLIES PARSIMONY?
//
// While `VariantBubble::NormalizeVcfParsimony` performs a simultaneous VT-style right-trim
// and left-alignment, it MUST apply the same trim across ALL alleles in a multi-allelic
// bubble block. Therefore, the global trimmer is restricted by the bounds of the WIDEST
// allele in the cluster.
//
// -> THE MULTI-ALLELIC SHIELDING PROBLEM:
// Consider a bubble traversing 3 paths:
//    [REF]     :  A T G C
//    [ALT1]    :  A T
//    [ALT2]    :  A G G C
//
// Let's trace Global VCF Parsimony on ALT1 vs REF:
// It WANTS to right-trim `G` and `C` to isolate the `GC` deletion completely (ATGC -> AT).
// However, `ALT2` completely blocks the `C` and `G` from being universally erased, because
// `ALT2` actively utilizes them for its own structural integrity.
//
// Thus, VCF Parsimony halts prematurely for ALT1!
// When ALT1 evaluates its payload conventionally: `REF="ATGC"`, `ALT="AT"`.
// If we naively utilized length logic (`diff = -2` and `length > 1`), we would classify this
// erroneously as a `CPX` (Complex) mutation because the length boundary is artificially inflated!
//
// -> THE SEQUENCE CORE SOLUTION:
// By aggressively symmetrically squeezing matching 5' prefixes and 3' suffixes exclusively
// between the 1-on-1 pairs immediately prior to classification, we computationally decouple
// the biological Core from the VCF-Padding constraints.
//
//    [REF]     :  A T (G C)  ---> "GC"
//    [ALT1]    :  A T        ---> ""     ===> Result: Pure `DEL`!
//
// This executes in O(N) using forward and reverse bounds pointers with zero memory allocation
// overhead.
// =========================================================================================
auto RawVariant::ClassifyVariant(std::string_view ref_seq, std::string_view alt_seq) -> AlleleType {
  usize start_match = 0;

  while (start_match < ref_seq.length() &&
         start_match < alt_seq.length() &&
         ref_seq[start_match] == alt_seq[start_match]) {
    start_match++;
  }

  if (start_match == ref_seq.length() && start_match == alt_seq.length()) {
    return AlleleType::REF;
  }

  // Ensures bounding pointers do not overlap or double-count characters
  usize end_match = 0;
  while (end_match < (ref_seq.length() - start_match) &&
         end_match < (alt_seq.length() - start_match) &&
         ref_seq[ref_seq.length() - 1 - end_match] == alt_seq[alt_seq.length() - 1 - end_match]) {
    end_match++;
  }

  auto const ref_core_len = ref_seq.length() - start_match - end_match;
  auto const alt_core_len = alt_seq.length() - start_match - end_match;

  if (ref_core_len == 0 && alt_core_len > 0) return AlleleType::INS;
  if (ref_core_len > 0 && alt_core_len == 0) return AlleleType::DEL;
  if (ref_core_len == 0 || alt_core_len == 0) return AlleleType::REF;

  // Mixed/Complex concurrent INS/DEL
  if (ref_core_len != alt_core_len) return AlleleType::CPX;

  // SNV/MNP classification
  return (ref_core_len == 1) ? AlleleType::SNV : AlleleType::MNP;
}

}  // namespace lancet::caller
