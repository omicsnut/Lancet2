#include "lancet/caller/variant_bubble.h"

#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"

#include "absl/container/flat_hash_map.h"

#include <algorithm>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::caller {

auto CalculateVariantLength(std::string_view ref_allele, std::string_view alt_allele,
                            AlleleType vtype) -> i64 {
  if (vtype == AlleleType::SNV) {
    return 1;
  }

  auto const ref_len = static_cast<i64>(ref_allele.length());
  auto const alt_len = static_cast<i64>(alt_allele.length());
  auto const diff = alt_len - ref_len;

  // Net structural variance equals the string length difference for INS/DEL/CPX
  if (vtype == AlleleType::INS || vtype == AlleleType::DEL || vtype == AlleleType::CPX) {
    return diff;
  }

  // For strict MNPs (diff == 0), the biological length IS exactly the Sequence Core!
  // E.g., REF="ATGC", ALT="ACCC". Trimming `start=1` ('A') and `end=1` ('C') extracts
  // the mutation core `TG`->`CC` (length 4 - 1 - 1 = 2).
  usize start_match = 0;
  while (std::cmp_less(start_match, ref_len) &&
         std::cmp_less(start_match, alt_len) &&
         ref_allele[start_match] == alt_allele[start_match]) {
    start_match++;
  }

  usize end_match = 0;
  while (end_match < (ref_len - start_match) &&
         end_match < (alt_len - start_match) &&
         ref_allele[ref_len - 1 - end_match] == alt_allele[alt_len - 1 - end_match]) {
    end_match++;
  }

  return alt_len - static_cast<i64>(start_match) - static_cast<i64>(end_match);
}

namespace {

// Generic trim loop: `can_trim` checks whether the boundary character matches across all
// alleles (using string_view to avoid copies), `do_trim` performs the string mutation.
template <typename CanTrimFunc, typename DoTrimFunc>
void ApplyUnifiedTrim(VariantBubble& bubble, CanTrimFunc can_trim, DoTrimFunc do_trim) {
  // Guard: REF must retain at least 1 base (VCF requires non-empty REF)
  while (bubble.mRefAllele.length() > 1) {
    std::string_view const ref_view(bubble.mRefAllele);

    // If any allele's boundary character differs, trimming stops for all alleles
    bool const is_trimmable = std::ranges::all_of(bubble.mAltAllelesToHaps, [&](auto const& pair) {
      return can_trim(ref_view, std::string_view(pair.first));
    });

    if (!is_trimmable) break;

    // All alleles match — apply the trim
    do_trim(bubble.mRefAllele);

    // Flat hash map keys are const (mutating them would corrupt the hash index).
    // Rebuild the map with trimmed keys, using std::move on the haplotype vectors
    // to avoid copying.
    absl::flat_hash_map<std::string, std::vector<usize>> rehashed_map;
    rehashed_map.reserve(bubble.mAltAllelesToHaps.size());

    for (auto& [alt_seq, haps] : bubble.mAltAllelesToHaps) {
      std::string new_alt = alt_seq;
      do_trim(new_alt);
      rehashed_map.try_emplace(std::move(new_alt), std::move(haps));
    }

    bubble.mAltAllelesToHaps = std::move(rehashed_map);
  }
}

}  // namespace

void VariantBubble::NormalizeVcfParsimony() {
  if (mAltAllelesToHaps.empty() || mRefAllele.empty()) return;

  static constexpr auto CAN_MATCH_RIGHT = [](std::string_view refseq,
                                             std::string_view altseq) -> bool {
    return altseq.length() > 1 && altseq.back() == refseq.back();
  };

  static constexpr auto CAN_MATCH_LEFT = [](std::string_view refseq,
                                            std::string_view altseq) -> bool {
    return altseq.length() > 1 && altseq.front() == refseq.front();
  };

  static constexpr auto DO_TRIM_RIGHT = [](std::string& tseq) -> void { tseq.pop_back(); };
  static constexpr auto DO_TRIM_LEFT = [](std::string& tseq) -> void { tseq.erase(0, 1); };

  // Right Trim (eg. REF: "ATCG", ALTS: ["ACCG", "AGGG"] => "ATC", ["ACC", "AGG"])
  // Erases rightward bloat first allowing indels to
  // left-align against leftward structural boundaries.
  ApplyUnifiedTrim(*this, CAN_MATCH_RIGHT, DO_TRIM_RIGHT);

  usize const initial_ref_len = mRefAllele.length();
  // Left Trim (e.g. REF: "TTC", ALTS: ["TGC"] => "TC", ["GC"])
  ApplyUnifiedTrim(*this, CAN_MATCH_LEFT, DO_TRIM_LEFT);

  // After left-trimming, advance the genomic start position by the number of characters removed.
  mGenomeStartPos += (initial_ref_len - mRefAllele.length());
}

}  // namespace lancet::caller
