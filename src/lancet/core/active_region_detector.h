#ifndef SRC_LANCET_CORE_ACTIVE_REGION_DETECTOR_H_
#define SRC_LANCET_CORE_ACTIVE_REGION_DETECTOR_H_

#include "lancet/core/read_collector.h"
#include "lancet/hts/reference.h"

#include <filesystem>

namespace lancet::core {

/// Returns true if the given BAM/CRAM file contains MD tags.
/// Peeks at the first 1000 reads to check for the MD auxiliary field.
/// Active region detection depends on MD tags — if absent, the caller
/// must disable active region filtering to avoid false negatives.
[[nodiscard]] auto HasMdTag(std::filesystem::path const& aln_path,
                            std::filesystem::path const& ref_path) -> bool;

/// Convenience negation of HasMdTag — avoids double negatives at call sites.
[[nodiscard]] inline auto IsMissingMdTag(std::filesystem::path const& aln_path,
                                         std::filesystem::path const& ref_path) -> bool {
  return !HasMdTag(aln_path, ref_path);
}

/// Returns true if the region shows evidence of mutations (≥2 reads with
/// mismatches, insertions, deletions, or soft-clips at the same position).
/// Uses a lightweight MD tag + CIGAR prescan to avoid full assembly of
/// inactive regions.
[[nodiscard]] auto IsActiveRegion(ReadCollector::Params const& params,
                                  hts::Reference::Region const& region) -> bool;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_ACTIVE_REGION_DETECTOR_H_
