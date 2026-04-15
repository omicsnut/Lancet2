#ifndef SRC_LANCET_CALLER_COMBINED_SCORER_H_
#define SRC_LANCET_CALLER_COMBINED_SCORER_H_

#include "lancet/base/types.h"
#include "lancet/caller/allele_scoring_types.h"
#include "lancet/caller/genotyper.h"

#include "absl/types/span.h"

#include <vector>

namespace lancet::caller {

// ============================================================================
// Combined scorer: bridges local scoring → allele assignment.
//
// Combines the global mm2 alignment score, the local variant-region score
// (from local_scorer), soft-clip penalty, and edit distance into a single
// combined_score that determines which allele a read best supports.
//
// combined = (global_score - sc_penalty - local_raw_score) + (local_pbq_score * local_identity)
// ============================================================================

/// Score one read-haplotype alignment at a variant site.
/// Caller must pre-validate: alignment overlaps the variant region
/// (via Genotyper::ExtractHapBounds + Genotyper::OverlapsAlignment).
/// Returns a fully-populated ReadAlleleAssignment.
[[nodiscard]] auto ScoreReadAtVariant(Mm2AlnResult const& aln,
                                      absl::Span<u8 const> encoded_haplotype,
                                      ReadAlnContext const& read_ctx,
                                      HapVariantBounds const& bounds) -> ReadAlleleAssignment;

/// Compute edit distance (NM) of a read against one specific haplotype.
/// Per SAM spec: excludes soft/hard clips from NM count.
/// Falls back to full read length as worst-case NM if no valid alignment exists.
[[nodiscard]] auto ComputeHaplotypeEditDistance(std::vector<Mm2AlnResult> const& alns,
                                                absl::Span<u8 const> encoded_haplotype,
                                                absl::Span<u8 const> qry_seq_encoded,
                                                usize qry_read_length, usize hap_idx) -> u32;

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_COMBINED_SCORER_H_
