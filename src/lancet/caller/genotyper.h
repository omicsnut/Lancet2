#ifndef SRC_LANCET_CALLER_GENOTYPER_H_
#define SRC_LANCET_CALLER_GENOTYPER_H_

#include "lancet/base/types.h"
#include "lancet/caller/allele_scoring_types.h"
#include "lancet/caller/scoring_constants.h"
#include "lancet/caller/support_array.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/read.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/cigar_utils.h"

extern "C" {
#include "minimap.h"
}

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::caller {

class VariantSet;
class RawVariant;
// ============================================================================
// Alignment result from mm_map for a single read-to-haplotype alignment
// ============================================================================
struct Mm2AlnResult {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<hts::CigarUnit> mCigar;  // 24B — CIGAR operations
  f64 mIdentity = 0.0;                 // 8B  — gap-compressed identity
  usize mHapIdx = 0;                   // 8B  — index of haplotype this alignment maps to
  // ── 4B Align ────────────────────────────────────────────────────────────
  i32 mScore = 0;     // 4B  — DP alignment score
  i32 mRefStart = 0;  // 4B  — 0-based start on haplotype
  i32 mRefEnd = 0;    // 4B  — 0-based end on haplotype
};

// ============================================================================
// ReadAlleleAssignment: per-read allele assignment result.
//
// Produced by ScoreReadAtVariant, consumed by AddToTable → ReadEvidence.
// ============================================================================
// Scoring components for allele assignment
// ============================================================================
//
// Each read-haplotype pair produces multiple independent signals. The combined
// score integrates them to assign the read to its best-matching allele:
//
//   combined = (global_score - local_raw_score - sc_penalty) + (local_pbq_score *
//   local_identity)
//
// Components:
//
//   global_score:    mm_map DP score of the full read→haplotype alignment.
//                    Captures how well the entire read fits this haplotype,
//                    including flanking contexts.
//
//   local_raw_score: The raw substitution matrix score of the variant overlap
//                    slice. We SUBTRACT this from global_score to prevent
//                    double-counting the variant when we add the PBQ score.
//
//   sc_penalty:      Penalization of soft-clipped read tails. Prevents noisy
//                    supplementary mappings from inflating their assignment scores
//                    over clean end-to-end alignments.
//
//   local_pbq_score: PBQ-weighted DP score within the variant region only.
//                    Scales substitution scores by Phred confidence
//                    (1 - 10^(-PBQ/10)), analogous to GATK's local PairHMM.
//
//   local_identity:  Fraction of exact matches inside the variant region.
//                    Acts as a confidence gate on the local PBQ score. A high
//                    local_score from a noisy alignment (low identity) is
//                    discounted, while one from a clean alignment is trusted.
//
//   - Why subtract `local_raw_score`?
//     global_score already includes the raw matrix alignment cost spanning the
//     variant sub-region. Naively appending local_pbq_score would double-count
//     the variant. Subtracting local_raw_score carves a "hole" out of the global
//     alignment path, allowing us to drop the PBQ-weighted score into that
//     specific locus.
//
//   - Why subtract `sc_penalty`?
//     Soft-clipped read tails are unaligned sequence. Minimap2 exempts them
//     from the primary DP score. Subtracting a soft-clip penalty suppresses
//     the combined global scores of chimeric supplementary mappings, preventing
//     partially-matching noise reads from winning allele assignments against
//     clean end-to-end alignments.
//
//   - Why (local_pbq_score * local_identity)?
//     This penalizes "lucky" alignments in low-complexity regions. A noisy
//     chimeric read might accumulate a high DP score by traversing a repetitive
//     STR locus. Multiplying by the exact-match fraction (identity) acts as a
//     confidence gate: perfect alignments (identity = 1.0) retain full PBQ weight;
//     fragmented, heavily-gapped alignments (e.g. identity < 0.7) are discounted,
//     suppressing repeat-region artifacts and structural chimeras.
//
// ============================================================================
// Folded read position
// ============================================================================
// Folded read position: min(p, 1−p) where p = variant_query_pos / read_length.
// 0.0 = variant at read edge, 0.5 = variant at read center.
//
// Why fold? Artifacts from 3' quality degradation and 5' soft-clip
// misalignment cluster at BOTH read ends. With raw positions, these
// bimodal clusters (e.g. positions 5 and 145 in a 150bp read) average
// to ~75 — indistinguishable from a true centered variant. Folding maps
// both ends to the same low-value space, converting the bimodal trap
// into a unidirectional signal: "Are ALT alleles systematically closer
// to read edges than REF alleles?" Used for RPCD FORMAT field.
//
// ============================================================================
// Edit distance (NM)
// ============================================================================
// Edit distance (NM) of this read against the REF haplotype (hap_idx=0).
// Mismatches (under M ops, comparing query vs encoded REF) + insertion
// bases + deletion bases. Soft clips, hard clips, N-skips excluded per
// SAM spec. Always measured against the reference haplotype regardless
// of allele assignment so that ASMD = mean(ALT NM) − mean(REF NM)
// cancels the variant's own contribution and isolates excess noise.
//
// ============================================================================
// Representative base quality
// ============================================================================
// Representative base quality at this variant for this read.
//
// For SNVs: the single base quality at the variant position.
// For indels: the MINIMUM base quality across all read positions spanning
//             the variant region (weakest-link summary).
//
// Why minimum (not mean/median)?
//   The confidence in observing a complete indel is bounded by the least
//   confident base in the region. A 10bp deletion where 9 bases are Q30
//   and 1 base is Q5 should not be treated as high-confidence.
//
// How other callers handle this:
//   - GATK HaplotypeCaller: PairHMM produces a single per-read likelihood
//     that integrates all base qualities across the alignment. One read
//     always contributes exactly one observation regardless of variant size.
//   - bcftools mpileup: uses the minimum base quality in the indel region
//     as the representative quality for that read.
//
// We follow the bcftools convention since we don't have a full PairHMM.
// The collapse to a single value here ensures that downstream PL and PBQ
// computations correctly treat each read as one independent observation.
struct ReadAlleleAssignment {
  // ── 8B Align ────────────────────────────────────────────────────────────
  f64 mLocalScore = 0.0;     // PBQ-weighted DP score within variant region
  f64 mLocalIdentity = 0.0;  // fraction of exact matches in variant region
  f64 mFoldedReadPos = 0.0;  // Used for RPCD FORMAT field

  // ── 4B Align ────────────────────────────────────────────────────────────
  i32 mGlobalScore = 0;                  // Adjusted score: mm_map DP − sc_penalty − local_raw_score
  u32 mRefNm = 0;                        // Edit distance (NM) against REF haplotype (for ASMD)
  u32 mOwnHapNm = 0;                     // Edit distance (NM) against assigned haplotype (for AHDD)
  u32 mAssignedHaplotypeId = 0;          // SPOA path index this read was assigned to (for HSE)
  AlleleIndex mAllele = REF_ALLELE_IDX;  // Allele enumerator, defaults to REF

  // ── 1B Align ────────────────────────────────────────────────────────────
  u8 mBaseQualAtVar = 0;  // Representative base quality at this variant for this read

  [[nodiscard]] auto CombinedScore() const -> f64 {
    return static_cast<f64>(mGlobalScore) + (mLocalScore * mLocalIdentity);
  }
};

/*
 * ============================================================================
 * ALIGNMENT PARADIGM SHIFT: MSA vs. Read-to-Haplotype Mapping
 * ============================================================================
 * There is an intentional paradox in the pipeline's alignment parameters:
 * We use highly FORGIVING parameters to build the Contig-Reference MSA, but
 * highly STRICT parameters when mapping raw reads back to the haplotypes.
 *
 * This is because the expected source of sequence divergence fundamentally
 * shifts between Variant Discovery and Genotyping.
 *
 * 1. THE MSA (Contig vs. Reference) -> Modeling BIOLOGY (The Mapmaker)
 *    - Assumption: The contig is high-confidence. Divergence is true mutation.
 *    - Strategy: FORGIVING. We use cheap gap extensions and high mismatch
 *      tolerance to force the algorithm to stretch across massive biological
 *      indels and dense MNVs without soft-clipping.
 *    - Goal: Discover and catalog the variant by keeping the alignment intact.
 *
 * 2. READ-TO-HAPLOTYPE -> Modeling PHYSICS (The GPS)
 *    - Assumption: Raw reads are noisy, but the biological variants are ALREADY
 *      baked into the assembled haplotype sequences.
 *    - Divergence: Sequencer error, adapter garbage, or chimeric artifacts.
 *    - Strategy: STRICT. A read should perfectly match its parent haplotype
 *      with zero biological gaps. We use heavy gap/mismatch penalties to
 *      punish sequencer noise, force soft-clipping of garbage read-tails, and
 *      prevent "allele bleeding" (noisy reads mapping to the wrong allele).
 *    - Goal: Accurately segregate read support to calculate clean VAFs.
 * ============================================================================
 */
// ============================================================================
// Genotyper: minimap2-based read-to-haplotype alignment for genotyping
//
// Scoring parameters are custom for Illumina read-to-contig realignment,
// NOT the standard 'sr' preset. See scoring_constants.h for scoring values.
//
// Isolation boundary: AssignReadToAlleles() encapsulates the alignment
// engine. Everything downstream (AddToTable, VariantSupport) is decoupled.
// ============================================================================
class Genotyper {
 public:
  Genotyper();

  using Reads = absl::Span<cbdg::Read const>;
  using Haplotypes = absl::Span<std::string const>;
  using Result = absl::flat_hash_map<RawVariant const*, SupportArray>;

  [[nodiscard]] auto Genotype(Haplotypes hap_seqs, Reads qry_reads, VariantSet const& variant_set)
      -> Result;

 private:
  // ============================================================================
  // Minimap2 RAII wrappers
  // ============================================================================
  struct MmIdxDeleter {
    void operator()(mm_idx_t* idx) noexcept { mm_idx_destroy(idx); }
  };

  struct MmTbufDeleter {
    void operator()(mm_tbuf_t* tbuf) noexcept { mm_tbuf_destroy(tbuf); }
  };

  using MappingOpts = std::unique_ptr<mm_mapopt_t>;
  using IndexingOpts = std::unique_ptr<mm_idxopt_t>;
  using ThreadBuffer = std::unique_ptr<mm_tbuf_t, MmTbufDeleter>;
  using Minimap2Index = std::unique_ptr<mm_idx_t, MmIdxDeleter>;

  static constexpr usize REF_HAP_IDX = 0;

  // ============================================================================
  // Outer Class Variables Block (Sorted by descending size: 24B -> 8B -> 4B)
  // ============================================================================
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<Minimap2Index> mIndices;  // 24B
  // numeric-encoded haplotypes for local scoring
  std::vector<std::vector<u8>> mEncodedHaplotypes;               // 24B
  MappingOpts mMappingOpts = std::make_unique<mm_mapopt_t>();    // 8B
  IndexingOpts mIndexingOpts = std::make_unique<mm_idxopt_t>();  // 8B
  ThreadBuffer mThreadBuffer = ThreadBuffer(mm_tbuf_init());     // 8B

  void ResetData(Haplotypes hap_seqs);

  using PerVariantAssignment = absl::flat_hash_map<RawVariant const*, ReadAlleleAssignment>;
  [[nodiscard]] auto AssignReadToAlleles(cbdg::Read const& qry_read, VariantSet const& variant_set)
      -> PerVariantAssignment;

  [[nodiscard]] auto AlignToAllHaplotypes(cbdg::Read const& qry_read) -> std::vector<Mm2AlnResult>;

  // ============================================================================
  // ExtractHapBounds: resolve a minimap2 alignment's haplotype index to the
  // variant's physical coordinates on that specific haplotype.
  //
  // Each assembled haplotype carries the variant at a DIFFERENT position.
  // The REF haplotype uses the variant's original reference position
  // (mLocalRefStart0Idx). Each ALT haplotype records its own position
  // via the per-ALT mLocalHapStart0Idxs map (populated during variant
  // extraction from the SPOA consensus paths).
  //
  //   Haplotype 0 (REF): [ ... var at pos 120, len=3 (REF allele) ... ]
  //   Haplotype 1 (ALT): [ ... var at pos 118, len=5 (ALT allele) ... ]
  //   Haplotype 2 (ALT): [ ... var at pos 119, len=4 (ALT allele) ... ]
  //
  // Returns std::nullopt if the haplotype index doesn't match any allele
  // (the haplotype doesn't carry this variant).
  // ============================================================================
  [[nodiscard]] static auto ExtractHapBounds(RawVariant const& variant, usize aln_hap_idx)
      -> std::optional<HapVariantBounds>;

  // ============================================================================
  // OverlapsAlignment: check if a minimap2 alignment overlaps a variant region.
  //
  // IMPORTANT: ref_start / ref_end in Mm2AlnResult are 0-based coordinates
  // on the haplotype the read was aligned to — NOT genomic reference
  // coordinates. Each haplotype is an independent assembled sequence;
  // coordinates are local to that haplotype string.
  //
  //   Haplotype:      [ ............. full assembled sequence .......... ]
  //   Alignment span: .....[aln.mRefStart ......... aln.mRefEnd).........
  //   Variant A:      [start..end)                              → skip
  //   Variant B:                    [start..end)                → score
  //   Variant C:                                         [start..end) → skip
  //
  // Returns true if the alignment overlaps, i.e. the read covers the
  // variant and can contribute meaningful evidence for allele assignment.
  // ============================================================================
  [[nodiscard]] static auto OverlapsAlignment(Mm2AlnResult const& aln,
                                              HapVariantBounds const& bounds) -> bool;

  static void AddToTable(Result& out_vars_table, cbdg::Read const& qry_read,
                         PerVariantAssignment const& allele_assignments);
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_GENOTYPER_H_
