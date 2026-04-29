#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/path.h"

#include "absl/types/span.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <memory>
#include <string>
#include <string_view>

namespace lancet::caller {

/*
 * ============================================================================
 * SPOA MSA Parameter Rationale for Lancet2 Variant Extraction
 * ============================================================================
 * Values: Match: 0, Mismatch: -6, Gap1: -6,-2, Gap2: -26,-1
 *
 * **SIMD Lane Width Note**: All classical parameters (typically +2/-4) have
 * been shifted downwards by 2.  Setting the Match score to 0 keeps all
 * runtime scores non-positive, which keeps SPOA's WorstCaseAlignmentScore()
 * above the int16 threshold — so the engine selects the faster int16 SIMD
 * path (16 lanes per AVX2 register) instead of falling back to int32
 * (8 lanes, half throughput).  SPOA 4.1.5 dynamically dispatches between
 * int16 and int32 via WorstCaseAlignmentScore().
 * With Match=0 and our gap parameters, the worst-case score for typical
 * alignments (~-2300) is safely above int16 minimum (-32768).
 *
 * Unlike minimap2's `asm5` preset (which aggressively splits contigs at major
 * divergences for whole-genome synteny filtering), these parameters are tuned
 * to force end-to-end global alignment within a specific micro-assembly window
 * to capture dense somatic mutations and large Insertions/Deletions.
 *
 * 1. Why Convex (Dual-Affine) vs. Affine or Linear Scoring:
 *    - Linear Scoring applies a flat penalty per gap base, which is biologically
 *      inaccurate (one 50bp deletion is one biological event, not fifty 1bp
 *      independent events).
 *    - Single Affine Scoring forces a compromise: tune for small variants (strict
 *      extension) and you penalize/clip large insertions/deletions; tune for large
 *      insertions/deletions (loose extension) and sequencer noise creates messy,
 *      spurious small gaps.
 *    - Convex (Dual-Affine) Scoring solves this by taking the minimum of two
 *      intersecting models. It is strict for short gaps to suppress sequencer
 *      noise, but switches to a cheap extension penalty for large
 *      biological gaps.
 *
 * 2. Mismatch Tolerance (Multi-Nucleotide Variants / MNVs):
 *    asm5 uses a +1 match / -19 mismatch, which shatters alignments at dense
 *    mutation clusters. We use 0 / -6. The geometrical difference (6) keeps
 *    the MSA robustly intact while globally forcing alignments through complex variants.
 *
 * 3. Micro-Indel Sensitivity (Convex Model 1: -6, -2):
 *    asm5's -39 gap open penalty prevents small indels, forcing them to misalign
 *    as false-positive SNPs. Our -6 open / -2 extend penalty allows true small
 *    biological indels to open naturally while still applying enough friction
 *    to prevent 1bp sequencing errors (e.g., homopolymer stutters) from opening gaps.
 *
 * 4. Large Insertion/Deletion Continuity (Convex Model 2: -26, -1):
 *    asm5's -81 penalty for large gaps will soft-clip contigs right at an insertion/deletion
 *    breakpoint. Our parameters mathematically intersect at exactly 20bp (6 + 2L = 26 + 1L).
 *    For gaps > 20bp, the algorithm switches to Model 2 where the extension cost
 *    drops to -1. This "cheap extension" forces the DP matrix into mapping massive
 *    insertions/deletions as single, contiguous blocks in the MSA rather than
 *    dropping the alignment entirely.
 *
 *    – Convex dual affine gap scoring -> min(g1+(i-1)*e1, g2+(i-1)*e2)
 */
constexpr i8 MSA_MATCH_SCORE = 0;
constexpr i8 MSA_MISMATCH_SCORE = -6;
constexpr i8 MSA_OPEN1_SCORE = -6;
constexpr i8 MSA_EXTEND1_SCORE = -2;
constexpr i8 MSA_OPEN2_SCORE = -26;
constexpr i8 MSA_EXTEND2_SCORE = -1;

class MsaBuilder {
 public:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::unique_ptr<spoa::AlignmentEngine> mEngine;
  spoa::Graph mGraph;

  MsaBuilder()
      : mEngine(spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, MSA_MATCH_SCORE, MSA_MISMATCH_SCORE, MSA_OPEN1_SCORE,
            MSA_EXTEND1_SCORE, MSA_OPEN2_SCORE, MSA_EXTEND2_SCORE)),
        mGraph(spoa::Graph()) {}

  void UpdateSpoaState(absl::Span<std::string_view const> sequences,
                       absl::Span<cbdg::Path::BaseWeights const> weights);

  /// Render the SPOA alignment graph as a GFA-1.0 document. Caller decides
  /// where the bytes go (typically enqueued for the background flusher).
  [[nodiscard]] auto BuildGfaString() const -> std::string;

  /// Render the multiple sequence alignment as a FASTA document. The
  /// `msa_alns` span comes from `mGraph.GenerateMultipleSequenceAlignment`.
  [[nodiscard]] static auto BuildFastaString(absl::Span<std::string const> msa_alns) -> std::string;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
