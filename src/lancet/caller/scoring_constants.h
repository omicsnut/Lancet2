#ifndef SRC_LANCET_CALLER_SCORING_CONSTANTS_H_
#define SRC_LANCET_CALLER_SCORING_CONSTANTS_H_

#include "lancet/base/types.h"

#include <array>

namespace lancet::caller {

// ── Scoring constants for Illumina read-to-contig realignment ────────────────
//
// These are intentionally STRICT: biological variation is already baked into
// the assembled haplotype sequences, so divergence here is sequencer noise.
// See genotyper.h ALIGNMENT PARADIGM SHIFT comment for full rationale.
inline constexpr int SCORING_MATCH = 1;
inline constexpr int SCORING_MISMATCH = 4;
inline constexpr int SCORING_GAP_OPEN = 12;
inline constexpr int SCORING_GAP_EXTEND = 3;
// ── WFA2 distance penalties (minimap2 score-model translation) ───────────────
//
// Minimap2's KSW2 scoring uses four parameters defined above:
//   a = SCORING_MATCH (1)       — reward per matching base
//   b = SCORING_MISMATCH (4)    — penalty per mismatched base
//   q = SCORING_GAP_OPEN (12)   — penalty to open a new gap
//   e = SCORING_GAP_EXTEND (3)  — penalty per additional gap base
//
// Minimap2 maximizes score: matches earn +a, errors subtract. WFA2 minimizes
// distance: matches cost 0, errors accumulate penalty. To produce identical
// optimal alignments, WFA2 must account for the forfeited match bonus on
// every error. A mismatch in minimap2 costs b plus the lost a (the match
// bonus the base would have earned), giving:
//
//   W_mismatch   = a + b = 1 + 4   = 5
//   W_gap_open   = q     = 12
//   W_gap_extend = a + e = 1 + 3   = 4
//
// Score recovery:  ksw2_score = read_len × a − wfa2_penalty
//   150bp, 0 errors → 150×1 − 0  = 150
//   150bp, 5mm      → 150×1 − 25 = 125  (matches minimap2: 145×1 − 5×4 = 125)
//   150bp, 5bp ins  → 150×1 − 32 = 118  (matches: 145×1 − 12 − 5×3 = 118)
//
// WFA2's coordinate-geometry algorithm strictly requires W_match = 0 to zip
// down the DP diagonal in O(1) time per matching base.
inline constexpr int WFA_MISMATCH = SCORING_MATCH + SCORING_MISMATCH;      // 5
inline constexpr int WFA_GAP_OPEN = SCORING_GAP_OPEN;                       // 12
inline constexpr int WFA_GAP_EXTEND = SCORING_MATCH + SCORING_GAP_EXTEND;   // 4

// ── WFA2 step threshold (fast-path circuit breaker) ──────────────────────────
//
// WFA2 computes alignments in O(s²) time where s = total accumulated penalty.
// MAX_WFA_PENALTY_STEPS caps s, bounding worst-case execution to O(75²) and
// memory to 2×75+1 = 151 integers (~600 bytes, fits entirely in L1 cache).
// WFA2 returns WF_STATUS_MAX_STEPS_REACHED when exceeded → minimap2 fallback.
//
// Biological basis: by Phase 2, variants are baked into assembled haplotypes.
// Read-to-haplotype divergence is sequencer noise, not biology. Illumina error
// rates (0.1−1%) mean >10% divergence indicates a chimera, paralog cross-map,
// or unresolved SV breakpoint — reads that need minimap2's seed-and-chain.
//
// Derivation for 150bp reads at 10% max divergence:
//   15 mismatches × W_mismatch(5) = 75 steps
//   15bp indel    → 12 + 15×4     = 72 steps  (also within budget)
//
// Tuning: S_max = floor(read_len × max_divergence) × W_mismatch.
//   250bp reads → floor(250 × 0.10) × 5 = 125 steps.
//   Increasing S_max trades L1 residency for broader fast-path coverage.
inline constexpr int MAX_WFA_PENALTY_STEPS = 75;

// clang-format off
// ┌───────────────────────────────────────────────┐
// │ 5×5 Scoring Matrix for ComputeLocalScore      │
// │ Target (R) × Query (C) | A=0 C=1 G=2 T=3 N=4 │
// ├───────┬───────┬───────┬───────┬───────┬───────┤
// │       │  A(0) │  C(1) │  G(2) │  T(3) │  N(4) │
// ├───────┼───────┼───────┼───────┼───────┼───────┤
// │  A(0) │    1  │   -4  │   -4  │   -4  │    0  │
// │  C(1) │   -4  │    1  │   -4  │   -4  │    0  │
// │  G(2) │   -4  │   -4  │    1  │   -4  │    0  │
// │  T(3) │   -4  │   -4  │   -4  │    1  │    0  │
// │  N(4) │    0  │    0  │    0  │    0  │    0  │
// └───────┴───────┴───────┴───────┴───────┴───────┘
inline constexpr std::array<i8, 25> SCORING_MATRIX = {
     1, -4, -4, -4,  0,
    -4,  1, -4, -4,  0,
    -4, -4,  1, -4,  0,
    -4, -4, -4,  1,  0,
     0,  0,  0,  0,  0
};

// ┌─────────────────────────────────────────────────────────────┐
// │ ASCII → Numeric Base Encoding Table                         │
// │ A/a → 0, C/c → 1, G/g → 2, T/t → 3, everything else → 4 (N) │
// │ Layout: 256 bytes total (16 rows × 16 hex columns)          │
// └─────────────────────────────────────────────────────────────┘
inline constexpr std::array<u8, 256> ENCODE_TABLE = {
    // 0x0_ (NUL .. SI)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x1_ (DLE .. US)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x2_ (SP .. /)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x3_ (0 .. ?)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x4_ (@ A B C D E F G H I J K L M N O)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x5_ (P Q R S T U V W X Y Z [ \ ] ^ _)
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x6_ (` a b c d e f g h i j k l m n o)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x7_ (p q r s t u v w x y z { | } ~ DEL)
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // 0x8_ - 0xF_ (Extended ASCII blocks / unused)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};
// clang-format on

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_SCORING_CONSTANTS_H_
