#ifndef SRC_LANCET_CALLER_SCORING_CONSTANTS_H_
#define SRC_LANCET_CALLER_SCORING_CONSTANTS_H_

#include "lancet/base/types.h"

#include <array>

namespace lancet::caller {

// ============================================================================
// Scoring constants for Illumina read-to-contig realignment
// ============================================================================
//
// These are intentionally STRICT: biological variation is already baked into
// the assembled haplotype sequences, so divergence here is sequencer noise.
// See genotyper.h ALIGNMENT PARADIGM SHIFT comment for full rationale.
inline constexpr int SCORING_MATCH = 1;
inline constexpr int SCORING_MISMATCH = 4;
inline constexpr int SCORING_GAP_OPEN = 12;
inline constexpr int SCORING_GAP_EXTEND = 3;

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
