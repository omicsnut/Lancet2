#ifndef SRC_LANCET_HTS_MATE_INFO_H_
#define SRC_LANCET_HTS_MATE_INFO_H_

#include "lancet/base/types.h"

namespace lancet::hts {

// ============================================================================
// MateInfo: mate-pair location for out-of-region mate retrieval.
//
// Stores the chromosome index and 0-based start position of a read's mate.
// Used by ReadCollector to fetch mates that fall outside the current assembly
// window but share a query name with an in-region read.
//
// Standalone — no dependency on Alignment or bam1_t. This allows ReadCollector
// to store mate locations without including the full Alignment header.
// ============================================================================
struct MateInfo {
  // ── 8B Alignment ──────────────────────────────────────────────────────
  i64 mMateStartPos0 = -1;  // 0-based mate start position

  // ── 4B Alignment ──────────────────────────────────────────────────────
  i32 mChromIndex = -1;  // reference sequence index (tid) of the mate
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_MATE_INFO_H_
