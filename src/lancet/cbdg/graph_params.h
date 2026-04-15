#ifndef SRC_LANCET_CBDG_GRAPH_PARAMS_H_
#define SRC_LANCET_CBDG_GRAPH_PARAMS_H_

#include "lancet/base/types.h"

#include <filesystem>

namespace lancet::cbdg {

/// Minimum kmer length. Must be odd (de Bruijn node identity requires canonical orientation).
/// 13 = smallest odd k where hash collisions are rare on typical 300bp reads.
static constexpr usize DEFAULT_MIN_KMER_LEN = 13;
/// Default maximum kmer length. 127 = practical ceiling where assembly quality plateaus.
static constexpr usize DEFAULT_MAX_KMER_LEN = 127;
/// Hard upper bound on kmer length. 255 = max value for a u8 index into kmer sequences.
static constexpr usize MAX_ALLOWED_KMER_LEN = 255;

/// Minimum total read support across all samples to retain a node.
/// Nodes with coverage < 2 are treated as sequencing errors and pruned.
static constexpr u32 DEFAULT_MIN_NODE_COV = 2;
/// Minimum total read support for a reference k-mer to qualify as a source/sink anchor.
/// Higher than DEFAULT_MIN_NODE_COV to ensure anchors are reliable alignment landmarks.
static constexpr u32 DEFAULT_MIN_ANCHOR_COV = 5;

/// Kmer step size for the outer retry loop. Even step ensures min_k + N*step stays odd
/// (given odd min_k). 6 = empirically good tradeoff between coverage and retry count.
static constexpr u16 DEFAULT_KMER_STEP_LEN = 6;

/// Parameters controlling de Bruijn graph construction and pruning.
struct GraphParams {
  std::filesystem::path mOutGraphsDir;  // 8B+

  usize mMinKmerLen = DEFAULT_MIN_KMER_LEN;  // 8B
  usize mMaxKmerLen = DEFAULT_MAX_KMER_LEN;  // 8B

  u32 mMinNodeCov = DEFAULT_MIN_NODE_COV;      // 4B
  u32 mMinAnchorCov = DEFAULT_MIN_ANCHOR_COV;  // 4B

  u16 mKmerStepLen = DEFAULT_KMER_STEP_LEN;  // 2B
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_PARAMS_H_
