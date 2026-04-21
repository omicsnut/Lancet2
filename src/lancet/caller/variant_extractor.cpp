#include "lancet/caller/variant_extractor.h"

#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_bubble.h"
#include "lancet/core/window.h"

#include "absl/container/btree_set.h"
#include "absl/types/span.h"
#include "spoa/graph.hpp"

#include <algorithm>
#include <limits>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

namespace lancet::caller {

VariantExtractor::VariantExtractor(spoa::Graph const& graph, core::Window const& win,
                                   usize anchor_start)
    : mGraph(graph),
      mWin(win),
      mRefAnchorStart(anchor_start),
      mCurrentRefPos(anchor_start),
      mNumSeqs(mGraph.sequences().size()) {
  // A graph with fewer than 2 sequences has only the REF — no variants to extract.
  if (mNumSeqs < 2) return;

  mCurrentHapPos.assign(mNumSeqs, 0);

  // ===========================================================================
  // RANK LOOKUP INITIALIZATION:  O(N) Inverse Topological Indexing
  // ============================================================================
  // WHY? In SPOA, a node's physical `id` is assigned based on
  // insertion order. A downstream reference base might be Node #5,
  // while an upstream variant inserted later could be Node #500.
  // Using raw `.id` to determine left-to-right ordering produces
  // incorrect results.
  //
  // SPOA provides a topologically sorted node array (rank_to_node).
  // We invert it here: given any node `.id`, `mNodeToRank[id]`
  // returns its true 5'-to-3' rank in O(1).
  // ===========================================================================
  auto const& topological_order = mGraph.rank_to_node();
  // Initialize pointer and rank lookup arrays
  mActivePtrs.assign(mNumSeqs, nullptr);
  mNodeToRank.assign(mGraph.nodes().size(), std::numeric_limits<u32>::max());

  for (u32 rank = 0; rank < topological_order.size(); ++rank) {
    mNodeToRank[topological_order[rank]->id] = rank;
  }

  // Set each haplotype's initial pointer to its first node in the graph
  for (usize i = 0; i < mNumSeqs; ++i) {
    mActivePtrs[i] = mGraph.sequences()[i];
  }
}

void VariantExtractor::SearchAndExtractTo(absl::btree_set<RawVariant>& out_variants) {
  if (mNumSeqs < 2) return;

  while (true) {
    if (AreAllPathsConverged()) {
      // All pointers are nullptr — graph traversal complete
      if (mActivePtrs[REF_HAP_IDX] == nullptr) break;
      AdvanceConvergedPaths();
    } else {
      EatTopologicalBubble(out_variants);
    }
  }
}

// O(M_paths): checks whether all haplotype pointers point to the same node
auto VariantExtractor::AreAllPathsConverged() const -> bool {
  auto const* ref_ptr = mActivePtrs[REF_HAP_IDX];
  return std::ranges::all_of(mActivePtrs | std::views::drop(1),
                             [ref_ptr](auto const* node_ptr) { return node_ptr == ref_ptr; });
}

// Advance all converged pointers to their next successor node
void VariantExtractor::AdvanceConvergedPaths() {
  mPrevMatchNode = mActivePtrs[REF_HAP_IDX];
  for (usize i = 0; i < mNumSeqs; ++i) {
    if (mActivePtrs[i]) {
      mActivePtrs[i] = mActivePtrs[i]->Successor(i);
      mCurrentHapPos[i]++;
    }
  }

  mCurrentRefPos++;
}

// Detect and resolve a single topological bubble, emitting its variants
void VariantExtractor::EatTopologicalBubble(absl::btree_set<RawVariant>& out_variants) {
  // Initialize empty base sequences uniformly
  std::vector<std::string> raw_alleles(mNumSeqs, "");
  std::vector<usize> bubble_hap_starts(mNumSeqs, 0);

  // 1. Anchor
  usize const exact_start_pos =
      InitializeBubbleAnchor(absl::MakeSpan(raw_alleles), bubble_hap_starts);

  // 2. Advance pointers until convergence
  SinkPointers(absl::MakeSpan(raw_alleles));

  // 3. Normalize into a trimmed VCF block
  VariantBubble bubble = CreateNormalizedBubble(exact_start_pos, std::move(raw_alleles));
  bubble.mHapStarts = std::move(bubble_hap_starts);

  if (!bubble.mAltAllelesToHaps.empty()) {
    out_variants.insert(AssembleMultiallelicVariant(std::move(bubble)));
  }
}

// VCF ANCHORING REQUIREMENT:
// Complex indels require a shared prefix match base as anchor.
auto VariantExtractor::InitializeBubbleAnchor(absl::Span<std::string> raw_alleles,
                                              std::vector<usize>& out_hap_starts) -> usize {
  bool const has_prev = (mPrevMatchNode != nullptr);
  auto const anchor_offset = static_cast<usize>(has_prev);

  // Shift genomic coordinate back by 1 if an anchor base exists
  usize const bubble_start_pos = mCurrentRefPos - anchor_offset;

  // Extract the last confirmed match base and prepend it to all alleles
  if (has_prev) {
    char const decoder_val = static_cast<char>(mGraph.decoder(mPrevMatchNode->code));
    std::ranges::for_each(raw_alleles,
                          [decoder_val](std::string& allele) { allele += decoder_val; });
  }

  // Record each haplotype's starting position within the bubble
  for (usize i = 0; i < mNumSeqs; ++i) {
    out_hap_starts[i] = mCurrentHapPos[i] - anchor_offset;
  }

  return bubble_start_pos;
}

// ===========================================================================
// TOPOLOGICAL SINK: advances pointers until convergence
// absl::Span ensures in-place mutation of the caller's raw_alleles vector
// without copying.
// ===========================================================================
// Pushes active pointers forward until all paths converge on the same node.
void VariantExtractor::SinkPointers(absl::Span<std::string> raw_alleles) {
  while (!AreAllPathsConverged()) {
    u32 const min_rank = FindLowestActiveRank();
    if (min_rank == std::numeric_limits<u32>::max()) break;

    ConsumePathsAtRank(min_rank, raw_alleles);
  }
}

// Find the minimum topological rank among all active (non-null) haplotype pointers.
auto VariantExtractor::FindLowestActiveRank() const -> u32 {
  u32 min_rank = std::numeric_limits<u32>::max();
  for (auto const* nptr : mActivePtrs) {
    if (nptr != nullptr) min_rank = std::min(min_rank, mNodeToRank.at(nptr->id));
  }
  return min_rank;
}

// Advance all haplotype pointers that sit at the given topological rank,
// appending their decoded base to the corresponding allele string.
void VariantExtractor::ConsumePathsAtRank(u32 target_rank, absl::Span<std::string> raw_alleles) {
  for (usize i = 0; i < mNumSeqs; ++i) {
    if (mActivePtrs[i] != nullptr && mNodeToRank.at(mActivePtrs[i]->id) == target_rank) {
      raw_alleles[i] += static_cast<char>(mGraph.decoder(mActivePtrs[i]->code));
      mActivePtrs[i] = mActivePtrs[i]->Successor(i);

      // Update per-haplotype position
      mCurrentHapPos[i]++;

      // Only the REF haplotype advances the shared genomic coordinate
      if (i == REF_HAP_IDX) mCurrentRefPos++;
    }
  }
}

// Group per-path sequences upon convergence and apply VCF parsimony trimming.
auto VariantExtractor::CreateNormalizedBubble(usize genome_start_pos,
                                              std::vector<std::string> raw_alleles) const
    -> VariantBubble {
  VariantBubble bubble;
  bubble.mGenomeStartPos = genome_start_pos;
  bubble.mRefAllele = std::move(raw_alleles[REF_HAP_IDX]);

  // Exclude the pure identical reference tracks. Everything else goes into the variant parser
  for (usize alt_idx = 1; alt_idx < mNumSeqs; ++alt_idx) {
    if (raw_alleles[alt_idx] != bubble.mRefAllele) {
      bubble.mAltAllelesToHaps[raw_alleles[alt_idx]].push_back(alt_idx);
    }
  }

  bubble.NormalizeVcfParsimony();
  return bubble;
}

// Classify each ALT allele and assemble the final multiallelic RawVariant.
auto VariantExtractor::AssembleMultiallelicVariant(VariantBubble bubble) -> RawVariant {
  // Build multiallelic VCF record
  RawVariant multi_var;
  multi_var.mChromIndex = mWin.ChromIndex();
  multi_var.mChromName = mWin.ChromName();
  multi_var.mGenomeChromPos1 = bubble.mGenomeStartPos;

  // Set REF coordinate
  multi_var.mLocalRefStart0Idx = bubble.mHapStarts[REF_HAP_IDX];
  multi_var.mRefAllele = std::move(bubble.mRefAllele);

  // Populate disjoint variant structures
  for (auto& [normalized_alt, haps] : bubble.mAltAllelesToHaps) {
    AltAllele sub_alt;
    sub_alt.mSequence = normalized_alt;
    sub_alt.mType = RawVariant::ClassifyVariant(multi_var.mRefAllele, sub_alt.mSequence);
    sub_alt.mLength =
        CalculateVariantLength(multi_var.mRefAllele, sub_alt.mSequence, sub_alt.mType);

    for (usize const hap_id : haps) {
      sub_alt.mLocalHapStart0Idxs.emplace(hap_id, bubble.mHapStarts[hap_id]);
    }

    multi_var.mAlts.push_back(std::move(sub_alt));
  }

  // Sort ALT alleles for deterministic equality comparison and hash stability in btree_set.
  std::ranges::sort(multi_var.mAlts);
  return multi_var;
}

}  // namespace lancet::caller
