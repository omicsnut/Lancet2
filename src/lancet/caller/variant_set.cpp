#include "lancet/caller/variant_set.h"

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "spoa/graph.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cstdint>

namespace {

using lancet::caller::RawVariant;

constexpr usize REF_HAP_IDX = 0;

// =========================================================================================
// STRICT SEQUENCE CORE LENGTH CALCULATOR
// -----------------------------------------------------------------------------------------
// Calculates the biological length of a mutation independent of VCF padding requirements.
// Without extracting the core, MNP lengths are inflated by the shared anchored padding
// bases bounding multi-allelic clusters.
// =========================================================================================
inline auto CalculateVariantLength(std::string_view ref, std::string_view alt,
                                   RawVariant::Type vtype) -> i64 {
  if (vtype == RawVariant::Type::SNV) {
    return 1;
  }

  auto const ref_len = static_cast<i64>(ref.length());
  auto const alt_len = static_cast<i64>(alt.length());
  auto const diff = alt_len - ref_len;

  // Net structural variance equals the string length difference for INS/DEL/CPX
  if (vtype == RawVariant::Type::INS ||
      vtype == RawVariant::Type::DEL ||
      vtype == RawVariant::Type::CPX) {
    return diff;
  }

  // For strict MNPs (diff == 0), the biological length IS exactly the Sequence Core!
  // E.g., REF="ATGC", ALT="ACCC". Trimming `start=1` ('A') and `end=1` ('C') extracts
  // the mutation core `TG`->`CC` (length 4 - 1 - 1 = 2).
  usize start_match = 0;
  while (std::cmp_less(start_match, ref_len) &&
         std::cmp_less(start_match, alt_len) &&
         ref[start_match] == alt[start_match]) {
    start_match++;
  }

  usize end_match = 0;
  while (end_match < (ref_len - start_match) &&
         end_match < (alt_len - start_match) &&
         ref[ref_len - 1 - end_match] == alt[alt_len - 1 - end_match]) {
    end_match++;
  }

  return alt_len - static_cast<i64>(start_match) - static_cast<i64>(end_match);
}

// =========================================================================================
// TOPOLOGICAL BUBBLE EXTRACTION — How Variants Are Discovered in the POA Graph
// =========================================================================================
// Rather than relying on 2D string matrices (MSA column diffing)—which scale poorly—or
// pairwise mapping—which fractures overlapping multiallelic loci—this algorithm tracks
// pointers traversing the SPOA directed acyclic graph (DAG).
//
// -> 1. BIALLELIC SNVs (Paths diverge at a single node, then reconverge):
//
//                                [REF]
//                          .--> (T)[3] --.
//                         /               \
//   Anchor: (A)[2] ------+                 +-----> Target: (G)[5] (CONVERGED!)
//                         \               /
//                          `--> (C)[4] --'
//                                [ALT]
//
// -> 2. DELETIONS (ALT path skips a stretch of REF nodes via a direct edge):
//
//                                            [REF]
//                          .--> (T)[3] --> (C)[4] --> (G)[5] --.
//                         /                                     \
//   Anchor: (A)[2] ------+                                       +-----> Target: (T)[6]
//   (CONVERGED!)
//                         \                                     /
//                          `-----------------------------------'
//                                            [ALT]
//
// -> 3. MULTIALLELIC COMPLEXES (Multiple ALT paths diverge simultaneously):
//       Resolving these requires advancing all path pointers in topological order.
//
//                               [ALT 1]
//                          .--> (C)[3] --> (A)[4] ---------.
//                         /                                 \
//   Anchor: (T)[2] ------+----> (T)[5] --> (G)[6] --> (C)[7] +-----> Target: (A)[10] (CONVERGED!)
//                         \     [REF]                       /
//                          `--> (G)[8] --> (A)[9] ---------'
//                               [ALT 2]
//
// HOW THE SWEEP EXTRACTOR ALGORITHM WORKS ("Greedy Sink"):
// 1. Maintain an array of active pointers — one per haplotype — tracking the
//    current node each path occupies in the DAG.
// 2. When pointers disagree (point to different nodes), a "Bubble" has opened.
// 3. Each pointer has a `Rank` — its topologically sorted 5'-to-3' position
//    in the graph. Rank is the linearized left-to-right ordering.
// 4. Repeatedly advance only the pointer(s) with the LOWEST rank, appending
//    their decoded bases to per-path string buffers.
// 5. This causes lagging paths to catch up until ALL pointers converge on
//    the same Rank — the bubble's convergence node.
// 6. Collect the buffered sequences, pass them through VCF multi-allelic
//    trimming (NormalizeVcfParsimony), and emit the RawVariant.
// ================================================================================================
// Encapsulates VCF string trimming and left-alignment, separate from graph traversal logic.
class VariantBubble {
 public:
  // Maps each distinct ALT allele string to the haplotype IDs that produced it.
  // O(1) lookup per allele via Abseil flat hash map.
  absl::flat_hash_map<std::string, std::vector<usize>> mAltAllelesToHaps;
  std::string mRefAllele;            // 24B
  std::vector<usize> mHapStarts;     // 24B
  usize mGenomeStartPos = SIZE_MAX;  // 8B

  // Multi-allelic normalization trims shared prefix/suffix across ALL alleles simultaneously.
  // Trimming must stop if any ALT would become empty (losing anchor integrity for INDELs).
  void NormalizeVcfParsimony() {
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
    ApplyUnifiedTrim(CAN_MATCH_RIGHT, DO_TRIM_RIGHT);

    usize const initial_ref_len = mRefAllele.length();
    // Left Trim (e.g. REF: "TTC", ALTS: ["TGC"] => "TC", ["GC"])
    ApplyUnifiedTrim(CAN_MATCH_LEFT, DO_TRIM_LEFT);

    // After left-trimming, advance the genomic start position by the number of characters removed.
    mGenomeStartPos += (initial_ref_len - mRefAllele.length());
  }

 private:
  // Generic trim loop: `can_trim` checks whether the boundary character matches across all
  // alleles (using string_view to avoid copies), `do_trim` performs the string mutation.
  template <typename CanTrimFunc, typename DoTrimFunc>
  void ApplyUnifiedTrim(CanTrimFunc can_trim, DoTrimFunc do_trim) {
    // Guard: REF must retain at least 1 base (VCF requires non-empty REF)
    while (mRefAllele.length() > 1) {
      std::string_view const ref_view(mRefAllele);

      // If any allele's boundary character differs, trimming stops for all alleles
      bool const is_trimmable = std::ranges::all_of(mAltAllelesToHaps, [&](auto const& pair) {
        return can_trim(ref_view, std::string_view(pair.first));
      });

      if (!is_trimmable) break;

      // All alleles match — apply the trim
      do_trim(mRefAllele);

      // Flat hash map keys are const (mutating them would corrupt the hash index).
      // Rebuild the map with trimmed keys, using std::move on the haplotype vectors
      // to avoid copying.
      absl::flat_hash_map<std::string, std::vector<usize>> rehashed_map;
      rehashed_map.reserve(mAltAllelesToHaps.size());

      for (auto& [alt_seq, haps] : mAltAllelesToHaps) {
        std::string new_alt = alt_seq;
        do_trim(new_alt);
        rehashed_map.try_emplace(std::move(new_alt), std::move(haps));
      }

      mAltAllelesToHaps = std::move(rehashed_map);
    }
  }
};

// =========================================================================================
// VariantExtractor — Finite State Machine for DAG Bubble Detection
// -----------------------------------------------------------------------------------------
// WHY A CLASS? The extraction loop maintains 6+ mutable state variables
// (active_ptrs, node_to_rank, current_ref_pos, etc.). Keeping them as class
// members avoids a single monolithic function and lets each step
// (InitializeBubbleAnchor, SinkPointers, CreateNormalizedBubble) read
// as self-contained operations.
// =========================================================================================
class VariantExtractor {
 public:
  VariantExtractor(spoa::Graph const& graph, lancet::core::Window const& win, usize anchor_start)
      : mGraph(graph),
        mWin(win),
        mRefAnchorStart(anchor_start),
        mCurrentRefPos(anchor_start),
        mNumSeqs(mGraph.sequences().size()) {
    // A graph with fewer than 2 sequences has only the REF — no variants to extract.
    if (mNumSeqs < 2) return;

    mCurrentHapPos.assign(mNumSeqs, 0);  // Initializes array accurately bounding paths

    // ===========================================================================
    // RANK LOOKUP INITIALIZATION:  O(N) Inverse Topological Indexing
    // ---------------------------------------------------------------------------
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

  // ===================================================================================================
  // VARIANT EXTRACTOR FLOWCHART
  // ===================================================================================================
  // Pushes active pointers 5'-to-3' through the POA graph. When pathways diverge
  // (a variant occurs), the extractor halts uniform sweeping and triggers a Bubble.
  //
  //                              (C)  (Prior Match Node — all pointers converged here)
  //                            /     \
  //                  (ALT)    /       \   (REF)
  //                 (C)->(T) .         . (T)->(G)
  //                           \       /
  //                            \     /
  //                              (G)      (Target Convergence Node — paths reunite here)
  //
  // Private helper call sequence:
  // 1. `InitializeBubbleAnchor`    : Prepends the last matched base (C) as the VCF anchor.
  // 2. `SinkPointers`              : The hot loop. Finds the lowest-rank pointer,
  //    |-- `FindLowestActiveRank`  :   advances it, and repeats until all pointers
  //    |-- `ConsumePathsAtRank`    :   converge on (G).
  // 3. `CreateNormalizedBubble`    : Groups per-path strings, runs VCF parsimony trimming.
  // 4. `AssembleMultiallelicVariant`: Classifies each ALT and emits a RawVariant.
  // ===================================================================================================
  // Single public interface endpoint parsing variants continuously
  // directly into the tracking payload container
  void SearchAndExtractTo(absl::btree_set<RawVariant>& out_variants) {
    if (mNumSeqs < 2) {
      return;
    }

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

 private:
  // Memory Alignment: 24B vector objects first, then 8B references/integers, then pointers.
  std::vector<u32> mNodeToRank;
  std::vector<spoa::Graph::Node const*> mActivePtrs;
  // Current position index within each haplotype's sequence
  std::vector<usize> mCurrentHapPos;
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  spoa::Graph const& mGraph;
  lancet::core::Window const& mWin;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
  usize mRefAnchorStart;
  usize mCurrentRefPos;
  usize mNumSeqs;
  // Last converged VCF anchor node (prepended to bubble sequences)
  spoa::Graph::Node const* mPrevMatchNode = nullptr;

  // O(M_paths): checks whether all haplotype pointers point to the same node
  [[nodiscard]] auto AreAllPathsConverged() const -> bool {
    auto const* target = mActivePtrs[REF_HAP_IDX];
    for (usize i = 1; i < mNumSeqs; ++i) {
      if (mActivePtrs[i] != target) return false;
    }
    return true;
  }

  // Advance all converged pointers to their next successor node
  void AdvanceConvergedPaths() {
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
  void EatTopologicalBubble(absl::btree_set<RawVariant>& out_variants) {
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
  auto InitializeBubbleAnchor(absl::Span<std::string> raw_alleles,
                              std::vector<usize>& out_hap_starts) -> usize {
    bool const has_prev = (mPrevMatchNode != nullptr);
    auto const anchor_offset = static_cast<usize>(has_prev);

    // Shift genomic coordinate back by 1 if an anchor base exists
    usize const bubble_start_pos = mCurrentRefPos - anchor_offset;

    // Extract the last confirmed match base and prepend it to all alleles
    if (has_prev) {
      char const decoder_val = static_cast<char>(mGraph.decoder(mPrevMatchNode->code));
      std::for_each(raw_alleles.begin(), raw_alleles.end(),
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
  void SinkPointers(absl::Span<std::string> raw_alleles) {
    while (!AreAllPathsConverged()) {
      u32 const min_rank = FindLowestActiveRank();
      if (min_rank == std::numeric_limits<u32>::max()) break;

      ConsumePathsAtRank(min_rank, raw_alleles);
    }
  }

  // Phase 1: Determine the lowest active topological boundary amongst traversing sets
  [[nodiscard]] auto FindLowestActiveRank() const -> u32 {
    u32 min_rank = std::numeric_limits<u32>::max();
    for (auto const* nptr : mActivePtrs) {
      if (nptr != nullptr) min_rank = std::min(min_rank, mNodeToRank.at(nptr->id));
    }
    return min_rank;
  }

  // Phase 2: Selectively consume and roll exclusively paths
  // pinned precisely against that lowest rank
  void ConsumePathsAtRank(u32 target_rank, absl::Span<std::string> raw_alleles) {
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
  [[nodiscard]] auto CreateNormalizedBubble(usize genome_start_pos,
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
  auto AssembleMultiallelicVariant(VariantBubble bubble) -> RawVariant {
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
      RawVariant::AltAllele sub_alt;
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
    std::sort(multi_var.mAlts.begin(), multi_var.mAlts.end());
    return multi_var;
  }
};

}  // namespace

namespace lancet::caller {

void VariantSet::ExtractVariantsFromGraph(spoa::Graph const& graph, core::Window const& win,
                                          usize ref_anchor_start) {
  if (graph.sequences().size() < 2) return;

  VariantExtractor extractor(graph, win, ref_anchor_start);
  extractor.SearchAndExtractTo(this->mResultVariants);
}

}  // namespace lancet::caller
