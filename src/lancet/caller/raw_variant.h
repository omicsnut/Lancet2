#ifndef SRC_LANCET_CALLER_RAW_VARIANT_H_
#define SRC_LANCET_CALLER_RAW_VARIANT_H_

#include "lancet/base/longdust_scorer.h"
#include "lancet/base/sequence_complexity.h"
#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"

#include "absl/strings/str_cat.h"

#include <compare>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::caller {

// ── Graph complexity metrics (coverage-invariant features) ──────────────
// Always populated. 3 fields matching the GRAPH_CX VCF INFO tag.
//
// Coverage stability: GEI uses CovCV (σ/μ, self-normalizing ratio),
// TipToPathCovRatio is a coverage ratio, MaxSingleDirDegree is pure
// topology. All three are coverage-stable above 20×.
//
// Raw topology metrics (CC, BP, EdgeDensity, UnitigRatio, CoverageCv) are
// mathematically correlated and compressed into the Graph Entanglement Index.
// Color-based metrics (UnsharedColorRatio, ColorDiscordantBranches) are not
// topological and are captured by other biologically relevant annotations.
struct GraphMetrics {
  /// GEI: log₁₀(1 + CC×BP×CovCV / UnitigRatio)
  f64 mGraphEntanglementIndex = 0.0;
  /// assembly tearing: tip cov / unitig cov (ratio, self-normalizing)
  f64 mTipToPathCovRatio = 0.0;
  /// hub k-mer detection: max outgoing edges (topology, invariant)
  usize mMaxSingleDirDegree = 0;

  /// Format as 3 comma-separated values for VCF GRAPH_CX INFO tag.
  [[nodiscard]] auto FormatVcfValue() const -> std::string {
    return absl::StrCat(base::FormatComplexityScore(mGraphEntanglementIndex), ",",
                        base::FormatComplexityScore(mTipToPathCovRatio), ",", mMaxSingleDirDegree);
  }
};

class RawVariant {
 public:
  RawVariant() = default;

  // ===========================================================================
  // STRICT SEQUENCE CORE CLASSIFICATION
  // ---------------------------------------------------------------------------
  // Extracts the Mutation Core by squeezing matching 5' prefixes and 3' suffixes.
  // Decouples variant classification from VCF padding constraints.
  // Implementation and full rationale in raw_variant.cpp.
  // ===========================================================================
  [[nodiscard]] static auto ClassifyVariant(std::string_view ref_seq, std::string_view alt_seq)
      -> AlleleType;

  // ── 8B Align ────────────────────────────────────────────────────────────
  usize mChromIndex = SIZE_MAX;

  // GLOBAL GENOMIC COORDINATE: tracks where the variant lands on the
  // reference genome. Used for VCF sorting and emitting.
  usize mGenomeChromPos1 = SIZE_MAX;

  // LOCAL MATRIX COORDINATE (REF): 0-indexed position within the Reference
  // string spanning this Micro-Assembly window. Binds CIGAR strings backwards.
  usize mLocalRefStart0Idx = SIZE_MAX;

  std::string mChromName;
  // Universal left-aligned bounding sequence encompassing ALL ALTs
  std::string mRefAllele;

  // Replaces rigid single scalar ALTs. Vectors avoid heap allocations for
  // small multiallelic blocks (most sites have 1 or 2 alts max).
  std::vector<AltAllele> mAlts;

  // ── Annotation fields (mutable — populated post-construction) ─────────
  // These do not participate in btree_set ordering and are annotated after
  // variant discovery, so they are mutable to allow modification through
  // const btree_set iterators without const_cast.

  mutable GraphMetrics mGraphMetrics;

  // Sequence complexity (11 coverage-invariant features). Distilled from
  // raw multi-scale metrics into 3 groups: Context (REF), Delta (ALT−REF),
  // TR Motif (ALT). Matches the SEQ_CX VCF INFO tag.
  mutable base::SequenceComplexity mSeqCx;

  // Total SPOA haplotype count for this graph component (for HSE
  // normalization). Populated in VariantBuilder::ProcessWindow().
  // Does NOT participate in operator==, operator<, or AbslHashValue.
  mutable usize mNumTotalHaps = 0;

  // Maximum path depth CV across all ALT de Bruijn graph paths.
  // Computed in VariantBuilder::ProcessWindow() from a single O(nhaps)
  // pass over comp_paths[1..]. std::optional because -ffast-math prohibits
  // NaN; nullopt when no ALT path has ≥ 2 nodes.
  mutable std::optional<f64> mMaxPathCv;

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, RawVariant const& var) -> HashState {
    return HashState::combine(std::move(hash_state), var.mChromIndex, var.mGenomeChromPos1,
                              var.mChromName, var.mRefAllele, var.mAlts);
  }

  friend auto operator==(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    // std::vector operator== handles multi-allelic array comparison element-by-element
    return lhs.mChromIndex == rhs.mChromIndex &&
           lhs.mGenomeChromPos1 == rhs.mGenomeChromPos1 &&
           lhs.mChromName == rhs.mChromName &&
           lhs.mRefAllele == rhs.mRefAllele &&
           lhs.mAlts == rhs.mAlts;
  }

  friend auto operator<(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    if (lhs.mChromIndex != rhs.mChromIndex) return lhs.mChromIndex < rhs.mChromIndex;

    if (lhs.mGenomeChromPos1 != rhs.mGenomeChromPos1) {
      return lhs.mGenomeChromPos1 < rhs.mGenomeChromPos1;
    }

    if (lhs.mRefAllele != rhs.mRefAllele) return lhs.mRefAllele < rhs.mRefAllele;
    // Lexical vector comparison behaves safely
    return lhs.mAlts < rhs.mAlts;
  }

  friend auto operator!=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(rhs == lhs);
  }

  friend auto operator>(RawVariant const& lhs, RawVariant const& rhs) -> bool { return rhs < lhs; }

  friend auto operator<=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(rhs < lhs);
  }

  friend auto operator>=(RawVariant const& lhs, RawVariant const& rhs) -> bool {
    return !(lhs < rhs);
  }
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_RAW_VARIANT_H_
