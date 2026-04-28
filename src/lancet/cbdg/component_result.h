#ifndef SRC_LANCET_CBDG_COMPONENT_RESULT_H_
#define SRC_LANCET_CBDG_COMPONENT_RESULT_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph_complexity.h"
#include "lancet/cbdg/path.h"

#include "absl/types/span.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// ComponentResult: per-component assembly output from the de Bruijn graph.
//
// Bundles haplotype paths, the underlying source→sink edge walks, reference
// anchor offset, and graph complexity metrics for one connected component.
// Provides computed accessors for downstream consumers (SPOA MSA, Genotyper,
// VariantAnnotator) without exposing raw path storage. The walks are read
// only by the DOT renderer; caller-layer code consumes Paths via the
// existing accessors.
//
// First haplotype is always the reference. Subsequent ALT haplotypes are
// sorted by descending MinWeight (weakest-link confidence), establishing
// structural priority in downstream SPOA MSA.
// ============================================================================
class ComponentResult {
 public:
  ComponentResult(std::vector<EnumeratedHaplotype> haplotypes, GraphComplexity metrics,
                  u32 anchor_start_offset);

  /// Non-owning views into each path's sequence string. Zero data copy.
  /// Use for SPOA alignment and sequence complexity scoring (both accept string_view).
  [[nodiscard]] auto HaplotypeSequenceViews() const -> std::vector<std::string_view>;

  /// Owning copies of each path's sequence string. Required by minimap2's
  /// mm_idx_str which needs null-terminated c_str() pointers.
  [[nodiscard]] auto HaplotypeSequences() const -> std::vector<std::string>;

  /// Per-haplotype per-base SPOA weights. Each inner vector has one u32 per
  /// sequence base, expanded lazily from Path's run-length-encoded node weights.
  [[nodiscard]] auto HaplotypeWeights() const -> Path::HapWeights;

  /// Max path depth CV across ALT haplotypes (index 1..N). Returns nullopt
  /// when only the reference path exists (no ALT haplotypes).
  [[nodiscard]] auto MaxAltPathCv() const -> std::optional<f64>;

  /// Total path count (REF + ALT). Used for HSE normalization.
  [[nodiscard]] auto NumPaths() const -> usize { return mPaths.size(); }

  /// ALT haplotype count: NumPaths() - 1 (excludes the leading REF path).
  [[nodiscard]] auto NumAltHaplotypes() const -> usize { return mPaths.size() - 1; }

  /// 0-based offset into the window reference where this component's source
  /// anchor node begins. Add to Window::StartPos1() to get genome position.
  [[nodiscard]] auto AnchorStartOffset() const -> u32 { return mAnchorStartOffset; }

  /// Graph topology metrics (cyclomatic complexity, branch points, etc.).
  [[nodiscard]] auto Metrics() const -> GraphComplexity const& { return mMetrics; }

  /// Underlying source→sink edge walks, indexed in lock-step with the
  /// haplotype paths. Walk index 0 (reference) may be empty when REF walk
  /// reconstruction fails. Read only by the DOT renderer.
  [[nodiscard]] auto Walks() const -> absl::Span<std::vector<Edge> const> {
    return absl::MakeConstSpan(mWalks);
  }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<Path> mPaths;
  std::vector<std::vector<Edge>> mWalks;
  GraphComplexity mMetrics;

  // ── 4B Align ────────────────────────────────────────────────────────────
  u32 mAnchorStartOffset = 0;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_COMPONENT_RESULT_H_
