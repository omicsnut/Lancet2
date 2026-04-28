#ifndef SRC_LANCET_CBDG_DOT_WALK_LAYERS_H_
#define SRC_LANCET_CBDG_DOT_WALK_LAYERS_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/dot_layers.h"
#include "lancet/cbdg/edge.h"

#include "absl/types/span.h"

#include <string_view>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// WalkColor — fixed accent for REF (index 0) + 64-entry LAB palette for ALTs.
//
// The 64-entry palette is pre-computed by scripts/gen_walk_palette.py via
// k-means clustering in CIE L*a*b* (perceptually uniform) color space, then
// reordered by farthest-first traversal so any prefix [1..K) has maximal
// pairwise separation (critical since most components have only 2-8 walks).
// Indices beyond 64 cycle via modulo through the palette.
// ============================================================================

/// Hex color (`#RRGGBB`) for walk index `walk_index`. Index 0 is the REF
/// accent (white); indices 1..N consume the LAB palette in farthest-first
/// order.
[[nodiscard]] auto WalkColor(usize walk_index) -> std::string_view;

/// Build one EdgeLayer per walk, indexed in lock-step with the input span.
/// Walk 0 (REF) gets the fixed accent color and penwidth=2; ALT walks
/// (1..N) consume the LAB palette starting at index 1. Walks with empty
/// edge sequences produce empty `EdgeLayer.mEdges` and are silently
/// skipped by the renderer (typical for unreconstructable REF backbone).
[[nodiscard]] auto MakeWalkLayers(absl::Span<std::vector<Edge> const> walks)
    -> std::vector<EdgeLayer>;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_WALK_LAYERS_H_
