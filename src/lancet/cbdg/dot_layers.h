#ifndef SRC_LANCET_CBDG_DOT_LAYERS_H_
#define SRC_LANCET_CBDG_DOT_LAYERS_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"

#include <string_view>
#include <utility>

namespace lancet::cbdg {

// ============================================================================
// LogicalEdge — bidirected-faithful edge identity for DOT overlays.
//
// Lancet's de Bruijn graph follows the BCALM2 bidirected model: every
// biological edge between two distinct nodes is stored twice in the
// adjacency lists (`A→B, ++` and the mirror `B→A, −−`), and the renderer
// emits both as separate directed DOT statements. Keying overlays by
// `(SrcId, DstId)` alone collapses two distinct identities — the cause of
// the half-colored walk-overlay bug. Keying by `(SrcId, DstId)` without
// `EdgeKind` also collides legitimate parallel edges between the same
// endpoints (e.g. inverted-repeat hairpins with `++` and `+−` siblings).
//
// LogicalEdge canonicalises both:
//   - endpoints are stored in `(lo, hi)` order so an Edge and its mirror
//     produce the SAME LogicalEdge value;
//   - `EdgeKind` is reversed via `RevEdgeKind` when the swap occurs so
//     `++/−−` mirror pairs collapse but `++/+−` siblings do not.
//
// The renderer expands every overlay LogicalEdge into BOTH `lo->hi` and
// `hi->lo` DOT directions, applying the same style — fixing mirror coverage
// as a renderer invariant rather than a per-call-site concern.
// ============================================================================
struct LogicalEdge {
  // ── 8B Align ────────────────────────────────────────────────────────────
  NodeID mLo = 0;
  NodeID mHi = 0;

  // ── 1B Align ────────────────────────────────────────────────────────────
  EdgeKind mKind = EdgeKind::PLUS_PLUS;

  /// Canonicalise an `Edge`. Mirrors collapse to the same LogicalEdge value;
  /// parallel edges with different `EdgeKind`s remain distinct.
  [[nodiscard]] static auto FromEdge(Edge const& edge) -> LogicalEdge {
    if (edge.SrcId() <= edge.DstId()) {
      return LogicalEdge{.mLo = edge.SrcId(), .mHi = edge.DstId(), .mKind = edge.Kind()};
    }
    return LogicalEdge{.mLo = edge.DstId(), .mHi = edge.SrcId(), .mKind = RevEdgeKind(edge.Kind())};
  }

  template <typename HashState>
  friend auto AbslHashValue(HashState state, LogicalEdge const& edge) -> HashState {
    return HashState::combine(std::move(state), edge.mLo, edge.mHi, static_cast<u64>(edge.mKind));
  }

  friend auto operator==(LogicalEdge const& lhs, LogicalEdge const& rhs) -> bool = default;
};

// ============================================================================
// Style primitives — orthogonal node and edge attribute axes.
//
// Each style field carries an "unset" sentinel (empty string_view, 0 numeric).
// The renderer composes layers in ascending z_order; each non-sentinel field
// in a higher-z layer overrides the same field from lower-z layers, so the
// axes stack independently:
//   • role (z=0)    → mFillColor
//   • anchor (z=10) → mBorderColor + mPenWidth + mPeripheries
//   • probe (z=20)  → mStripeAccent
//
// Edge layers compose `mColor` differently: all layers' colors are
// concatenated into a graphviz colorList (`color="A:B:C"`) so multi-walk
// edges render as parallel-color stripes — no last-writer-wins data loss.
// ============================================================================
struct NodeStyle {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string_view mFillColor;     // empty = no override
  std::string_view mBorderColor;   // empty = no override
  std::string_view mStripeAccent;  // empty = no stripe; non-empty toggles `style="filled,striped"`

  // ── 4B Align ────────────────────────────────────────────────────────────
  u32 mPenWidth = 0;     // 0 = no override
  u32 mPeripheries = 0;  // 0 = no override (1 = default border, 2 = double ring)
};

struct EdgeStyle {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string_view mColor;  // single hex color (#RRGGBB); composed into colorList for shared edges

  // ── 4B Align ────────────────────────────────────────────────────────────
  u32 mPenWidth = 0;  // 0 = no override

  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mDashed = false;
};

// ============================================================================
// Layer — a named set of identifiers + a style + a z-order priority.
//
// The renderer iterates layers in ascending mZOrder; each layer's style
// fields override prior layers' same-field values for nodes/edges in its
// id set (with the colorList composition rule for edge colors). Layer name
// is for debugging and DOT comments only — not part of identity.
// ============================================================================
struct NodeLayer {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string_view mName;
  absl::flat_hash_set<NodeID> mIds;
  NodeStyle mStyle;

  // ── 4B Align ────────────────────────────────────────────────────────────
  i32 mZOrder = 0;
};

struct EdgeLayer {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string_view mName;
  absl::flat_hash_set<LogicalEdge> mEdges;
  EdgeStyle mStyle;

  // ── 4B Align ────────────────────────────────────────────────────────────
  i32 mZOrder = 0;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_LAYERS_H_
