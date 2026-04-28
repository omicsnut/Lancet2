#ifndef SRC_LANCET_CBDG_DOT_PLAN_H_
#define SRC_LANCET_CBDG_DOT_PLAN_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/dot_layers.h"
#include "lancet/cbdg/probe_tracker.h"

#include <string>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// DotSnapshotKind — orthogonal axis to PruneStage.
//
// PruneStage describes which pruning boundary a snapshot was captured at;
// DotSnapshotKind describes whether the snapshot is an intermediate stage
// (carries a PruneStage tag) or the final post-prune snapshot (carries an
// optional walk overlay). Splitting these axes is what `GraphState`
// previously failed to do — the old enum jammed both axes into one.
// ============================================================================
enum class DotSnapshotKind : u8 {
  PRUNE_STAGE,  ///< intermediate stage; mPruneStage names the pruning boundary
  FINAL,        ///< post-prune; walk overlay present iff mEdgeLayers is non-empty
};

// ============================================================================
// GraphSnapshotMode — runtime control of DOT snapshot verbosity.
//
// Plumbed from the `--graph-snapshots` CLI flag through GraphParams. FINAL
// (default) emits one DOT per component per window. VERBOSE additionally
// emits PRUNE_STAGE snapshots after each pruning boundary that survives the
// pre-compression cull (compress1, lowcov2, compress2, tips). The verbose
// mode replaces the old `LANCET_DEVELOP_MODE` compile-time gate; no rebuild
// required.
// ============================================================================
enum class GraphSnapshotMode : u8 {
  FINAL = 0,
  VERBOSE = 1,
};

// ============================================================================
// DotPlan — fully-resolved render input.
//
// Bundles everything the renderer needs to emit a single DOT file: the
// component to render, a stable subgraph name (used as the DOT `subgraph`
// identifier and matched to the eventual filename stem), the snapshot kind
// + optional PruneStage tag, and the ordered layer stacks. Construction
// happens in graph.cpp (BufferFinalSnapshot / BufferStageSnapshot); the
// renderer is policy-free and consumes the assembled plan.
// ============================================================================
struct DotPlan {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<NodeLayer> mNodeLayers;
  std::vector<EdgeLayer> mEdgeLayers;
  std::string mSubgraphName;
  usize mCompId = 0;

  // ── 1B Align ────────────────────────────────────────────────────────────
  DotSnapshotKind mKind = DotSnapshotKind::FINAL;
  PruneStage mPruneStage =
      PruneStage::PRUNED_AT_BUILD;  ///< meaningful only when mKind == PRUNE_STAGE
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_PLAN_H_
