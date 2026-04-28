#ifndef SRC_LANCET_CBDG_DOT_OVERLAY_FACTORIES_H_
#define SRC_LANCET_CBDG_DOT_OVERLAY_FACTORIES_H_

#include "lancet/cbdg/dot_layers.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/probe_tracker.h"

namespace lancet::cbdg {

// ============================================================================
// Layer factories — assemble the standard node-layer stack for the renderer.
//
// Each factory returns a single `NodeLayer` describing one orthogonal axis
// of node styling; the caller composes them into the layer stack passed to
// the renderer. Layer z-orders are conventional:
//   • role layer    z=0   (built per-component by the renderer itself)
//   • anchor layer  z=10  (heavy goldenrod border + peripheries=2)
//   • probe layer   z=20  (striped fill overlay)
//
// All three may apply to the same node simultaneously — they touch
// different attribute slots, so a probe-marked anchor in a CASE-only
// k-mer renders all three signals at once.
// ============================================================================

/// Heavy goldenrod border + double peripheries on the source and sink
/// anchor nodes. The role layer's fillcolor still shows through.
[[nodiscard]] auto MakeAnchorLayer(NodeIDPair const& source_and_sink) -> NodeLayer;

/// Striped-fill overlay on probe-marked nodes for the given component.
/// Stripe accent colour is orchid (matches the legacy `COLOR_HIGHLIGHT`).
/// Returns a layer with empty mIds (no-op) when the tracker has no probe
/// hits in this component.
[[nodiscard]] auto MakeProbeLayer(ProbeTracker const& tracker, usize comp_id,
                                  ProbeTracker::NodeTable const& nodes) -> NodeLayer;

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_OVERLAY_FACTORIES_H_
