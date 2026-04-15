#include "lancet/caller/variant_set.h"

#include "lancet/caller/variant_extractor.h"

#include "spoa/graph.hpp"

namespace lancet::caller {

// ============================================================================
// Variant Set — Extracted Architecture
//
// ┌──────────────────────────────────────────────────────┐
// │ variant_bubble.h/.cpp                                │
// │   CalculateVariantLength()   — sequence core length  │
// │   VariantBubble              — VCF parsimony trimmer │
// ├──────────────────────────────────────────────────────┤
// │ variant_extractor.h/.cpp                             │
// │   VariantExtractor           — DAG bubble FSM        │
// │   SearchAndExtractTo()       — public entry point    │
// │   AreAllPathsConverged()     — convergence check     │
// │   AdvanceConvergedPaths()    — uniform sweep         │
// │   EatTopologicalBubble()     — bubble resolution     │
// │   SinkPointers()             — greedy sink loop      │
// │   CreateNormalizedBubble()   — allele grouping       │
// │   AssembleMultiallelicVariant() — VCF record build   │
// ├──────────────────────────────────────────────────────┤
// │ variant_set.cpp (this file)                          │
// │   VariantSet(graph, win, start) — constructor        │
// └──────────────────────────────────────────────────────┘
// ============================================================================
VariantSet::VariantSet(spoa::Graph const& graph, core::Window const& win, usize ref_anchor_start) {
  if (graph.sequences().size() < 2) return;

  VariantExtractor extractor(graph, win, ref_anchor_start);
  extractor.SearchAndExtractTo(this->mResultVariants);
}

}  // namespace lancet::caller
