#include "lancet/cbdg/graph_complexity.h"

#include "lancet/base/compute_stats.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"

#include <algorithm>

namespace lancet::cbdg {

auto ComputeGraphComplexity(Graph const& graph, usize const component_id) -> GraphComplexity {
  GraphComplexity cplx;

  // Local accumulators — V and E are intermediates for CC, not stored on the class.
  usize num_nodes = 0;
  usize num_edges = 0;
  usize unitig_nodes = 0;

  // Streaming statistics — avoid heap-allocated coverage vectors by computing
  // mean and standard deviation in a single pass via Welford's online algorithm.
  // Separate accumulators per node category (all / tip / unitig) for the
  // tip-to-path coverage ratio.
  lancet::base::OnlineStats cov_stats;
  lancet::base::OnlineStats tip_stats;
  lancet::base::OnlineStats unitig_stats;

  for (auto const& [nid, node_ptr] : graph.Nodes()) {
    if (node_ptr->GetComponentId() != component_id) continue;

    num_nodes++;
    auto const dflt_sign = node_ptr->SignFor(Kmer::Ordering::DEFAULT);
    usize dflt_dir_edges = 0;
    usize oppo_dir_edges = 0;

    for (Edge const& edge : *node_ptr) {
      if (edge.SrcSign() == dflt_sign) {
        dflt_dir_edges++;
      } else {
        oppo_dir_edges++;
      }
    }

    // Total edges (including mirrors stored at both endpoints); halved below
    num_edges += dflt_dir_edges + oppo_dir_edges;
    auto const max_dir = std::max(dflt_dir_edges, oppo_dir_edges);
    cplx.mMaxSingleDirDegree = std::max(cplx.mMaxSingleDirDegree, max_dir);

    bool const is_branch = (dflt_dir_edges >= 2 || oppo_dir_edges >= 2);
    if (is_branch) cplx.mNumBranchPoints++;

    // Unitig ratio: nodes with exactly 1-in and 1-out (linear chain)
    if (dflt_dir_edges == 1 && oppo_dir_edges == 1) unitig_nodes++;

    // Feed coverage into streaming accumulators
    auto const cov = static_cast<f64>(node_ptr->TotalReadSupport());
    cov_stats.Add(cov);

    // Categorize for tip-to-path ratio
    if (dflt_dir_edges == 0 || oppo_dir_edges == 0) {
      tip_stats.Add(cov);
    } else if (dflt_dir_edges == 1 && oppo_dir_edges == 1) {
      unitig_stats.Add(cov);
    }
  }

  // Each edge stored at both endpoints (forward + mirror) → divide by 2
  num_edges /= 2;
  cplx.mCyclomaticComplexity = (num_edges >= num_nodes) ? (num_edges - num_nodes + 1) : 0;

  // Unitig ratio
  cplx.mUnitigRatio =
      num_nodes > 0 ? static_cast<f64>(unitig_nodes) / static_cast<f64>(num_nodes) : 0.0;

  // Coverage CV (σ/μ) — Welford's online algorithm avoids the numerical
  // instability of the naive two-pass formula (Var = E[x²] − E[x]²), which
  // gives wrong answers when the mean is large relative to the spread
  // (e.g., coverage values around 30× with σ ≈ 2).
  if (!cov_stats.IsEmpty() && cov_stats.Mean() > 0.0) {
    cplx.mCoverageCv = cov_stats.StdDev() / cov_stats.Mean();
  }

  // Tip-to-path coverage ratio
  if (!tip_stats.IsEmpty() && !unitig_stats.IsEmpty() && unitig_stats.Mean() > 0.0) {
    cplx.mTipToPathCovRatio = tip_stats.Mean() / unitig_stats.Mean();
  }

  return cplx;
}

}  // namespace lancet::cbdg
