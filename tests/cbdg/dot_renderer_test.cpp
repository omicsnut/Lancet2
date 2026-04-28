#include "lancet/cbdg/dot_renderer.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/dot_layers.h"
#include "lancet/cbdg/dot_overlay_factories.h"
#include "lancet/cbdg/dot_plan.h"
#include "lancet/cbdg/dot_snapshot_buffer.h"
#include "lancet/cbdg/dot_walk_layers.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/node.h"

#include "absl/container/flat_hash_map.h"
#include "catch_amalgamated.hpp"

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using lancet::cbdg::DotPlan;
using lancet::cbdg::DotSnapshotBuffer;
using lancet::cbdg::DotSnapshotKind;
using lancet::cbdg::Edge;
using lancet::cbdg::EdgeKind;
using lancet::cbdg::Kmer;
using lancet::cbdg::Label;
using lancet::cbdg::LogicalEdge;
using lancet::cbdg::MakeAnchorLayer;
using lancet::cbdg::MakeFwdEdgeKind;
using lancet::cbdg::MakeWalkLayers;
using lancet::cbdg::Node;
using lancet::cbdg::NodeID;
using lancet::cbdg::NodeIDPair;
using lancet::cbdg::ReconstructRefWalk;
using lancet::cbdg::RevEdgeKind;
using lancet::cbdg::SerializeToDotString;
using lancet::cbdg::WalkColor;

namespace {

// ── Test scaffold: hand-built NodeTable for renderer unit tests ───────────
// Mirrors the helper in graph_test.cpp; copied locally to keep test files
// independent.
struct TestGraph {
  absl::flat_hash_map<NodeID, std::unique_ptr<Node>> mNodes;
  std::vector<NodeID> mNodeIds;

  auto AddNode(std::string_view seq, Label label = Label(Label::REFERENCE)) -> NodeID {
    auto mer = Kmer(seq);
    auto const node_id = mer.Identifier();
    mNodes.try_emplace(node_id, std::make_unique<Node>(std::move(mer), label));
    mNodeIds.push_back(node_id);
    return node_id;
  }

  // Add a forward edge from src to dst using the canonical signs of both nodes.
  // Stores both the forward edge and its mirror, matching how Graph::AddNodes
  // populates adjacency lists.
  auto AddEdge(NodeID src, NodeID dst) -> Edge {
    auto& src_node = mNodes.at(src);
    auto& dst_node = mNodes.at(dst);
    auto const src_sign = src_node->SignFor(Kmer::Ordering::DEFAULT);
    auto const dst_sign = dst_node->SignFor(Kmer::Ordering::DEFAULT);
    auto const kind = MakeFwdEdgeKind({src_sign, dst_sign});
    Edge const fwd(NodeIDPair{src, dst}, kind);
    src_node->EmplaceEdge(NodeIDPair{src, dst}, kind);
    dst_node->EmplaceEdge(NodeIDPair{dst, src}, RevEdgeKind(kind));
    return fwd;
  }

  // Like AddEdge but with an explicit EdgeKind (for hairpin / EdgeKind tests).
  auto AddEdgeKind(NodeID src, NodeID dst, EdgeKind kind) -> Edge {
    Edge const fwd(NodeIDPair{src, dst}, kind);
    mNodes.at(src)->EmplaceEdge(NodeIDPair{src, dst}, kind);
    mNodes.at(dst)->EmplaceEdge(NodeIDPair{dst, src}, RevEdgeKind(kind));
    return fwd;
  }

  void SetAllComponentId(usize const comp_id) {
    for (auto& [node_id, node_ptr] : mNodes) {
      node_ptr->SetComponentId(comp_id);
    }
  }
};

// Count occurrences of `needle` in `haystack`.
auto CountOccurrences(std::string_view haystack, std::string_view needle) -> usize {
  if (needle.empty()) return 0;
  usize count = 0;
  usize pos = 0;
  while ((pos = haystack.find(needle, pos)) != std::string_view::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

}  // namespace

// ============================================================================
//  LogicalEdge canonicalization
// ============================================================================

TEST_CASE("LogicalEdge collapses anti-parallel mirrors", "[lancet][cbdg][LogicalEdge]") {
  // Build two NodeIDs that are guaranteed distinct so the canonical (lo, hi)
  // ordering swap actually fires.
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  REQUIRE(nid_a != nid_b);
  auto const lo_id = std::min(nid_a, nid_b);
  auto const hi_id = std::max(nid_a, nid_b);

  Edge const fwd(NodeIDPair{lo_id, hi_id}, EdgeKind::PLUS_PLUS);
  Edge const mirror = fwd.MirrorEdge();
  REQUIRE(mirror.Kind() == EdgeKind::MINUS_MINUS);

  auto const le_fwd = LogicalEdge::FromEdge(fwd);
  auto const le_mirror = LogicalEdge::FromEdge(mirror);

  CHECK(le_fwd == le_mirror);
  CHECK(le_fwd.mLo == lo_id);
  CHECK(le_fwd.mHi == hi_id);
  CHECK(le_fwd.mKind == EdgeKind::PLUS_PLUS);
}

TEST_CASE("LogicalEdge keeps EdgeKind siblings distinct", "[lancet][cbdg][LogicalEdge]") {
  // Two edges A→B with different EdgeKinds (the hairpin bug cause). They
  // share endpoints but must remain distinct LogicalEdge values so a walk
  // overlay over only one EdgeKind doesn't recolor the other.
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");

  Edge const edge_pp(NodeIDPair{nid_a, nid_b}, EdgeKind::PLUS_PLUS);
  Edge const edge_pm(NodeIDPair{nid_a, nid_b}, EdgeKind::PLUS_MINUS);

  auto const le_pp = LogicalEdge::FromEdge(edge_pp);
  auto const le_pm = LogicalEdge::FromEdge(edge_pm);

  CHECK(le_pp != le_pm);
  CHECK(le_pp.mLo == le_pm.mLo);
  CHECK(le_pp.mHi == le_pm.mHi);
  CHECK(le_pp.mKind != le_pm.mKind);
}

// ============================================================================
//  WalkColor palette
// ============================================================================

TEST_CASE("WalkColor index 0 is the REF accent", "[lancet][cbdg][WalkColor]") {
  CHECK(WalkColor(0) == "#FFFFFF");
}

TEST_CASE("WalkColor ALT walks consume the LAB palette", "[lancet][cbdg][WalkColor]") {
  // Adjacent ALT colors must differ (farthest-first ordering of the palette).
  // This isn't an exhaustive Delta-E check; just a sanity guarantee that
  // sequential ALTs aren't accidentally identical.
  for (usize idx = 1; idx < 10; ++idx) {
    CHECK(WalkColor(idx) != WalkColor(idx + 1));
  }
}

TEST_CASE("WalkColor wraps via modulo past 64 ALT entries", "[lancet][cbdg][WalkColor]") {
  // 64-entry palette + 1 REF index = walk index 65 should wrap to ALT[0]
  // (which is walk index 1's color).
  CHECK(WalkColor(1) == WalkColor(1 + 64));
}

// ============================================================================
//  MakeWalkLayers
// ============================================================================

TEST_CASE("MakeWalkLayers builds one layer per walk in lock-step", "[lancet][cbdg][WalkLayers]") {
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  auto const edge_ab = tgraph.AddEdge(nid_a, nid_b);
  auto const edge_bc = tgraph.AddEdge(nid_b, nid_c);

  std::vector<std::vector<Edge>> const walks{
      {edge_ab, edge_bc},  // REF walk
      {edge_ab},           // ALT walk
  };

  auto const layers = MakeWalkLayers(walks);
  REQUIRE(layers.size() == 2);

  // REF layer: walk index 0, fixed accent, two LogicalEdges.
  CHECK(layers[0].mStyle.mColor == WalkColor(0));
  CHECK(layers[0].mEdges.size() == 2);
  CHECK(layers[0].mEdges.contains(LogicalEdge::FromEdge(edge_ab)));
  CHECK(layers[0].mEdges.contains(LogicalEdge::FromEdge(edge_bc)));

  // ALT layer: walk index 1, palette[0] color, one LogicalEdge.
  CHECK(layers[1].mStyle.mColor == WalkColor(1));
  CHECK(layers[1].mEdges.size() == 1);
  CHECK(layers[1].mEdges.contains(LogicalEdge::FromEdge(edge_ab)));

  // ALT walk must layer above REF so its color appends later in the
  // colorList, keeping run-to-run rendering deterministic.
  CHECK(layers[1].mZOrder > layers[0].mZOrder);
}

// ============================================================================
//  MakeAnchorLayer
// ============================================================================

TEST_CASE("MakeAnchorLayer carries border + double peripheries", "[lancet][cbdg][AnchorLayer]") {
  NodeIDPair const anchors{42, 99};
  auto const layer = MakeAnchorLayer(anchors);

  CHECK(layer.mIds.contains(42));
  CHECK(layer.mIds.contains(99));
  CHECK(layer.mIds.size() == 2);
  CHECK_FALSE(layer.mStyle.mBorderColor.empty());
  CHECK(layer.mStyle.mPenWidth > 0);
  CHECK(layer.mStyle.mPeripheries == 2);
  // Anchor layer sits between role (z=0) and probe (z=20).
  CHECK(layer.mZOrder > 0);
  CHECK(layer.mZOrder < 20);
}

// ============================================================================
//  ReconstructRefWalk
// ============================================================================

TEST_CASE("ReconstructRefWalk traces surviving backbone", "[lancet][cbdg][RefWalk]") {
  // Build a 3-node REF chain A→B→C in component 1. mRefNodeIds carries
  // [A, B, C]; all three survive, so the reconstructed walk should be
  // [edge(A,B), edge(B,C)].
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  auto const edge_ab = tgraph.AddEdge(nid_a, nid_b);
  auto const edge_bc = tgraph.AddEdge(nid_b, nid_c);
  tgraph.SetAllComponentId(1);

  std::vector<NodeID> const ref_ids{nid_a, nid_b, nid_c};
  auto const walk = ReconstructRefWalk(absl::MakeConstSpan(ref_ids), tgraph.mNodes, 1);

  REQUIRE(walk.size() == 2);
  CHECK(walk[0] == edge_ab);
  CHECK(walk[1] == edge_bc);
}

TEST_CASE("ReconstructRefWalk skips absorbed k-mers", "[lancet][cbdg][RefWalk]") {
  // Simulate compression: B was absorbed into A so it no longer exists in
  // the node table. mRefNodeIds still mentions B, but the surviving REF
  // walk should be [A → C] (with the direct edge created during compression).
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b_phantom = Kmer("CGTACGTACGA").Identifier();  // B never added
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  auto const edge_ac = tgraph.AddEdge(nid_a, nid_c);
  tgraph.SetAllComponentId(1);

  std::vector<NodeID> const ref_ids{nid_a, nid_b_phantom, nid_c};
  auto const walk = ReconstructRefWalk(absl::MakeConstSpan(ref_ids), tgraph.mNodes, 1);

  REQUIRE(walk.size() == 1);
  CHECK(walk[0] == edge_ac);
}

TEST_CASE("ReconstructRefWalk returns empty on fragmented backbone", "[lancet][cbdg][RefWalk]") {
  // Two surviving REF nodes with no connecting edge → backbone fragmented.
  // The function must abandon and return empty rather than half a walk.
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  // No AddEdge between A and C
  tgraph.SetAllComponentId(1);

  std::vector<NodeID> const ref_ids{nid_a, nid_c};
  auto const walk = ReconstructRefWalk(absl::MakeConstSpan(ref_ids), tgraph.mNodes, 1);

  CHECK(walk.empty());
}

// ============================================================================
//  Renderer — bidirected mirror coverage
// ============================================================================

TEST_CASE("Renderer overlays both anti-parallel splines for a walked edge",
          "[lancet][cbdg][Renderer]") {
  // The headline bug fix: a walk edge A→B must color BOTH the A->B and the
  // B->A DOT statements (graphviz strict-mode merges attributes per direction;
  // neato draws each as a separate spline).
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const edge_ab = tgraph.AddEdge(nid_a, nid_b);
  tgraph.SetAllComponentId(7);

  std::vector<std::vector<Edge>> const walks{{edge_ab}};
  DotPlan plan;
  plan.mKind = DotSnapshotKind::FINAL;
  plan.mCompId = 7;
  plan.mSubgraphName = "test_mirror";
  plan.mEdgeLayers = MakeWalkLayers(walks);

  auto const dot = SerializeToDotString(tgraph.mNodes, plan);

  // Both directions of the walked edge must appear with the WalkColor(0) value.
  // The renderer emits an overlay statement of the form
  //   `<src> -> <dst> [color="#FFFFFF" ...]`.
  std::string const fwd_overlay = std::to_string(nid_a) +
                                  " -> " +
                                  std::to_string(nid_b) +
                                  R"raw( [color=")raw" +
                                  std::string(WalkColor(0));
  std::string const rev_overlay = std::to_string(nid_b) +
                                  " -> " +
                                  std::to_string(nid_a) +
                                  R"raw( [color=")raw" +
                                  std::string(WalkColor(0));

  CHECK(dot.find(fwd_overlay) != std::string::npos);
  CHECK(dot.find(rev_overlay) != std::string::npos);
}

// ============================================================================
//  Renderer — multi-walk colorList composition
// ============================================================================

TEST_CASE("Renderer composes shared edges into a colorList stripe", "[lancet][cbdg][Renderer]") {
  // Three walk layers all covering the same edge → emitted as
  // `color="A:B:C"` with all three colors concatenated. No last-writer-wins.
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const edge_ab = tgraph.AddEdge(nid_a, nid_b);
  tgraph.SetAllComponentId(7);

  std::vector<std::vector<Edge>> const walks{{edge_ab}, {edge_ab}, {edge_ab}};
  DotPlan plan;
  plan.mKind = DotSnapshotKind::FINAL;
  plan.mCompId = 7;
  plan.mSubgraphName = "test_stripe";
  plan.mEdgeLayers = MakeWalkLayers(walks);

  auto const dot = SerializeToDotString(tgraph.mNodes, plan);

  // The colorList contains 3 entries separated by colons.
  std::string const expected_color_list =
      std::string(WalkColor(0)) + ":" + std::string(WalkColor(1)) + ":" + std::string(WalkColor(2));
  CHECK(dot.find(expected_color_list) != std::string::npos);
}

// ============================================================================
//  Renderer — empty walks produce no overlay block
// ============================================================================

TEST_CASE("Renderer emits no overlay block when walks vector is empty",
          "[lancet][cbdg][Renderer]") {
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  tgraph.AddEdge(nid_a, nid_b);
  tgraph.SetAllComponentId(7);

  DotPlan plan;
  plan.mKind = DotSnapshotKind::FINAL;
  plan.mCompId = 7;
  plan.mSubgraphName = "test_empty";
  // No edge layers added.

  auto const dot = SerializeToDotString(tgraph.mNodes, plan);

  // No `color="#` should appear since neither walks nor probe-strip overlays
  // were added — only the per-node fillcolor literals.
  CHECK(dot.find(R"raw( [color=")raw") == std::string::npos);
}

// ============================================================================
//  Renderer — node layer composition (anchor stacks on role)
// ============================================================================

TEST_CASE("Renderer stacks anchor border above role fillcolor",
          "[lancet][cbdg][Renderer][NodeLayers]") {
  TestGraph tgraph;
  auto const nid_a = tgraph.AddNode("ACGTACGTACG");
  auto const nid_b = tgraph.AddNode("CGTACGTACGA");
  auto const nid_c = tgraph.AddNode("GTACGTACGAC");
  tgraph.AddEdge(nid_a, nid_b);
  tgraph.AddEdge(nid_b, nid_c);
  tgraph.SetAllComponentId(7);

  DotPlan plan;
  plan.mKind = DotSnapshotKind::FINAL;
  plan.mCompId = 7;
  plan.mSubgraphName = "test_stack";
  plan.mNodeLayers.emplace_back(MakeAnchorLayer(NodeIDPair{nid_a, nid_c}));

  auto const dot = SerializeToDotString(tgraph.mNodes, plan);

  // Both anchors must have the anchor-layer attributes attached. The
  // renderer emits `peripheries=2` only when set by a layer.
  CHECK(CountOccurrences(dot, "peripheries=2") == 2);
  CHECK(dot.find("goldenrod") != std::string::npos);
  // Non-anchor B keeps its role fillcolor (lightblue, the REF default).
  // Both anchors also keep their role fillcolor — the anchor layer touches
  // a different attribute slot.
  CHECK(dot.find("fillcolor=\"lightblue\"") != std::string::npos);
}

// ============================================================================
//  DotSnapshotBuffer lifecycle
// ============================================================================

TEST_CASE("DotSnapshotBuffer Discard drops pending entries", "[lancet][cbdg][SnapshotBuffer]") {
  DotSnapshotBuffer buf;
  CHECK(buf.IsEmpty());

  buf.Buffer("a.dot", "graph A {}");
  buf.Buffer("b.dot", "graph B {}");
  CHECK_FALSE(buf.IsEmpty());

  buf.Discard();
  CHECK(buf.IsEmpty());
}

TEST_CASE("DotSnapshotBuffer Commit flushes all entries to disk",
          "[lancet][cbdg][SnapshotBuffer]") {
  // Catch2 runs tests sequentially per process; no contention on a fixed name.
  auto const tmp_dir = std::filesystem::temp_directory_path() / "lancet_dot_buffer_test";
  std::filesystem::remove_all(tmp_dir);

  DotSnapshotBuffer buf;
  buf.Buffer("first.dot", "first contents\n");
  buf.Buffer("second.dot", "second contents\n");

  buf.Commit(tmp_dir);
  CHECK(buf.IsEmpty());

  // Both files should exist and contain the buffered bytes verbatim.
  auto const read_all = [](std::filesystem::path const& path) -> std::string {
    std::ifstream const handle(path);
    std::stringstream buf;
    buf << handle.rdbuf();
    return buf.str();
  };
  CHECK(read_all(tmp_dir / "first.dot") == "first contents\n");
  CHECK(read_all(tmp_dir / "second.dot") == "second contents\n");

  std::filesystem::remove_all(tmp_dir);
}
