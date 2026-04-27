#ifndef SRC_LANCET_CBDG_PROBE_TRACKER_H_
#define SRC_LANCET_CBDG_PROBE_TRACKER_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/probe_results_writer.h"
#include "lancet/cbdg/read.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// PruneStage: typed identifier for each graph pruning boundary.
//
// Used by ProbeTracker::LogStatus to record the surviving k-mer count at
// each stage. Replaces string-based dispatch with a compile-time enum.
// ============================================================================
enum class PruneStage : u8 {
  PRUNED_AT_BUILD = 0,
  PRUNED_AT_LOWCOV1,
  PRUNED_AT_COMPRESS1,
  PRUNED_AT_LOWCOV2,
  PRUNED_AT_COMPRESS2,
  PRUNED_AT_TIPS,
};

static constexpr usize NUM_PRUNE_STAGES = 6;

/// Human-readable names for each pruning stage (log output and survival counting).
static constexpr std::array<std::string_view, NUM_PRUNE_STAGES> PRUNE_STAGE_NAMES = {{
    "pruned_at_build",
    "pruned_at_lowcov1",
    "pruned_at_compress1",
    "pruned_at_lowcov2",
    "pruned_at_compress2",
    "pruned_at_tips",
}};

// ============================================================================
// ProbeHit: a single k-mer tag stored in the node side-table.
//
// Records which probe variant and which positional offset within that
// probe's ALT k-mer chain this node carries. Multiple probes can tag
// the same node if their ALT k-mers share a canonical hash.
// ============================================================================
struct ProbeHit {
  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mProbeId = 0;           // 2B — which probe variant this tag belongs to
  u16 mAltContextOffset = 0;  // 2B — 0-indexed position within the probe's ALT context k-mer chain
};

// ============================================================================
// SurvivalProfile: offset-aware k-mer survival statistics per probe.
//
// Extends the simple surviving count with positional analysis of which
// k-mers in the ALT context chain survived: edge k-mers (offset 0 or max),
// interior k-mers, and contiguous gaps. This enables distinguishing between
// uniform attrition (random pruning) and systematic edge erosion (anchor
// overlap, chain-end effects).
// ============================================================================
struct SurvivalProfile {
  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mCount = 0;          // 2B — total surviving k-mers
  u16 mEdgeCount = 0;      // 2B — surviving at chain edges (offset 0 or max)
  u16 mInteriorCount = 0;  // 2B — surviving in chain interior
  u16 mGapCount = 0;       // 2B — contiguous gaps in surviving offset sequence
};

// ============================================================================
// ProbeKRecord: per-probe, per-attempt lifecycle record.
//
// Each record captures the complete observed state for a single assembly
// attempt, uniquely identified by the 4-tuple (probe_id, window, comp_id, k).
// This granularity isolates each (window × component × k-value) attempt so
// that structural failures in one component do not contaminate records from
// another component.
//
// Record lifecycle:
//   1. Created by FindOrCreateRecord when any probe tracking method first
//      accesses a (probe_id, region, comp_id, k) combination.
//   2. Survival counts populated by LogStatus at each pruning stage.
//   3. Structural flags set by SetNoAnchor/SetGraphCycle/etc., filtered
//      to only probes with tagged nodes in the failing component.
//   4. Haplotype indices set by CheckPaths if variant found in walks.
//   5. MSA/genotyper fields populated by ProbeDiagnostics, which annotates
//      ALL records with non-empty mHapIndices via FindRecordsWithPaths.
//   6. Written to TSV as raw facts — no attribution logic in C++.
//
// CRITICAL: comp_id = 0 is the sentinel for pre-component stages
// (PRUNED_AT_BUILD, PRUNED_AT_LOWCOV1) that run before
// MarkConnectedComponents. All real component IDs are >= 1.
// ============================================================================
struct ProbeKRecord {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::string mRegion;      // 24B — window region string (e.g. "chr1:1000-2000"), owned
  usize mKmerSize = 0;      // 8B  — kmer length for this assembly attempt
  usize mCompId = 0;        // 8B  — component ID (0 = pre-component sentinel, >= 1 = real)
  usize mCompNumNodes = 0;  // 8B  — number of nodes in that component

  /// Indices of haplotype paths containing the variant's ALT context.
  /// Empty means variant was not found in any enumerated path.
  /// InlinedVector<u8, 8>: 8 inline slots fit within minimum InlinedVector
  /// footprint (24B) — no wasted space. Typical count is 1–3 paths.
  absl::InlinedVector<u8, 8> mHapIndices;  // 24B — haplotype indices carrying variant

  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mProbeId = 0;             // 2B — index into the variants vector
  u16 mExpectedAltKmers = 0;    // 2B — expected ALT-unique k-mers at this k
  u16 mAltKmersInReads = 0;     // 2B — ALT-unique k-mers found in the raw read set
  u16 mSurvivingBuild = 0;      // 2B — surviving k-mers after graph construction
  u16 mSurvivingLowcov1 = 0;    // 2B — surviving after global low-cov removal
  u16 mSurvivingCompress1 = 0;  // 2B — surviving after first compression
  u16 mSurvivingLowcov2 = 0;    // 2B — surviving after second low-cov removal
  u16 mSurvivingCompress2 = 0;  // 2B — surviving after second compression
  u16 mSurvivingTips = 0;       // 2B — surviving after tip removal
  u16 mSplitAcrossComps = 0;    // 2B — >0 means probe split across components

  // Offset analysis — positional survival profile at last pruning stage
  u16 mLastEdgeCount = 0;      // 2B — surviving edge k-mers (offset 0 or max)
  u16 mLastInteriorCount = 0;  // 2B — surviving interior k-mers
  u16 mLastGapCount = 0;       // 2B — contiguous gaps in surviving offsets

  // MSA extraction classification
  i16 mMsaShiftBp = 0;  // 2B — signed shift from truth position (0 = exact or unmatched)

  // Genotyper read assignment analysis
  u16 mGenoTrueAltReads = 0;      // 2B — reads correctly assigned to truth ALT allele
  u16 mGenoTotalRefReads = 0;     // 2B — total reads assigned to REF at this variant
  u16 mGenoStolenToRef = 0;       // 2B — ALT-carrying reads misassigned to REF
  u16 mGenoStolenToWrongAlt = 0;  // 2B — ALT-carrying reads misassigned to wrong ALT
  u16 mGenoNonOverlapping = 0;    // 2B — ALT-carrying reads that didn't overlap alignment

  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mIsTraversalLimited = false;  // 1B — BFS walk-tree budget exhausted during path enum
  bool mIsGraphCycle = false;        // 1B — de Bruijn graph had a cycle in this component
  bool mIsGraphComplex = false;      // 1B — graph exceeded complexity threshold for this component
  bool mIsNoAnchor = false;          // 1B — component had no valid source/sink anchor pair
  bool mIsShortAnchor = false;       // 1B — anchor region too short (<150bp)
  bool mIsVariantInAnchor = false;   // 1B — overlaps source/sink k-mer range (diagnostic only)

  // MSA extraction flags
  bool mIsMsaExactMatch = false;  // 1B — exact position + allele match in VariantSet
  bool mIsMsaShifted = false;     // 1B — same allele but shifted position
  bool mIsMsaSubsumed = false;    // 1B — alleles absorbed into a larger MNV by MSA alignment

  // Genotyper outcome flags
  bool mIsGenoHasAltSupport = false;  // 1B — genotyper found >0 reads supporting truth ALT
  bool mIsGenoNoOverlap = false;  // 1B — zero read alignments overlapped variant's haplotype coords
  bool mIsNotProcessed = false;   // 1B — probe never activated in any window
};

// ============================================================================
// ProbeTracker: zero-overhead diagnostic for tracing variant k-mer survival.
//
// Pure data emitter — captures every observed fact per (probe_id, window,
// comp_id, k) attempt. No attribution logic. All stage attribution is
// performed downstream in the Python analysis script.
//
// Architecture:
//   - Always compiled, no #ifdef. All methods are no-ops when mVariants is
//     empty (i.e., --probe-variants was not specified on the CLI).
//   - Owned by ProbeDiagnostics, which provides a ProbeTracker* to
//     Graph via SetProbeTracker. Graph holds a non-owning pointer.
//   - Node tagging uses a side-table (mNodeTags) keyed by NodeID, avoiding
//     any modifications to the Node class itself.
//
// Record key: (probe_id, window, comp_id, k)
//   A single probe may generate multiple records across overlapping windows,
//   different connected components, and k-value iterations. Each record
//   independently tracks survival counts, structural flags, path matches,
//   MSA extraction results, and genotyper outcomes.
//
// Lifecycle per window × k:
//   1. GenerateAndTag — clear side-table, find probes in this window,
//      generate ALT k-mers at current k, tag matching nodes.
//   2. CountInReads — count ALT-unique k-mers in raw reads.
//   3. LogStatus(PRUNED_AT_BUILD) — count surviving tags (comp_id=0).
//   4. Graph mutations (RemoveNode, CompressNode) trigger OnNodeRemove /
//      OnNodeMerge to keep the side-table synchronized.
//   5. LogStatus at each pruning boundary (comp_id from Context).
//   6. CheckAnchorOverlap — flag probes whose variant overlaps anchors.
//   7. CheckPaths — substring-check for variant in enumerated haplotypes.
//   8. SubmitCompleted — flush records to writer after each window.
// ============================================================================
class ProbeTracker {
 public:
  using NodeTable = absl::flat_hash_map<NodeID, std::unique_ptr<Node>>;

  /// Component metadata — matches Graph::ComponentInfo layout.
  struct ComponentInfo {
    // ── 8B Align ──────────────────────────────────────────────────────────
    usize mCompId = 0;    // 8B
    usize mNumNodes = 0;  // 8B
  };

  /// Source/sink anchor coordinates for anchor overlap detection.
  struct AnchorInfo {
    // ── 8B Align ──────────────────────────────────────────────────────────
    usize mRefOffset = 0;  // 8B — reference offset of the anchor k-mer
    // ── 1B Align ──────────────────────────────────────────────────────────
    bool mFound = false;  // 1B — whether the anchor was discovered
  };

  /// Assembly context shared across probe tracking methods.
  /// Constructed once per k-iteration and updated per component.
  struct Context {
    // ── 8B Align ────────────────────────────────────────────────────────────
    std::string_view mChrom;   // 16B — chromosome name
    std::string_view mRefSeq;  // 16B — reference sequence for the window
    std::string_view mRegStr;  // 16B — region string (e.g. "chr1:100-200")
    usize mRegionStart = 0;    // 8B  — 0-based region start position
    usize mKmerSize = 0;       // 8B  — current k-mer length
    usize mCompId = 0;         // 8B  — current component ID
  };

  ProbeTracker() = default;

  /// Submit any remaining probe records to the results writer on destruction (safety net).
  ~ProbeTracker() { SubmitCompleted(); }

  /// Non-copyable: accumulator semantics — copying would duplicate records and
  /// double-write on destruction. Movable: moved-from containers are empty, so
  /// the destructor's WriteResults() safely no-ops (early return on empty records).
  ProbeTracker(ProbeTracker const&) = delete;
  auto operator=(ProbeTracker const&) -> ProbeTracker& = delete;
  ProbeTracker(ProbeTracker&&) noexcept = default;
  auto operator=(ProbeTracker&&) noexcept -> ProbeTracker& = default;

  /// Load truth variants from a pre-parsed vector.
  void LoadVariants(std::vector<ProbeVariant> variants);

  /// Set the shared results writer for thread-safe TSV output.
  void SetResultsWriter(std::shared_ptr<ProbeResultsWriter> writer);

  /// Flush current records to the writer. Called after each window and from destructor.
  void SubmitCompleted();

  /// Set the precomputed global k-mer index (shared across threads).
  void SetProbeIndex(std::shared_ptr<ProbeIndex const> probe_index);

  /// Clear side-table, find probes in this window, generate ALT k-mers
  /// at the given kmer_size, and tag matching nodes in the graph.
  void GenerateAndTag(NodeTable const& nodes, Context const& ctx);

  /// Count how many ALT-unique k-mers appear in the raw read set (before graph
  /// construction). Uses the precomputed ProbeIndex for O(1) lookups. Also
  /// builds the per-probe read ownership index (mProbeReadIndex) for downstream
  /// stolen-read analysis.
  void CountInReads(absl::Span<Read const> reads, Context const& ctx);

  /// After MarkConnectedComponents: determine which component each probe's
  /// tagged nodes are in. Detects fragmentation (k-mers split across components).
  void RecordComponentInfo(NodeTable const& nodes, absl::Span<ComponentInfo const> components,
                           Context const& ctx);

  /// Transfer probe tags from absorbed_id to surviving_id during node compression.
  void OnNodeMerge(NodeID absorbed_id, NodeID surviving_id);

  /// Remove probe tags for a deleted node.
  void OnNodeRemove(NodeID node_id);

  /// Count surviving tagged k-mers per probe and log + record at a pipeline stage.
  /// Requires the node table for per-component filtering of survival counts.
  void LogStatus(PruneStage stage, Context const& ctx, NodeTable const& nodes);

  /// Flag probes whose variant position overlaps the source/sink anchor k-mers.
  void CheckAnchorOverlap(AnchorInfo const& source, AnchorInfo const& sink, Context const& ctx);

  /// Check if any enumerated haplotype path contains each probe's ALT context.
  void CheckPaths(absl::Span<Path const> haplotypes, Context const& ctx);

  /// Flag the current k-attempt as a graph cycle / graph too complex.
  void SetGraphCycle(Context const& ctx);
  void SetGraphComplex(Context const& ctx);

  /// Flag probes with tagged nodes in this component as having no valid anchor.
  /// Checks component membership before setting — only probes with tagged nodes
  /// in ctx.mCompId are flagged, preventing cross-component contamination.
  void SetNoAnchor(Context const& ctx, NodeTable const& nodes);
  void SetShortAnchor(Context const& ctx, NodeTable const& nodes);

  /// Mark all active probes with the traversal limit flag.
  void SetTraversalLimit(Context const& ctx);

  /// True when probe variants have been loaded (--probe-variants was specified).
  [[nodiscard]] auto IsActive() const noexcept -> bool { return !mVariants.empty(); }

  /// Return NodeIDs of tagged nodes in a specific component for DOT highlighting.
  [[nodiscard]] auto GetHighlightNodeIds(NodeTable const& nodes, usize comp_id) const
      -> absl::flat_hash_set<NodeID>;

  // ── Public accessors for ProbeDiagnostics ─────────────────────────────

  /// Read access to loaded truth variants.
  [[nodiscard]] auto Variants() const noexcept -> absl::Span<ProbeVariant const> {
    return absl::MakeConstSpan(mVariants);
  }

  /// Active probe IDs for the current window.
  [[nodiscard]] auto ActiveProbeIds() const noexcept -> absl::Span<u16 const> {
    return absl::MakeConstSpan(mActiveProbeIds);
  }

  /// All records for this probe that have non-empty haplotype path indices.
  /// Returns pointers to each matching record for MSA/genotyper annotation.
  /// InlinedVector<*, 4>: typical probes succeed in ≤4 (window × component)
  /// combinations, so 4 inline slots avoid heap allocation in the common case.
  [[nodiscard]] auto FindRecordsWithPaths(u16 probe_id) -> absl::InlinedVector<ProbeKRecord*, 4>;

  /// Per-probe read ownership index (populated by CountInReads).
  /// Maps probe_id → set of qname hashes that carry ALT-unique k-mers.
  [[nodiscard]] auto ReadOwnershipIndex() const noexcept
      -> absl::flat_hash_map<u16, absl::flat_hash_set<u32>> const& {
    return mProbeReadIndex;
  }

  /// True if this probe has any tagged graph nodes in the specified component.
  /// Used by SetNoAnchor/SetShortAnchor to filter structural flags to only
  /// probes that actually participate in the failing component.
  [[nodiscard]] auto HasTaggedNodesInComponent(u16 probe_id, usize comp_id,
                                               NodeTable const& nodes) const -> bool;

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<ProbeVariant> mVariants;
  absl::flat_hash_map<NodeID, absl::InlinedVector<ProbeHit, 2>> mNodeTags;
  std::vector<ProbeKRecord> mRecords;
  std::vector<u16> mActiveProbeIds;
  absl::flat_hash_map<u16, u16> mAltKmerCounts;
  std::shared_ptr<ProbeResultsWriter> mResultsWriter;
  std::shared_ptr<ProbeIndex const> mProbeIndex;
  absl::flat_hash_map<u16, absl::flat_hash_set<u32>> mProbeReadIndex;
  absl::flat_hash_set<u16> mActiveProbeIdSet;

  // ============================================================================
  // FindOrCreateRecord: 4-key lookup for probe lifecycle records.
  //
  // Searches mRecords for a record matching (probe_id, region, comp_id, kmer_size).
  // If found, returns a reference to the existing record. If not found, creates
  // a new record with those key fields and appends it to mRecords.
  //
  // O(n) linear scan — acceptable for typical probe counts (<100 per window).
  // If profiling shows a hot spot, switch to a hash map of indices.
  // ============================================================================
  [[nodiscard]] auto FindOrCreateRecord(u16 probe_id, std::string_view region, usize comp_id,
                                        usize kmer_size) -> ProbeKRecord&;

  /// Count surviving k-mers per probe from the side-table with offset analysis.
  /// When nodes is non-null, only counts tags on nodes belonging to comp_id.
  [[nodiscard]] auto CountSurvivingKmers(NodeTable const* nodes, usize comp_id) const
      -> absl::flat_hash_map<u16, SurvivalProfile>;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PROBE_TRACKER_H_
