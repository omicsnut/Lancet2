#ifndef SRC_LANCET_CBDG_PROBE_TRACKER_H_
#define SRC_LANCET_CBDG_PROBE_TRACKER_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/read.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <array>
#include <filesystem>
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

/// Human-readable names for each pruning stage (log output and lost_at_stage TSV values).
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
// ProbeKRecord: per-probe, per-k lifecycle record.
//
// Tracks the survival count of ALT-unique k-mers at each pruning stage,
// plus boolean flags for structural failure modes (anchor overlap, cycle
// retry, traversal limit, etc.). One record is created for each
// (probe_id, k) pair encountered during assembly.
//
// Extended fields (offset analysis, MSA extraction, genotyper read
// assignment) are populated by ProbeDiagnostics through the
// FindFinalKRecord accessor. All fields remain in this struct so
// WriteResults can emit a single unified TSV line.
// ============================================================================
struct ProbeKRecord {
  // ── 8B Align ────────────────────────────────────────────────────────────
  usize mKmerSize = 0;      // 8B — kmer length for this assembly attempt
  usize mCompId = 0;        // 8B — component ID the probe's nodes are in
  usize mCompNumNodes = 0;  // 8B — number of nodes in that component

  /// Indices of haplotype paths containing the variant's ALT context.
  /// Empty means variant was not found in any enumerated path.
  absl::InlinedVector<u8, 8> mHapIndices;  // 8B inline — haplotype indices carrying variant

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
  bool mIsTraversalLimited = false;  // 1B — BFS budget exhausted during path enum
  bool mIsCycleRetry = false;        // 1B — cycle detected, k was retried
  bool mIsComplexRetry = false;      // 1B — graph too complex, k was retried
  bool mIsNoAnchor = false;          // 1B — component had no valid source/sink
  bool mIsShortAnchor = false;       // 1B — anchor region too short (<150bp)
  bool mIsVariantInAnchor = false;   // 1B — variant overlaps source/sink k-mer range

  // MSA extraction flags
  bool mIsMsaExactMatch = false;      // 1B — exact position + allele match in VariantSet
  bool mIsMsaShifted = false;         // 1B — same allele but shifted position
  bool mIsMsaRepresentation = false;  // 1B — variant subsumed by larger MNV representation

  // Genotyper outcome flags
  bool mIsGenoHasAltSupport = false;  // 1B — genotyper found >0 reads supporting truth ALT
  bool mIsGenoNoResult = false;       // 1B — variant not found in genotyper result map
};

// ============================================================================
// ProbeTracker: zero-overhead diagnostic for tracing variant k-mer survival.
//
// Architecture:
//   - Always compiled, no #ifdef. All methods are no-ops when mVariants is
//     empty (i.e., --probe-variants was not specified on the CLI).
//   - Owned by ProbeDiagnostics, which provides a ProbeTracker* to
//     Graph via SetProbeTracker. Graph holds a non-owning pointer.
//   - Node tagging uses a side-table (mNodeTags) keyed by NodeID, avoiding
//     any modifications to the Node class itself.
//
// Lifecycle per window × k:
//   1. GenerateAndTag — clear side-table, find probes in this window,
//      generate ALT k-mers at current k, tag matching nodes.
//   2. LogStatus("after_build") — count surviving tags.
//   3. Graph mutations (RemoveNode, CompressNode) trigger OnNodeRemove /
//      OnNodeMerge to keep the side-table synchronized.
//   4. LogStatus at each pruning boundary.
//   5. CheckAnchorOverlap — flag probes whose variant overlaps anchors.
//   6. CheckPaths — substring-check for variant in enumerated haplotypes.
//   7. WriteResults — emit probe_results.tsv at pipeline completion.
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

  /// Flush accumulated probe records to the results TSV on destruction (RAII).
  ~ProbeTracker() { WriteResults(); }

  /// Non-copyable: accumulator semantics — copying would duplicate records and
  /// double-write on destruction. Movable: moved-from containers are empty, so
  /// the destructor's WriteResults() safely no-ops (early return on empty records).
  ProbeTracker(ProbeTracker const&) = delete;
  auto operator=(ProbeTracker const&) -> ProbeTracker& = delete;
  ProbeTracker(ProbeTracker&&) noexcept = default;
  auto operator=(ProbeTracker&&) noexcept -> ProbeTracker& = default;

  /// Load truth variants from a pre-parsed vector.
  void LoadVariants(std::vector<ProbeVariant> variants);

  /// Set the output path for probe results TSV.
  void SetResultsPath(std::filesystem::path results_path);

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
  void LogStatus(PruneStage stage, Context const& ctx);

  /// Flag probes whose variant position overlaps the source/sink anchor k-mers.
  void CheckAnchorOverlap(AnchorInfo const& source, AnchorInfo const& sink, Context const& ctx);

  /// Check if any enumerated haplotype path contains each probe's ALT context.
  void CheckPaths(absl::Span<Path const> haplotypes, Context const& ctx);

  /// Flag the current k-attempt as a cycle/complex/no-anchor/short-anchor retry.
  void SetCycleRetry(Context const& ctx);
  void SetComplexRetry(Context const& ctx);
  void SetNoAnchor(Context const& ctx);
  void SetShortAnchor(Context const& ctx);

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

  /// Record for the final k-value assembly attempt for a given probe.
  /// Searches records in reverse order (highest k first) for the last entry.
  [[nodiscard]] auto FindFinalKRecord(u16 probe_id) -> ProbeKRecord&;

  /// Per-probe read ownership index (populated by CountInReads).
  /// Maps probe_id → set of qname hashes that carry ALT-unique k-mers.
  [[nodiscard]] auto ReadOwnershipIndex() const noexcept
      -> absl::flat_hash_map<u16, absl::flat_hash_set<u32>> const& {
    return mProbeReadIndex;
  }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<ProbeVariant> mVariants;
  absl::flat_hash_map<NodeID, absl::InlinedVector<ProbeHit, 2>> mNodeTags;
  std::vector<ProbeKRecord> mRecords;
  std::vector<u16> mActiveProbeIds;
  absl::flat_hash_map<u16, u16> mAltKmerCounts;
  std::filesystem::path mResultsPath;
  std::shared_ptr<ProbeIndex const> mProbeIndex;
  absl::flat_hash_map<u16, absl::flat_hash_set<u32>> mProbeReadIndex;
  absl::flat_hash_set<u16> mActiveProbeIdSet;

  /// Write all accumulated probe records to the results TSV file.
  /// Called automatically by the destructor (RAII).
  void WriteResults() const;

  /// Find or create a ProbeKRecord for the given (probe_id, kmer_size) pair.
  [[nodiscard]] auto FindOrCreateRecord(u16 probe_id, usize kmer_size) -> ProbeKRecord&;

  /// Count surviving k-mers per probe from the side-table with offset analysis.
  [[nodiscard]] auto CountSurvivingKmers() const -> absl::flat_hash_map<u16, SurvivalProfile>;

  /// Derive the lost_at attribution string from a record's state.
  [[nodiscard]] static auto DeriveLostAt(ProbeKRecord const& record) -> std::string_view;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_PROBE_TRACKER_H_
