#include "lancet/cbdg/probe_tracker.h"

#include "lancet/base/logging.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/path.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ranges.h"

#include <algorithm>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {
namespace {

// ============================================================================
// VariantFlankRange: the reference coordinate window around a variant that
// captures every k-mer overlapping any base of the REF allele.
//
// Geometry (non-boundary case, ref_allele_len = R, k = kmer_size):
//
//   ref_seq:  ... [ctx_start          var_offset     var_offset+R        ctx_end) ...
//                  |<--- k-1 bases --->|<--- R bases --->|<--- k-1 bases --->|
//                  |                   |                 |                   |
//   leftmost k-mer:|===================|                 |  (starts at ctx_start,
//                  |                   |                 |   last base at var_offset)
//                  |                   |  rightmost k-mer:|==================|
//                  |                   |  (starts at var_offset+R-1,
//                  |                   |   last base at ctx_end-1)
//
// REF context length: (k-1) + R + (k-1) = 2k + R - 2
// ALT context length: (k-1) + A + (k-1) = 2k + A - 2  (A = alt_allele_len)
//
// Examples with k=25:
//   SNV  (R=1, A=1):  REF = 49bp → 25 k-mers,  ALT = 49bp → 25 k-mers
//   DEL  (R=3, A=1):  REF = 51bp → 27 k-mers,  ALT = 49bp → 25 junction k-mers
//   INS  (R=1, A=4):  REF = 49bp → 25 k-mers,  ALT = 52bp → 28 k-mers
//   BigDEL (R=500):   REF = 548bp → 524 k-mers, ALT = 49bp → 25 breakpoint k-mers
//   BigINS (A=500):   REF = 49bp → 25 k-mers,  ALT = 548bp → 524 k-mers
//
// Boundary clipping: ctx_start is clamped to 0, ctx_end to ref_seq.length().
// If the resulting context is shorter than kmer_size, no k-mers are produced
// (callers guard with `if (context.length() >= kmer_size)`).
// ============================================================================
struct VariantFlankRange {
  usize mCtxStart = 0;
  usize mCtxEnd = 0;
};

[[nodiscard]] auto ComputeFlankRange(usize ref_seq_len, usize var_offset, usize ref_allele_len,
                                     usize kmer_size) -> VariantFlankRange {
  auto const flank_len = kmer_size - 1;
  return {
      .mCtxStart = var_offset > flank_len ? var_offset - flank_len : usize{0},
      .mCtxEnd = std::min(ref_seq_len, var_offset + ref_allele_len + flank_len),
  };
}

// ============================================================================
// BuildAltContext: splice a variant's ALT allele into the reference flanking
// context, producing the ALT haplotype substring for k-mer extraction.
//
// Output: ref_seq[ctx_start..var_offset] + alt_allele + ref_seq[var_offset+R..ctx_end]
// Length: left_flank + alt_allele_len + right_flank
// Each k-mer from SlidingView(result, kmer_size) is exactly kmer_size bases.
// ============================================================================
[[nodiscard]] auto BuildAltContext(std::string_view ref_seq, usize var_offset, usize ref_allele_len,
                                   std::string_view alt_allele, usize kmer_size) -> std::string {
  auto const [ctx_start, ctx_end] =
      ComputeFlankRange(ref_seq.length(), var_offset, ref_allele_len, kmer_size);

  std::string alt_context;
  alt_context.reserve(ctx_end - ctx_start + alt_allele.length());
  alt_context.append(ref_seq.substr(ctx_start, var_offset - ctx_start));
  alt_context.append(alt_allele);

  auto const right_flank_start = var_offset + ref_allele_len;
  if (right_flank_start < ctx_end) {
    auto const right_flank_len = ctx_end - right_flank_start;
    alt_context.append(ref_seq.substr(right_flank_start, right_flank_len));
  }

  return alt_context;
}

// ============================================================================
// CollectRefKmerHashes: hash all canonical k-mers from the reference flanking
// context around a variant. Used to identify ALT-unique k-mers by set
// difference with the ALT context k-mers.
//
// Context: ref_seq[ctx_start..ctx_end], length = 2k + R - 2 (non-boundary).
// Each k-mer from SlidingView is exactly kmer_size bases.
// ============================================================================
[[nodiscard]] auto CollectRefKmerHashes(std::string_view ref_seq, usize var_offset,
                                        usize ref_allele_len, usize kmer_size)
    -> absl::flat_hash_set<u64> {
  auto const [ctx_start, ctx_end] =
      ComputeFlankRange(ref_seq.length(), var_offset, ref_allele_len, kmer_size);

  auto const ref_context = ref_seq.substr(ctx_start, ctx_end - ctx_start);
  absl::flat_hash_set<u64> ref_hashes;
  if (ref_context.length() < kmer_size) return ref_hashes;

  auto const sliding = lancet::base::SlidingView(ref_context, kmer_size);
  ref_hashes.reserve(sliding.size());
  std::ranges::transform(sliding, std::inserter(ref_hashes, ref_hashes.end()),
                         [](auto const& kmer_seq) { return Kmer(kmer_seq).Identifier(); });
  return ref_hashes;
}

// ============================================================================
// ReadContainsKmer: pure predicate — does a single read contain a k-mer
// matching the given canonical hash? Uses std::ranges::any_of over the
// sliding window of exactly kmer_size bases.
// ============================================================================
[[nodiscard]] auto ReadContainsKmer(Read const& read, u64 kmer_hash, usize kmer_size) -> bool {
  if (read.SeqView().length() < kmer_size) return false;
  return std::ranges::any_of(
      lancet::base::SlidingView(read.SeqView(), kmer_size),
      [kmer_hash](std::string_view mer_seq) { return Kmer(mer_seq).Identifier() == kmer_hash; });
}

// ============================================================================
// CountTaggedNodesPerComponent: for a given probe, count how many of its
// tagged graph nodes landed in each connected component.
// ============================================================================
using NodeTagMap = absl::flat_hash_map<NodeID, absl::InlinedVector<ProbeHit, 2>>;
using NodeTable = absl::flat_hash_map<NodeID, std::unique_ptr<Node>>;

[[nodiscard]] auto CountTaggedNodesPerComponent(NodeTagMap const& node_tags, NodeTable const& nodes,
                                                u16 probe_id) -> absl::flat_hash_map<usize, usize> {
  absl::flat_hash_map<usize, usize> counts;
  for (auto const& [node_id, hits] : node_tags) {
    auto const matches = [probe_id](ProbeHit const& hit) { return hit.mProbeId == probe_id; };
    if (!std::ranges::any_of(hits, matches)) continue;

    auto const node_it = nodes.find(node_id);
    if (node_it == nodes.end()) continue;
    counts[node_it->second->GetComponentId()]++;
  }
  return counts;
}

// ============================================================================
// LookupComponentSize: find a component's total node count by ID.
// Returns 0 if the component is not found in the span.
// ============================================================================
[[nodiscard]] auto LookupComponentSize(absl::Span<ProbeTracker::ComponentInfo const> components,
                                       usize comp_id) -> usize {
  auto const* const comp_it =
      std::ranges::find_if(components, [comp_id](ProbeTracker::ComponentInfo const& comp_info) {
        return comp_info.mCompId == comp_id;
      });

  return comp_it != components.end() ? comp_it->mNumNodes : usize{0};
}

}  // namespace

void ProbeTracker::LoadVariants(std::vector<ProbeVariant> variants) {
  mVariants = std::move(variants);
}

void ProbeTracker::SetResultsWriter(std::shared_ptr<ProbeResultsWriter> writer) {
  mResultsWriter = std::move(writer);
}

void ProbeTracker::SubmitCompleted() {
  if (mRecords.empty() || !mResultsWriter) return;
  mResultsWriter->Append(absl::MakeConstSpan(mRecords));
  mRecords.clear();
}

void ProbeTracker::SetProbeIndex(std::shared_ptr<ProbeIndex const> probe_index) {
  mProbeIndex = std::move(probe_index);
}

auto ProbeTracker::FindFinalKRecord(u16 const probe_id) -> ProbeKRecord& {
  for (auto& record : std::views::reverse(mRecords)) {
    if (record.mProbeId == probe_id) return record;
  }

  mRecords.emplace_back(ProbeKRecord{.mProbeId = probe_id});
  return mRecords.back();
}

// ============================================================================
// GenerateAndTag: produce ALT-unique k-mers and tag matching graph nodes.
//
// For each probe variant in the current window:
//   1. Build a local ALT sequence by splicing the ALT allele into the
//      reference context (k-1 bases of flanking on each side).
//   2. Extract all k-mers from both REF and ALT contexts.
//   3. K-mers present in ALT but absent from REF are "ALT-unique" —
//      they represent the variant signal in the graph.
//   4. Look up each ALT-unique k-mer in the node table by canonical hash.
//      If found, tag the node in the side-table.
//
// For a SNV with k=25: exactly 25 ALT-unique k-mers (each position in
// the sliding window covers the variant base).
// For INS of length L: up to k+L ALT-unique k-mers.
// For DEL of length D: up to k junction k-mers.
// ============================================================================
void ProbeTracker::GenerateAndTag(NodeTable const& nodes, Context const& ctx) {
  if (mVariants.empty()) return;

  mNodeTags.clear();
  mActiveProbeIds.clear();
  mActiveProbeIdSet.clear();
  mAltKmerCounts.clear();

  auto const region_end = ctx.mRegionStart + ctx.mRefSeq.length();

  for (auto const& probe : mVariants) {
    if (probe.mChrom != ctx.mChrom) continue;
    if (probe.mGenomeStart0 < ctx.mRegionStart || probe.mGenomeStart0 >= region_end) continue;

    mActiveProbeIds.push_back(probe.mProbeId);
    mActiveProbeIdSet.insert(probe.mProbeId);

    auto const var_offset = probe.mGenomeStart0 - ctx.mRegionStart;
    auto const ref_allele_len = probe.mRef.length();
    auto const ref_hashes =
        CollectRefKmerHashes(ctx.mRefSeq, var_offset, ref_allele_len, ctx.mKmerSize);
    auto const alt_context =
        BuildAltContext(ctx.mRefSeq, var_offset, ref_allele_len, probe.mAlt, ctx.mKmerSize);

    // Extract ALT-unique k-mers and tag matching nodes.
    // Early-continue flattens the loop body: skip REF k-mers, then check graph.
    u16 alt_kmer_offset = 0;
    if (alt_context.length() >= ctx.mKmerSize) {
      for (auto const& kmer_seq : lancet::base::SlidingView(alt_context, ctx.mKmerSize)) {
        auto const kmer_hash = Kmer(kmer_seq).Identifier();
        if (ref_hashes.contains(kmer_hash) || !nodes.contains(kmer_hash)) continue;

        mNodeTags[kmer_hash].push_back(
            {.mProbeId = probe.mProbeId, .mAltContextOffset = alt_kmer_offset});
        ++alt_kmer_offset;
      }
    }

    mAltKmerCounts[probe.mProbeId] = alt_kmer_offset;
  }
}

// ============================================================================
// CountInReads: scan the raw read set for ALT-unique k-mers using ProbeIndex.
//
// Uses the precomputed ProbeIndex for O(1) hash lookups instead of re-deriving
// ALT k-mers per window. Also builds the per-probe read ownership index
// (mProbeReadIndex) mapping probe_id → set of qname hashes for reads that
// carry ALT-unique k-mers. This enables S6 stolen-read analysis.
// ============================================================================
void ProbeTracker::CountInReads(absl::Span<Read const> reads, Context const& ctx) {
  if (mVariants.empty() || mActiveProbeIds.empty()) return;
  if (!mProbeIndex) return;
  auto const* kmer_index = mProbeIndex->ForKmerSize(ctx.mKmerSize);
  if (kmer_index == nullptr) return;

  mProbeReadIndex.clear();
  absl::flat_hash_map<u16, u16> kmer_read_counts;

  for (auto const& read : reads) {
    if (read.SeqView().length() < ctx.mKmerSize) continue;
    auto const read_qname_hash = static_cast<u32>(absl::HashOf(read.QnameView()));

    for (auto const& kmer_seq : lancet::base::SlidingView(read.SeqView(), ctx.mKmerSize)) {
      auto const canonical_hash = Kmer(kmer_seq).Identifier();
      auto const iter = kmer_index->find(canonical_hash);
      if (iter == kmer_index->end()) continue;

      for (auto const& entry : iter->second) {
        if (!mActiveProbeIdSet.contains(entry.mProbeId)) continue;
        mProbeReadIndex[entry.mProbeId].insert(read_qname_hash);
        ++kmer_read_counts[entry.mProbeId];
      }
    }
  }

  for (auto const probe_id : mActiveProbeIds) {
    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
    auto const iter = kmer_read_counts.find(probe_id);
    record.mAltKmersInReads = (iter != kmer_read_counts.end()) ? iter->second : u16{0};
  }
}

void ProbeTracker::RecordComponentInfo(NodeTable const& nodes,
                                       absl::Span<ComponentInfo const> components,
                                       Context const& ctx) {
  if (mActiveProbeIds.empty()) return;

  for (auto const probe_id : mActiveProbeIds) {
    auto const comp_counts = CountTaggedNodesPerComponent(mNodeTags, nodes, probe_id);
    if (comp_counts.empty()) continue;

    static constexpr auto BY_COUNT = [](auto const& lhs, auto const& rhs) {
      return lhs.second < rhs.second;
    };
    auto const majority_comp = std::ranges::max_element(comp_counts, BY_COUNT)->first;
    auto const comp_num_nodes = LookupComponentSize(components, majority_comp);

    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
    record.mCompId = majority_comp;
    record.mCompNumNodes = comp_num_nodes;
    record.mSplitAcrossComps =
        comp_counts.size() > 1 ? static_cast<u16>(comp_counts.size()) : u16{0};

    LOG_DEBUG("KMER_PROBE {} kmer_size={} probe={} {}:{} in comp={} ({} nodes) split={}",
              ctx.mRegStr, ctx.mKmerSize, probe_id, mVariants[probe_id].mChrom,
              mVariants[probe_id].mGenomeStart0, majority_comp, comp_num_nodes,
              record.mSplitAcrossComps)
  }
}

void ProbeTracker::OnNodeMerge(NodeID absorbed_id, NodeID surviving_id) {
  if (mNodeTags.empty()) return;

  auto absorbed_it = mNodeTags.find(absorbed_id);
  if (absorbed_it == mNodeTags.end()) return;

  // Move tags out before modifying the map. operator[] on surviving_id may
  // insert a new entry, triggering a rehash that invalidates absorbed_it.
  auto absorbed_tags = std::move(absorbed_it->second);
  mNodeTags.erase(absorbed_it);

  auto& surviving_tags = mNodeTags[surviving_id];
  surviving_tags.insert(surviving_tags.end(), absorbed_tags.begin(), absorbed_tags.end());
}

void ProbeTracker::OnNodeRemove(NodeID node_id) {
  if (mNodeTags.empty()) return;
  mNodeTags.erase(node_id);
}

void ProbeTracker::LogStatus(PruneStage stage, Context const& ctx) {
  if (mActiveProbeIds.empty()) return;

  static constexpr std::array<u16 ProbeKRecord::*, NUM_PRUNE_STAGES> STAGE_FIELDS = {{
      &ProbeKRecord::mSurvivingBuild,
      &ProbeKRecord::mSurvivingLowcov1,
      &ProbeKRecord::mSurvivingCompress1,
      &ProbeKRecord::mSurvivingLowcov2,
      &ProbeKRecord::mSurvivingCompress2,
      &ProbeKRecord::mSurvivingTips,
  }};

  auto const stage_name = PRUNE_STAGE_NAMES[static_cast<usize>(stage)];
  auto const profiles = CountSurvivingKmers();

  for (auto const probe_id : mActiveProbeIds) {
    auto const prof_it = profiles.find(probe_id);
    auto const found = prof_it != profiles.end() ? prof_it->second.mCount : u16{0};
    auto const total_it = mAltKmerCounts.find(probe_id);
    auto const total = total_it != mAltKmerCounts.end() ? total_it->second : u16{0};
    auto const& probe = mVariants[probe_id];
    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);

    LOG_DEBUG("KMER_PROBE {} comp={} kmer_size={} probe={} {}:{}:{}>{} found={}/{} stage={} "
              "probe_comp={} comp_nodes={}",
              ctx.mRegStr, ctx.mCompId, ctx.mKmerSize, probe_id, probe.mChrom, probe.mGenomeStart0,
              probe.mRef, probe.mAlt, found, total, stage_name, record.mCompId,
              record.mCompNumNodes)

    record.mExpectedAltKmers = total;
    record.*STAGE_FIELDS[static_cast<usize>(stage)] = found;

    if (prof_it != profiles.end()) {
      record.mLastEdgeCount = prof_it->second.mEdgeCount;
      record.mLastInteriorCount = prof_it->second.mInteriorCount;
      record.mLastGapCount = prof_it->second.mGapCount;
    }
  }
}

// ============================================================================
// CheckAnchorOverlap: flag probes whose variant position falls within the
// source or sink anchor k-mer's genomic range.
//
// A variant overlapping an anchor is a structural impossibility — the anchor
// is the fixed reference endpoint of the walk, so the ALT allele can never
// appear in an assembled path. This is distinct from a pruning loss.
//
//   Source k-mer range: [region_start + source.mRefOffset,
//                        region_start + source.mRefOffset + k)
//   Sink k-mer range:   [region_start + sink.mRefOffset,
//                        region_start + sink.mRefOffset + k)
// ============================================================================
void ProbeTracker::CheckAnchorOverlap(AnchorInfo const& source, AnchorInfo const& sink,
                                      Context const& ctx) {
  if (mActiveProbeIds.empty()) return;

  auto const src_start = ctx.mRegionStart + source.mRefOffset;
  auto const src_end = src_start + ctx.mKmerSize;
  auto const snk_start = ctx.mRegionStart + sink.mRefOffset;
  auto const snk_end = snk_start + ctx.mKmerSize;

  for (auto const probe_id : mActiveProbeIds) {
    auto const& probe = mVariants[probe_id];
    auto const var_start = probe.mGenomeStart0;
    auto const var_end = var_start + probe.mRef.length();

    bool const overlaps_source = source.mFound && var_start < src_end && var_end > src_start;
    bool const overlaps_sink = sink.mFound && var_start < snk_end && var_end > snk_start;

    if (overlaps_source || overlaps_sink) {
      auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
      record.mIsVariantInAnchor = true;
      LOG_DEBUG("KMER_PROBE {} comp={} kmer_size={} probe={} {}:{} VARIANT_IN_ANCHOR", ctx.mRegStr,
                ctx.mCompId, ctx.mKmerSize, probe_id, probe.mChrom, probe.mGenomeStart0)
    }
  }
}

// ============================================================================
// CheckPaths: verify whether each probe's ALT context substring appears in
// any enumerated haplotype path.
//
// Builds a minimal ALT context string (the variant allele with k-1 flanking
// bases on each side) and searches each haplotype sequence for it. This is
// the definitive test: if the ALT context is in a path, the variant signal
// survived the entire graph pipeline.
// ============================================================================
void ProbeTracker::CheckPaths(absl::Span<Path const> haplotypes, Context const& ctx) {
  if (mActiveProbeIds.empty()) return;

  for (auto const probe_id : mActiveProbeIds) {
    auto const& probe = mVariants[probe_id];
    if (probe.mGenomeStart0 < ctx.mRegionStart) continue;

    auto const var_offset = probe.mGenomeStart0 - ctx.mRegionStart;
    auto const ref_allele_len = probe.mRef.length();
    auto const alt_context =
        BuildAltContext(ctx.mRefSeq, var_offset, ref_allele_len, probe.mAlt, ctx.mKmerSize);

    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
    for (usize hap_idx = 0; hap_idx < haplotypes.size(); ++hap_idx) {
      if (haplotypes[hap_idx].Sequence().find(alt_context) != std::string_view::npos) {
        record.mHapIndices.push_back(static_cast<u8>(hap_idx));
      }
    }

    LOG_DEBUG("KMER_PROBE {} comp={} kmer_size={} probe={} {}:{}:{}>{} haps=[{}] num_paths={}",
              ctx.mRegStr, ctx.mCompId, ctx.mKmerSize, probe_id, probe.mChrom, probe.mGenomeStart0,
              probe.mRef, probe.mAlt, fmt::join(record.mHapIndices, ","), haplotypes.size())
  }
}

void ProbeTracker::SetCycleRetry(Context const& ctx) {
  for (auto const probe_id : mActiveProbeIds) {
    FindOrCreateRecord(probe_id, ctx.mKmerSize).mIsCycleRetry = true;
  }
}

void ProbeTracker::SetComplexRetry(Context const& ctx) {
  for (auto const probe_id : mActiveProbeIds) {
    FindOrCreateRecord(probe_id, ctx.mKmerSize).mIsComplexRetry = true;
  }
}

void ProbeTracker::SetNoAnchor(Context const& ctx) {
  for (auto const probe_id : mActiveProbeIds) {
    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
    record.mIsNoAnchor = true;
    LOG_DEBUG("KMER_PROBE {} comp={} kmer_size={} probe={} {}:{} NO_ANCHOR", ctx.mRegStr,
              ctx.mCompId, ctx.mKmerSize, probe_id, mVariants[probe_id].mChrom,
              mVariants[probe_id].mGenomeStart0)
  }
}

void ProbeTracker::SetShortAnchor(Context const& ctx) {
  for (auto const probe_id : mActiveProbeIds) {
    auto& record = FindOrCreateRecord(probe_id, ctx.mKmerSize);
    record.mIsShortAnchor = true;
    LOG_DEBUG("KMER_PROBE {} comp={} kmer_size={} probe={} {}:{} SHORT_ANCHOR", ctx.mRegStr,
              ctx.mCompId, ctx.mKmerSize, probe_id, mVariants[probe_id].mChrom,
              mVariants[probe_id].mGenomeStart0)
  }
}

void ProbeTracker::SetTraversalLimit(Context const& ctx) {
  for (auto const probe_id : mActiveProbeIds) {
    FindOrCreateRecord(probe_id, ctx.mKmerSize).mIsTraversalLimited = true;
  }
}

auto ProbeTracker::GetHighlightNodeIds(NodeTable const& nodes, usize comp_id) const
    -> absl::flat_hash_set<NodeID> {
  absl::flat_hash_set<NodeID> result;
  if (!IsActive()) return result;

  for (auto const& [node_id, hits] : mNodeTags) {
    auto const node_it = nodes.find(node_id);
    if (node_it != nodes.end() && node_it->second->GetComponentId() == comp_id) {
      result.insert(node_id);
    }
  }
  return result;
}

auto ProbeTracker::FindOrCreateRecord(u16 probe_id, usize kmer_size) -> ProbeKRecord& {
  for (auto& record : mRecords) {
    if (record.mProbeId == probe_id && record.mKmerSize == kmer_size) return record;
  }

  mRecords.push_back(ProbeKRecord{.mKmerSize = kmer_size, .mProbeId = probe_id});
  return mRecords.back();
}

auto ProbeTracker::CountSurvivingKmers() const -> absl::flat_hash_map<u16, SurvivalProfile> {
  // Collect surviving offsets per probe
  absl::flat_hash_map<u16, std::vector<u16>> offsets_per_probe;
  for (auto const& [node_id, hits] : mNodeTags) {
    for (auto const& hit : hits) {
      offsets_per_probe[hit.mProbeId].push_back(hit.mAltContextOffset);
    }
  }

  absl::flat_hash_map<u16, SurvivalProfile> profiles;
  for (auto& [probe_id, offsets] : offsets_per_probe) {
    std::ranges::sort(offsets);
    auto const count = static_cast<u16>(offsets.size());

    // Find expected max offset from mAltKmerCounts
    auto const total_it = mAltKmerCounts.find(probe_id);
    auto const max_offset = (total_it != mAltKmerCounts.end() && total_it->second > 0)
                                ? static_cast<u16>(total_it->second - 1)
                                : u16{0};

    u16 edge_count = 0;
    u16 interior_count = 0;
    for (auto const offset : offsets) {
      if (offset == 0 || offset == max_offset) {
        ++edge_count;
      } else {
        ++interior_count;
      }
    }

    // Count contiguous gaps: positions where consecutive offsets differ by > 1
    u16 gap_count = 0;
    for (usize idx = 1; idx < offsets.size(); ++idx) {
      if (offsets[idx] - offsets[idx - 1] > 1) ++gap_count;
    }

    profiles[probe_id] = SurvivalProfile{
        .mCount = count,
        .mEdgeCount = edge_count,
        .mInteriorCount = interior_count,
        .mGapCount = gap_count,
    };
  }
  return profiles;
}

}  // namespace lancet::cbdg
