#include "lancet/core/active_region_detector.h"

#include "lancet/base/types.h"
#include "lancet/core/read_collector.h"
#include "lancet/core/sample_header_reader.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/iterator.h"
#include "lancet/hts/reference.h"
#include "lancet/hts/sam_flag.h"

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "absl/types/span.h"

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <cstdlib>

using CountMap = absl::flat_hash_map<u32, u32>;

namespace lancet::core {

auto HasMdTag(std::filesystem::path const& aln_path, std::filesystem::path const& ref_path)
    -> bool {
  static constexpr usize NUM_READS_TO_PEEK = 1000;
  static std::vector<std::string> const TAGS{"MD"};

  using hts::Alignment::Fields::AUX_RGAUX;
  hts::Reference const ref(ref_path);
  usize peeked_read_count = 0;
  hts::Extractor extractor(aln_path, ref, AUX_RGAUX, TAGS, true);

  for (auto const& aln : extractor) {
    if (peeked_read_count > NUM_READS_TO_PEEK) break;
    if (aln.HasTag("MD")) return true;
    peeked_read_count++;
  }
  return false;
}

}  // namespace lancet::core

// ---------------------------------------------------------------------------
// Anonymous namespace: helpers used exclusively by IsActiveRegion below.
// ---------------------------------------------------------------------------

namespace {

// ============================================================================
// ParseMd — extract mismatch positions from the MD:Z auxiliary tag
//
// PURPOSE: Detects positions where ≥2 reads disagree with the reference,
// which signals an active region worth assembling. This is a lightweight
// pre-filter before the full de Bruijn graph assembly.
//
// MD:Z FORMAT (SAM spec §1.5):
//   Concatenation of: [0-9]+ (matching run) | [ACGT] (mismatch base) | ^[ACGT]+ (deletion)
//   Example: "10A5^AC6"  →  10 matches, A→? mismatch, 5 matches, 2bp deletion, 6 matches
//
// STATE MACHINE:
//   ┌──────────┐  digit   ┌──────────────┐
//   │ scanning ├─────────►│ accumulating │
//   │          │◄─────────┤ match-run len│
//   └────┬─────┘  letter  └──────────────┘
//        │ (mismatch base A/C/G/T)
//        ▼
//   record genome_pos in result map; if count reaches 2 → return true
//
// Returns true as soon as any genome position has ≥2 mismatches (early exit).
// ============================================================================
inline auto ParseMd(std::string_view md_val, absl::Span<u8 const> quals, i64 const start,
                    CountMap* result) -> bool {
  if (start < 0) return false;

  std::string token;
  token.reserve(md_val.length());
  auto genome_pos = static_cast<u32>(start);

  for (auto const& character : md_val) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(character))) {
      token += character;
      continue;
    }

    auto const step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genome_pos += static_cast<u32>(step);
    token.clear();

    auto const base_pos = static_cast<usize>(genome_pos - start);
    static constexpr u8 MIN_BASE_QUAL = 20;
    if (quals.at(base_pos) < MIN_BASE_QUAL) continue;

    auto const base = absl::ascii_toupper(static_cast<unsigned char>(character));
    if (base == 'A' || base == 'C' || base == 'T' || base == 'G') {
      if (++(*result)[genome_pos] == 2) return true;
    }
  }

  return false;
}

// ---------------------------------------------------------------------------
// IncrementHitCount — returns true when a genome position reaches ≥2 hits.
// Active region detection uses a threshold of 2: a single read with a
// mismatch/indel is noise; two reads at the same position is signal.
// ---------------------------------------------------------------------------
inline auto IncrementHitCount(CountMap& counts, u32 const genome_pos) -> bool {
  return ++counts[genome_pos] == 2;
}

// ============================================================================
// MutationAccumulator — per-sample evidence tracker for active region detection
//
// Owns four CountMaps (mismatches, insertions, deletions, softclips) keyed
// by genome position. CheckAlignment runs three independent checks per read:
//   1. MD tag mismatches (via ParseMd)
//   2. CIGAR-based indels and mismatches
//   3. Soft-clip positions
// Returns true the moment any position accumulates ≥2 supporting reads.
//
// Lifetime: one instance per IsActiveRegion call, cleared between samples.
// ============================================================================
class MutationAccumulator {
 public:
  CountMap mMismatches;
  CountMap mInsertions;
  CountMap mDeletions;
  CountMap mSoftclips;
  std::vector<u32> mSoftclipPositions;

  void ClearAll() {
    mMismatches.clear();
    mInsertions.clear();
    mDeletions.clear();
    mSoftclips.clear();
    mSoftclipPositions.clear();
  }

  /// Entry point: filter QC-fail/dup/unmapped/mapq0 reads, then check
  /// MD mismatches → CIGAR events → soft-clips in order of cost.
  [[nodiscard]] auto CheckAlignment(lancet::hts::Alignment const& aln) -> bool {
    auto const bflag = aln.Flag();
    if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() == 0) {
      return false;
    }

    if (CheckMdTag(aln)) return true;
    if (CheckCigarEvents(aln)) return true;
    return CheckSoftClips(aln);
  }

 private:
  /// Check MD:Z tag for reference mismatches. BuildQualities() is called
  /// on-demand only when the MD tag is present (avoids deep copy otherwise).
  [[nodiscard]] auto CheckMdTag(lancet::hts::Alignment const& aln) -> bool {
    if (!aln.HasTag("MD")) return false;
    auto const md_tag = aln.GetTag<std::string_view>("MD");
    // BuildQualities performs on-demand deep copy of quality values
    auto const quals = aln.BuildQualities();
    return ParseMd(md_tag.value(), absl::MakeConstSpan(quals), aln.StartPos0(), &mMismatches);
  }

  /// Scan CIGAR operations for insertions, deletions, and
  /// explicit mismatch ops (X). Skips alignment-consuming ops.
  [[nodiscard]] auto CheckCigarEvents(lancet::hts::Alignment const& aln) -> bool {
    auto const cigar_units = aln.CigarData();
    auto curr_genome_pos = static_cast<u32>(aln.StartPos0());

    for (auto const& cig_unit : cigar_units) {
      if (cig_unit.ConsumesReference()) {
        curr_genome_pos += cig_unit.Length();
      }

      switch (cig_unit.Operation()) {
        case lancet::hts::CigarOp::INSERTION:
          if (IncrementHitCount(mInsertions, curr_genome_pos)) return true;
          break;
        case lancet::hts::CigarOp::DELETION:
          if (IncrementHitCount(mDeletions, curr_genome_pos)) return true;
          break;
        case lancet::hts::CigarOp::SEQUENCE_MISMATCH:
          if (IncrementHitCount(mMismatches, curr_genome_pos)) return true;
          break;
        default:
          break;
      }
    }
    return false;
  }

  /// Extract soft-clip genome positions and check for ≥2 reads clipped
  /// at the same position — a common signal for structural variant edges.
  [[nodiscard]] auto CheckSoftClips(lancet::hts::Alignment const& aln) -> bool {
    mSoftclipPositions.clear();
    if (!aln.GetSoftClips(nullptr, nullptr, &mSoftclipPositions, false)) return false;
    return std::ranges::any_of(mSoftclipPositions,
                               [this](u32 gpos) { return IncrementHitCount(mSoftclips, gpos); });
  }
};

}  // namespace

namespace lancet::core {

auto IsActiveRegion(ReadCollector::Params const& params, hts::Reference::Region const& region)
    -> bool {
  MutationAccumulator accumulator;

  auto const sample_list = core::MakeSampleList(params);
  for (auto const& sinfo : sample_list) {
    accumulator.ClearAll();

    using hts::Alignment::Fields::AUX_RGAUX;
    hts::Extractor extractor(sinfo.Path(), hts::Reference(params.mRefPath), AUX_RGAUX, {"MD"},
                             params.mNoCtgCheck);
    extractor.SetRegionToExtract(region.ToSamtoolsRegion());

    for (auto const& aln : extractor) {
      if (accumulator.CheckAlignment(aln)) return true;
    }
  }

  return false;
}

}  // namespace lancet::core
