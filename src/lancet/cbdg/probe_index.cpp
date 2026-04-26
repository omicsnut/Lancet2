#include "lancet/cbdg/probe_index.h"

#include "lancet/base/logging.h"
#include "lancet/base/sliding.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
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
//
// REF context length: (k-1) + R + (k-1) = 2k + R - 2
// ALT context length: (k-1) + A + (k-1) = 2k + A - 2  (A = alt_allele_len)
//
// Boundary clipping: ctx_start is clamped to 0, ctx_end to ref_seq.length().
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
// ============================================================================
[[nodiscard]] auto BuildAltContext(std::string_view ref_seq, usize var_offset, usize ref_allele_len,
                                   std::string_view alt_allele, usize kmer_size) -> std::string {
  auto const [ctx_start, ctx_end] =
      ComputeFlankRange(ref_seq.length(), var_offset, ref_allele_len, kmer_size);

  std::string alt_context;
  alt_context.reserve(ctx_end - ctx_start + alt_allele.length());
  alt_context.append(ref_seq.substr(ctx_start, var_offset - ctx_start));
  alt_context.append(alt_allele);
  if (var_offset + ref_allele_len < ctx_end) {
    alt_context.append(
        ref_seq.substr(var_offset + ref_allele_len, ctx_end - var_offset - ref_allele_len));
  }
  return alt_context;
}

// ============================================================================
// CollectRefKmerHashes: hash all canonical k-mers from the reference flanking
// context around a variant. Used to identify ALT-unique k-mers by set
// difference with the ALT context k-mers.
// ============================================================================
[[nodiscard]] auto CollectRefKmerHashes(std::string_view ref_seq, usize var_offset,
                                        usize ref_allele_len, usize kmer_size)
    -> absl::flat_hash_set<u64> {
  auto const [ctx_start, ctx_end] =
      ComputeFlankRange(ref_seq.length(), var_offset, ref_allele_len, kmer_size);
  auto const ref_context = ref_seq.substr(ctx_start, ctx_end - ctx_start);

  absl::flat_hash_set<u64> ref_hashes;
  if (ref_context.length() >= kmer_size) {
    auto const sliding = lancet::base::SlidingView(ref_context, kmer_size);
    ref_hashes.reserve(sliding.size());
    for (auto const& kmer_seq : sliding) {
      ref_hashes.insert(Kmer(kmer_seq).Identifier());
    }
  }
  return ref_hashes;
}

// ============================================================================
// IndexProbeAtK: fetch flanking reference context for a single probe at a
// given k-mer size, splice in the ALT allele, and insert ALT-unique k-mers
// into the output map.
// ============================================================================
void IndexProbeAtK(ProbeVariant const& probe, hts::Reference const& reference, usize kmer_size,
                   ProbeIndex::KmerToProbes& kmer_map) {
  auto const flank_len = kmer_size - 1;
  auto const ref_allele_len = probe.mRef.length();

  // probe.mGenomeStart0 is 0-based, so start1 = mGenomeStart0 + 1.
  auto const var_start1 = probe.mGenomeStart0 + 1;
  auto const fetch_start1 = var_start1 > flank_len ? var_start1 - flank_len : u64{1};
  auto const fetch_end1 = var_start1 + ref_allele_len - 1 + flank_len;

  auto const chrom_result = reference.FindChromByName(probe.mChrom);
  if (!chrom_result.ok()) return;
  auto const chrom_len = chrom_result->Length();
  auto const clamped_end1 = std::min<u64>(fetch_end1, chrom_len);
  if (fetch_start1 > clamped_end1) return;

  hts::Reference::OneBasedClosedOptional const interval = {fetch_start1, clamped_end1};
  auto const region = reference.MakeRegion(probe.mChrom, interval);
  auto const ref_seq = region.SeqView();

  auto const var_offset = static_cast<usize>(var_start1 - fetch_start1);
  auto const ref_hashes = CollectRefKmerHashes(ref_seq, var_offset, ref_allele_len, kmer_size);
  auto const alt_context =
      BuildAltContext(ref_seq, var_offset, ref_allele_len, probe.mAlt, kmer_size);

  // Insert ALT-unique k-mers (present in ALT context, absent from REF context).
  if (alt_context.length() < kmer_size) return;

  u16 alt_kmer_offset = 0;
  for (auto const& kmer_seq : lancet::base::SlidingView(alt_context, kmer_size)) {
    auto const kmer_hash = Kmer(kmer_seq).Identifier();
    if (ref_hashes.contains(kmer_hash)) continue;

    kmer_map[kmer_hash].emplace_back(
        ProbeIndex::KmerEntry{.mProbeId = probe.mProbeId, .mAltContextOffset = alt_kmer_offset});
    ++alt_kmer_offset;
  }
}

}  // namespace

auto ProbeIndex::LoadVariantsFromFile(std::filesystem::path const& path)
    -> std::vector<ProbeVariant> {
  std::vector<ProbeVariant> variants;
  std::ifstream infile(path);
  if (!infile.is_open()) {
    LOG_WARN("ProbeIndex could not open probe variants file: {}", path.string())
    return variants;
  }

  std::string line;
  // Skip header line
  std::getline(infile, line);

  u16 probe_idx = 0;
  usize file_line = 1;  // header is line 1
  while (std::getline(infile, line)) {
    ++file_line;
    // missed_variants.txt columns (tab-separated):
    // chrom, start, end, ref, alt, type, length, classification,
    // tier1_alt_count, tier2_alt_count, max_alt_count, engine, source
    std::vector<absl::string_view> const fields = absl::StrSplit(line, '\t');
    static constexpr usize MIN_FIELDS = 9;
    if (fields.size() < MIN_FIELDS) continue;

    static constexpr usize CHROM_IDX = 0;
    static constexpr usize START_IDX = 1;
    static constexpr usize REF_IDX = 3;
    static constexpr usize ALT_IDX = 4;
    static constexpr usize TIER1_IDX = 8;

    u16 tier1_count = 0;
    if (!absl::SimpleAtoi(fields[TIER1_IDX], &tier1_count)) {
      LOG_WARN("ProbeIndex could not parse tier1_alt_count '{}' at line {} from file {}",
               fields[TIER1_IDX], file_line, path.string());
      continue;
    }

    usize genome_start = 0;
    if (!absl::SimpleAtoi(fields[START_IDX], &genome_start)) {
      LOG_WARN("ProbeIndex could not parse start position '{}' at line {} from file {}",
               fields[START_IDX], file_line, path.string());
      continue;
    }

    variants.emplace_back(ProbeVariant{
        .mChrom = std::string(fields[CHROM_IDX]),
        .mRef = std::string(fields[REF_IDX]),
        .mAlt = std::string(fields[ALT_IDX]),
        .mGenomeStart0 = genome_start,
        .mProbeId = probe_idx,
        .mTier1AltCount = tier1_count,
    });

    ++probe_idx;
  }

  LOG_INFO("ProbeIndex loaded {} probe variants from {}", variants.size(), path.string())
  return variants;
}

auto ProbeIndex::ForKmerSize(usize const kmer_size) const -> KmerToProbes const* {
  auto const iter = mPerKIndex.find(kmer_size);
  return iter != mPerKIndex.end() ? &iter->second : nullptr;
}

// ============================================================================
// ProbeIndex::Build — precompute ALT-unique canonical k-mer hashes for all
// probe variants across every k in [min_kmer_len..max_kmer_len, step].
//
// For each (variant, k) pair:
//   1. Fetch the variant's flanking reference context from the FASTA.
//   2. Splice the ALT allele into the reference context.
//   3. Hash all canonical k-mers from both REF and ALT contexts.
//   4. K-mers present in ALT but absent from REF are "ALT-unique."
//   5. Store each ALT-unique hash → {probe_id, offset} in mPerKIndex[k].
//
// The reference context requires fetching a region of length:
//   (k-1) + max(R, A) + (k-1) bases centered on the variant.
// We fetch the full flanking region once per probe per chromosome and
// compute the local offset from the fetched region.
// ============================================================================
auto ProbeIndex::Build(absl::Span<ProbeVariant const> variants, hts::Reference const& reference,
                       usize const min_kmer_len, usize const max_kmer_len,
                       usize const kmer_step_len) -> ProbeIndex {
  ProbeIndex index;
  if (variants.empty()) return index;

  for (usize kmer_size = min_kmer_len; kmer_size <= max_kmer_len; kmer_size += kmer_step_len) {
    auto& kmer_map = index.mPerKIndex[kmer_size];
    for (auto const& probe : variants) {
      IndexProbeAtK(probe, reference, kmer_size, kmer_map);
    }
  }

  return index;
}

auto ProbeIndex::IsEmpty() const noexcept -> bool {
  return mPerKIndex.empty();
}

auto ProbeIndex::TotalEntries() const noexcept -> usize {
  usize total = 0;
  for (auto const& [_kmer_size, kmer_map] : mPerKIndex) {
    for (auto const& [_hash, entries] : kmer_map) {
      total += entries.size();
    }
  }
  return total;
}

}  // namespace lancet::cbdg
