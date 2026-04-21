#ifndef SRC_LANCET_CORE_READ_COLLECTOR_H_
#define SRC_LANCET_CORE_READ_COLLECTOR_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/mate_info.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lancet::core {

class ReadCollector {
 public:
  static constexpr f64 DEFAULT_MAX_WINDOW_COVERAGE = 1000.0;

  using ExtractorPtr = std::unique_ptr<hts::Extractor>;
  using SampleExtractors =
      absl::flat_hash_map<SampleInfo, ExtractorPtr, SampleInfo::Hash, SampleInfo::Equal>;

  struct Params {
    // ── 8B Align ────────────────────────────────────────────────────────────
    std::filesystem::path mRefPath;                 // 8B+
    std::vector<std::filesystem::path> mCtrlPaths;  // 8B+
    std::vector<std::filesystem::path> mCasePaths;  // 8B+
    /// Forward-facing unified sample input. Each entry is "<path>:<role>"
    /// where role is control or case. If role is omitted, defaults to control.
    std::vector<std::string> mSampleSpecs;            // 8B+
    f64 mMaxSampleCov = DEFAULT_MAX_WINDOW_COVERAGE;  // 8B
    // ── 1B Align ────────────────────────────────────────────────────────────
    bool mNoCtgCheck = false;    // 1B
    bool mExtractPairs = false;  // 1B

    /// Number of input file paths (NOT unique logical samples).
    /// Unique sample count is determined after sorting by MakeSampleList().
    [[nodiscard]] auto InputPathCount() const -> usize {
      return mCtrlPaths.size() + mCasePaths.size() + mSampleSpecs.size();
    }
  };

  using Read = cbdg::Read;
  using Region = hts::Reference::Region;

  ReadCollector() = delete;
  ReadCollector(Params params, absl::Span<SampleInfo const> sample_list);

  struct Result {
    // ── 8B Align ────────────────────────────────────────────────────────────
    std::vector<Read> mSampleReads;       // 8B+
    std::vector<SampleInfo> mSampleList;  // 8B+
  };

  [[nodiscard]] auto CollectRegionResult(Region const& region) -> Result;
  [[nodiscard]] auto IsCaseCtrlMode() const noexcept -> bool { return mIsCaseCtrlMode; }

  /// Expose cached sample list for IsActiveRegion (per-thread, not shared).
  [[nodiscard]] auto SampleList() const noexcept -> absl::Span<SampleInfo const> {
    return absl::MakeConstSpan(mSampleList);
  }

  /// Expose per-sample extractors for IsActiveRegion (per-thread, not shared).
  [[nodiscard]] auto Extractors() noexcept -> SampleExtractors& { return mExtractors; }

 private:
  using AlnAndRefPaths = std::array<std::filesystem::path, 2>;

  /// Maps qname hash (u64) -> mate location info for out-of-region mate retrieval.
  /// Using u64 hashes instead of std::string keys avoids string copies during Pass 1.
  using MateRegionsMap = absl::flat_hash_map<u64, hts::MateInfo>;
  using MateHashAndLocation = std::pair<u64, hts::MateInfo>;

  /// Pass 1 output: downsampling decisions + mate locations.
  struct ProfileResult {
    // ── 8B Align ────────────────────────────────────────────────────────────
    absl::flat_hash_set<u64> mKeepQnames;  // 8B+ — qname hashes kept after downsampling
    MateRegionsMap mExpectedMates;         // 8B+ — mate locations for out-of-region retrieval
    u64 mSampledReadCount = 0;             // 8B  — number of reads kept
  };

  // ── 8B Align ────────────────────────────────────────────────────────────
  Params mParams;                       // 8B+ — immutable construction params
  SampleExtractors mExtractors;         // 8B+ — per-sample HTSlib extractors
  std::vector<SampleInfo> mSampleList;  // 8B+ — sorted sample metadata
  std::vector<Read> mSampledReads;      // 8B+ — reusable per-region accumulator
  u64 mSampledBaseCount = 0;            // 8B  — tracks bases across Pass 2 + 3
  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mIsCaseCtrlMode{false};  // 1B

  /// Computes a deterministic u64 hash for a query name string_view.
  [[nodiscard]] static auto HashQname(std::string_view qname) -> u64;

  /// Pass 1: zero-copy profiling + deterministic downsampling.
  /// Scans all alignments in the region for a single sample, computing
  /// coverage statistics and selecting which reads survive downsampling.
  [[nodiscard]] auto ProfileAndDownsample(hts::Extractor& extractor, std::string const& region_spec,
                                          f64 max_sample_bases) const -> ProfileResult;

  /// Pass 2: re-iterate region, deep-copy only kept reads into mSampledReads.
  /// Each read gets tagged with sample metadata (name, index, label).
  void ExtractKeptReads(hts::Extractor& extractor, std::string const& region_spec,
                        absl::flat_hash_set<u64> const& keep_qnames, SampleInfo const& sinfo);

  /// Pass 3: fetch out-of-region mates for reads with distant mates.
  /// Walks mate locations in reverse-sorted genomic order for cache efficiency.
  void RecaptureMates(hts::Extractor& extractor, absl::flat_hash_set<u64> const& keep_qnames,
                      MateRegionsMap& expected_mates, SampleInfo const& sinfo);

  /// Sort mate regions in descending genomic order (largest chrom + position first).
  /// RecaptureMates processes via pop_back(), yielding ascending order for
  /// sequential BAM/CRAM disk access with O(1) removal per element.
  [[nodiscard]] static auto ReverseSortMateRegions(MateRegionsMap const& data)
      -> std::vector<MateHashAndLocation>;

  [[nodiscard]] static auto MakeRegionSpec(hts::MateInfo const& info, hts::Extractor const* ext)
      -> std::string;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_READ_COLLECTOR_H_
