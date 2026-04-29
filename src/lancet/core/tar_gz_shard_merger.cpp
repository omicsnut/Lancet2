#include "lancet/core/tar_gz_shard_merger.h"

#include "lancet/base/eta_timer.h"
#include "lancet/base/gzip_ostream.h"
#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"

#include "absl/time/time.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace lancet::core {

// 64 KiB buffer for the bytewise shard-to-final copy loop. Big enough that
// each `read`/`write` call amortises syscall overhead, small enough to
// stay in L2 cache on every CPU.
static constexpr usize SHARD_COPY_BUFFER_BYTES = 65'536;

// Throttle progress logging during the copy loop. Workers may have produced
// hundreds of shards so per-shard logging is noisy; 5 seconds between
// updates matches the cadence of the per-window progress logs without
// flooding the log file.
static constexpr auto MERGE_PROGRESS_LOG_INTERVAL = std::chrono::seconds(5);

namespace {

// Helper: enumerate `worker_<idx>.tar.gz` files under `shards_dir`,
// parse the index, and return the paths sorted by index ascending.
[[nodiscard]] auto CollectShardPathsInOrder(std::filesystem::path const& shards_dir)
    -> std::vector<std::filesystem::path> {
  static std::regex const SHARD_FILENAME_PATTERN(R"(^worker_(\d+)\.tar\.gz$)");

  std::vector<std::pair<u32, std::filesystem::path>> indexed_shards;
  for (auto const& entry : std::filesystem::directory_iterator(shards_dir)) {
    if (!entry.is_regular_file()) continue;

    auto const filename_str = entry.path().filename().string();
    std::smatch parsed_match;
    if (!std::regex_match(filename_str, parsed_match, SHARD_FILENAME_PATTERN)) continue;

    auto const worker_index = static_cast<u32>(std::stoul(parsed_match[1].str()));
    indexed_shards.emplace_back(worker_index, entry.path());
  }

  std::ranges::sort(indexed_shards,
                    [](auto const& lhs, auto const& rhs) { return lhs.first < rhs.first; });

  std::vector<std::filesystem::path> shard_paths_in_order;
  shard_paths_in_order.reserve(indexed_shards.size());
  for (auto& [worker_index, shard_path] : indexed_shards) {
    shard_paths_in_order.push_back(std::move(shard_path));
  }
  return shard_paths_in_order;
}

// Copy bytes from `source_path` into the open `output_stream` via a
// fixed-size scratch buffer. Throws on read or write failure.
void CopyShardBytesIntoOutput(std::filesystem::path const& source_path,
                              std::ofstream& output_stream) {
  std::ifstream input_stream(source_path, std::ios::binary);
  if (!input_stream.is_open()) {
    throw std::runtime_error(
        fmt::format("TarGzShardMerger: failed to open shard for read: {}", source_path.string()));
  }

  std::array<char, SHARD_COPY_BUFFER_BYTES> copy_scratch_buffer{};
  while (input_stream) {
    input_stream.read(copy_scratch_buffer.data(),
                      static_cast<std::streamsize>(copy_scratch_buffer.size()));
    auto const bytes_read = static_cast<usize>(input_stream.gcount());
    if (bytes_read == 0) break;

    output_stream.write(copy_scratch_buffer.data(), static_cast<std::streamsize>(bytes_read));
    if (!output_stream.good()) {
      throw std::runtime_error(fmt::format("TarGzShardMerger: write failure while merging shard {}",
                                           source_path.string()));
    }
  }

  if (input_stream.bad()) {
    throw std::runtime_error(
        fmt::format("TarGzShardMerger: read failure on shard {}", source_path.string()));
  }
}

// Append a final gzip-wrapped end-of-archive marker (1024 zero bytes,
// gzipped) to the open output stream. We do this by spinning up a
// temporary GzipOstream that writes to its own scratch file, then copying
// that file's bytes onto the end of the final archive. Slightly wasteful
// (one tiny temp file) but uses the same well-tested gzip path as
// everything else; the overhead is microscopic.
//
// Alternative: pre-compute the gzipped 1024-zero bytes and embed as a
// constexpr byte array. Skipped because the runtime cost is negligible
// and embedded bytes are brittle if zlib's framing ever changes.
void AppendEndOfArchiveMarker(std::ofstream& output_stream,
                              std::filesystem::path const& shards_dir) {
  auto const trailer_scratch_path = shards_dir / ".end_of_archive_trailer.tar.gz";
  static constexpr usize TAR_END_OF_ARCHIVE_BYTES = 1024;
  static constexpr std::array<char, TAR_END_OF_ARCHIVE_BYTES> TAR_END_ZEROS{};

  {
    base::GzipOstream trailer_gzip_stream(trailer_scratch_path);
    trailer_gzip_stream.Write(std::string_view(TAR_END_ZEROS.data(), TAR_END_ZEROS.size()));
  }  // Close() runs in destructor; gzip trailer is written.

  CopyShardBytesIntoOutput(trailer_scratch_path, output_stream);

  std::error_code remove_error;
  std::filesystem::remove(trailer_scratch_path, remove_error);
  if (remove_error) {
    LOG_WARN("TarGzShardMerger: failed to remove temporary trailer file {}: {}",
             trailer_scratch_path.string(), remove_error.message())
  }
}

}  // namespace

void TarGzShardMerger::Merge() {
  auto const shard_paths_in_order = CollectShardPathsInOrder(mShardsDir);
  if (shard_paths_in_order.empty()) {
    LOG_WARN("TarGzShardMerger: no worker shards found under {}; producing an empty archive",
             mShardsDir.string())
  }

  std::ofstream output_stream(mFinalArchivePath, std::ios::binary | std::ios::trunc);
  if (!output_stream.is_open()) {
    throw std::runtime_error(
        fmt::format("TarGzShardMerger: failed to open final archive for write: {}",
                    mFinalArchivePath.string()));
  }

  // EtaTimer treats one shard as one unit. With ~96 shards on chr1 and
  // ~250ms per shard at ~200 MB/s, the ETA stabilises after a handful of
  // shards which is plenty for the user-facing progress format.
  base::EtaTimer merge_eta_timer(shard_paths_in_order.size());
  base::Timer merge_wallclock_timer;
  base::Timer last_log_throttle_timer;
  static constexpr auto ELAPSED_PRECISION = absl::Seconds(1);
  bool first_iteration = true;

  for (usize shard_index = 0; shard_index < shard_paths_in_order.size(); ++shard_index) {
    auto const& shard_path = shard_paths_in_order[shard_index];
    CopyShardBytesIntoOutput(shard_path, output_stream);
    merge_eta_timer.Increment();

    auto const should_log_progress =
        first_iteration ||
        last_log_throttle_timer.Runtime() >= absl::FromChrono(MERGE_PROGRESS_LOG_INTERVAL);
    if (should_log_progress) {
      auto const completed_count = shard_index + 1;
      auto const percent_complete =
          100.0 * static_cast<f64>(completed_count) / static_cast<f64>(shard_paths_in_order.size());
      // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
      auto const elapsed_str =
          absl::FormatDuration(absl::Trunc(merge_wallclock_timer.Runtime(), ELAPSED_PRECISION));
      // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
      auto const eta_str =
          absl::FormatDuration(absl::Trunc(merge_eta_timer.EstimatedEta(), ELAPSED_PRECISION));
      LOG_INFO("Merging graph output shards: {:>8.4f}% | Elapsed: {} | ETA: {} @ {:.2f} shards/s | "
               "{} of {} shards merged",
               percent_complete, elapsed_str, eta_str, merge_eta_timer.RatePerSecond(),
               completed_count, shard_paths_in_order.size())
      last_log_throttle_timer.Reset();
      first_iteration = false;
    }
  }

  AppendEndOfArchiveMarker(output_stream, mShardsDir);

  output_stream.close();
  if (output_stream.fail() && !output_stream.eof()) {
    throw std::runtime_error(fmt::format("TarGzShardMerger: ofstream close reported failure for {}",
                                         mFinalArchivePath.string()));
  }

  // Successful merge: remove the shards directory. If removal fails we
  // log but don't propagate — the merged archive is good and orphan
  // shards are recoverable.
  std::error_code remove_error;
  std::filesystem::remove_all(mShardsDir, remove_error);
  if (remove_error) {
    LOG_WARN("TarGzShardMerger: merge succeeded but failed to remove shards dir {}: {}",
             mShardsDir.string(), remove_error.message())
  }
}

}  // namespace lancet::core
