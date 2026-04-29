#ifndef SRC_LANCET_CORE_TAR_GZ_SHARD_MERGER_H_
#define SRC_LANCET_CORE_TAR_GZ_SHARD_MERGER_H_

#include <filesystem>
#include <utility>

namespace lancet::core {

// ============================================================================
// TarGzShardMerger: Bytewise concatenation of per-worker tar.gz shards into
// a single final archive, plus end-of-archive marker append and shards-dir cleanup.
//
// Workers each write their own gzipped TAR shard during the compute phase
// (each shard omits the trailing TAR end-of-archive zero blocks per RFC
// 1952's multi-member-gzip rule). After all workers join, this class reads
// every `worker_<idx>.tar.gz` shard under `mShardsDir` in deterministic
// index order and copies their bytes verbatim into `mFinalArchivePath`,
// appends a single gzipped end-of-archive marker once at the end, and
// removes the shards directory.
//
// Throws on I/O failure with the partial state preserved (shards intact)
// so the caller can log + exit non-zero without losing data.
// ============================================================================
class TarGzShardMerger {
 public:
  TarGzShardMerger(std::filesystem::path shards_dir, std::filesystem::path final_archive_path)
      : mShardsDir(std::move(shards_dir)), mFinalArchivePath(std::move(final_archive_path)) {}

  TarGzShardMerger(TarGzShardMerger const&) = delete;
  auto operator=(TarGzShardMerger const&) -> TarGzShardMerger& = delete;
  TarGzShardMerger(TarGzShardMerger&&) noexcept = delete;
  auto operator=(TarGzShardMerger&&) noexcept -> TarGzShardMerger& = delete;
  ~TarGzShardMerger() = default;

  /// Concatenate every `worker_<idx>.tar.gz` shard into the final archive,
  /// append the gzipped end-of-archive marker, remove the shards
  /// directory. Logs progress every 5 seconds via `EtaTimer` in the
  /// `Progress: X% | Elapsed: ... | ETA: ... @ ... shards/s` format.
  /// Throws std::runtime_error on I/O failure (shards preserved).
  void Merge();

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::filesystem::path mShardsDir;
  std::filesystem::path mFinalArchivePath;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_TAR_GZ_SHARD_MERGER_H_
