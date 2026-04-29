#ifndef SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_
#define SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::base {
class TarGzWriter;  // forward decl — full def in lancet/base/tar_gz_writer.h
}  // namespace lancet::base

namespace lancet::cbdg {

// ============================================================================
// DotSnapshotBuffer — defer DOT writes until the successful k-attempt.
//
// `Graph::BuildComponentResults` retries assembly across k-values when the
// pruned graph has cycles or excessive complexity. Per-k-attempt buffering
// is what preserves the abandoned-k cleanup invariant: Discard() drops
// everything in memory on retry; Commit() appends each buffered DOT as a
// regular-file entry into the worker's per-thread `base::TarGzWriter`
// shard. Failed k-attempts never reach the shard, so they cannot leak DOT
// content into the final merged archive.
//
// Per-snapshot memory is small (~10–100 KB for post-compression graphs).
// Typical worst case per window worker holds a few MB in flight.
// ============================================================================
class DotSnapshotBuffer {
 public:
  DotSnapshotBuffer() = default;

  /// Set the per-window subdirectory name (e.g. "chr1_38506673_38507173")
  /// that Commit will use as the per-window directory under the stream
  /// root. Set once per window before any Buffer() call.
  void SetWindowSubdir(std::string subdir) noexcept { mWindowSubdir = std::move(subdir); }

  /// Buffer one rendered DOT file in memory. `filename` is the basename
  /// (no directory component); `dot_contents` is the full file body.
  void Buffer(std::string filename, std::string dot_contents) {
    mPending.emplace_back(std::move(filename), std::move(dot_contents));
  }

  /// Drop all buffered snapshots without touching disk. Called when the
  /// current k-attempt is abandoned (cycle / complexity retry).
  void Discard() noexcept { mPending.clear(); }

  /// Append all buffered snapshots into `shard_writer` as regular-file
  /// TAR entries. Each entry's archive path is
  /// `<stream_root>/<window_subdir>/<filename>`. Empties the buffer on
  /// return. No-op when the buffer is empty.
  void Commit(base::TarGzWriter& shard_writer, std::string_view stream_root);

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mPending.empty(); }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<std::pair<std::string, std::string>> mPending;
  std::string mWindowSubdir;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_
