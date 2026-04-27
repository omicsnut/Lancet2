#ifndef SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_
#define SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace lancet::cbdg {

// ============================================================================
// DotSnapshotBuffer — defer DOT writes until the successful k-attempt.
//
// Graph::BuildComponentResults retries assembly across k-values when the
// pruned graph has cycles or excessive complexity. Eager per-k DOT writes
// leave abandoned-k artifacts on disk. This buffer accumulates rendered DOT
// content in memory; the outer loop calls Discard() on retry and Commit()
// once on success, so only the final attempt's snapshots reach disk.
//
// Per-snapshot memory is small (~10–100 KB for post-compression graphs).
// Typical worst case per window worker holds a few MB in flight.
// ============================================================================
class DotSnapshotBuffer {
 public:
  DotSnapshotBuffer() = default;

  /// Buffer one rendered DOT file in memory. `filename` is the basename
  /// (no directory component); `dot_contents` is the full file body.
  void Buffer(std::string filename, std::string dot_contents) {
    mPending.emplace_back(std::move(filename), std::move(dot_contents));
  }

  /// Drop all buffered snapshots without touching disk. Called when the
  /// current k-attempt is abandoned (cycle / complexity retry).
  void Discard() noexcept { mPending.clear(); }

  /// Write all buffered snapshots to `out_dir` (creates the directory if
  /// missing). Empties the buffer on return. No-op when the buffer is empty.
  void Commit(std::filesystem::path const& out_dir);

  [[nodiscard]] auto IsEmpty() const noexcept -> bool { return mPending.empty(); }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::vector<std::pair<std::string, std::string>> mPending;
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_DOT_SNAPSHOT_BUFFER_H_
