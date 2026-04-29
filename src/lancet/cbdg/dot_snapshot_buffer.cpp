#include "lancet/cbdg/dot_snapshot_buffer.h"

#include "lancet/base/tar_gz_writer.h"

#include <filesystem>
#include <string_view>

namespace lancet::cbdg {

void DotSnapshotBuffer::Commit(base::TarGzWriter& shard_writer, std::string_view stream_root) {
  if (mPending.empty()) return;

  // Build entry paths as `<stream_root>/<window_subdir>/<filename>` so
  // the merged archive's extracted layout has one directory per window
  // (e.g. `dbg_graph/chr1_X_Y/<file>.dot`).
  auto const window_subdir_path = std::filesystem::path(stream_root) / mWindowSubdir;

  for (auto const& [entry_filename, entry_contents] : mPending) {
    auto const entry_path = window_subdir_path / entry_filename;
    shard_writer.AddRegularFileEntry(entry_path.string(), entry_contents);
  }

  mPending.clear();
}

}  // namespace lancet::cbdg
