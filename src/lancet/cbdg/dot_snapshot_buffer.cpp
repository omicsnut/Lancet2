#include "lancet/cbdg/dot_snapshot_buffer.h"

#include <filesystem>
#include <fstream>

namespace lancet::cbdg {

void DotSnapshotBuffer::Commit(std::filesystem::path const& out_dir) {
  if (mPending.empty()) return;

  std::filesystem::create_directories(out_dir);

  for (auto& [filename, contents] : mPending) {
    std::ofstream out_handle(out_dir / filename, std::ios::trunc);
    out_handle << contents;
  }

  mPending.clear();
}

}  // namespace lancet::cbdg
