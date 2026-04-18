#include "lancet/base/crash_handler.h"
#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"

#include "absl/cleanup/cleanup.h"
#include "mimalloc.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <memory>

namespace {

// ============================================================================
// MimallocConfigurator — Tuned mimalloc defaults for Lancet2
//
// This struct's constructor executes during C++ static initialization, before
// main() starts. It configures mimalloc via mi_option_set() to embed the
// performance profile directly in the binary, so every user gets the tuned
// allocator without needing environment variables.
//
// Timing: mi_option_set() is safe to call before mi_process_init(). Values
// are stored internally and respected when mimalloc fully initializes. All
// options below are checked dynamically via atomics, so even if a pre-main
// arena is created before this constructor runs, all subsequent arenas (the
// hundreds of GiBs requested by worker threads) will use these settings.
//
// See: mimalloc_tuning_analysis.md for profiling data and rationale.
// ============================================================================
struct MimallocConfigurator {
  MimallocConfigurator() noexcept {
    // Disable memory purging — eliminates madvise(MADV_DONTNEED) syscall
    // overhead that accounts for ~8.3% CPU in profiling.
    // mi_option_set(mi_option_purge_delay, -1);

    // Fallback alternative to purge_delay=-1: use MADV_FREE instead of
    // MADV_DONTNEED. ~10x cheaper — kernel marks pages reclaimable but
    // doesn't unmap them. Use if RSS control is needed on shared nodes.
    mi_option_set(mi_option_purge_decommits, 0);

    // Enable 2MB large OS pages — reduces TLB pressure for large hash table
    // backing arrays (~512KB-1MB per NodeTable). Requires kernel THP support.
    mi_option_set(mi_option_allow_large_os_pages, 1);

    // Increase arena reservation from 1 GiB to 4 GiB — reduces mmap syscalls
    // and improves NUMA locality. API expects KiB: 4 GiB = 4L * 1024 * 1024.
    mi_option_set(mi_option_arena_reserve, 4L * 1024 * 1024);

    // Diagnostic: print mimalloc statistics on process exit.
    // Useful for verifying purge count drops to 0 and tracking peak RSS.
    mi_option_set(mi_option_show_stats, 1);

    // Diagnostic: print active mimalloc configuration on startup.
    // mi_option_set(mi_option_verbose, 1);
  }
};

MimallocConfigurator const MIMALLOC_CONFIG;

}  // namespace

auto main(int const argc, char const** argv) -> int {
  lancet::base::InstallCrashHandler();
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  lancet::RegisterLancetLogger();
  absl::Cleanup const cleanup = []() -> void {
    spdlog::shutdown();
    mi_collect(true);
  };

  lancet::cli::CliInterface cli;
  return cli.RunMain(argc, argv);
}
