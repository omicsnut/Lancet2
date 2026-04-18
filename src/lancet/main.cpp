#include "lancet/base/crash_handler.h"
#include "lancet/base/logging.h"
#include "lancet/cli/cli_interface.h"

#include "absl/cleanup/cleanup.h"
#include "mimalloc.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <memory>

auto main(int const argc, char const** argv) -> int {
  // Disable mimalloc's arena purge cycle entirely. Without this, mimalloc
  // calls madvise(MADV_DONTNEED) to decommit freed pages every ~100ms,
  // costing ~200s of CPU time over a full run. Since graph node count is
  // bounded, per-thread heap memory plateaus early and is reused across
  // windows — purging just adds decommit/recommit churn with no benefit.
  // mi_collect(true) at exit releases everything.
  mi_option_set(mi_option_purge_delay, -1);

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
