#include "benchmark/benchmark.h"
#ifndef __APPLE__
#include "mimalloc-override.h"  // IWYU pragma: keep
#endif

// NOLINTNEXTLINE(misc-include-cleaner, cppcoreguidelines-pro-type-vararg, cppcoreguidelines-owning-memory, cert-err58-cpp)
BENCHMARK_MAIN();
