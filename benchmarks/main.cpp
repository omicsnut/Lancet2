#include "benchmark/benchmark.h"
#if !defined(__APPLE__) && !defined(LANCET_SANITIZE_BUILD)
#include "mimalloc-override.h"  // IWYU pragma: keep
#endif

BENCHMARK_MAIN();
