#include "lancet/core/pipeline_executor.h"

#include <blockingconcurrentqueue.h>

#ifdef LANCET_PROFILE_MODE
#include "gperftools/profiler.h"
#endif

#include "lancet/base/eta_timer.h"
#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/core/async_worker.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/variant_store.h"
#include "lancet/core/window.h"
#include "lancet/core/window_builder.h"

#include "absl/container/fixed_array.h"
#include "absl/hash/hash.h"
#include "absl/time/time.h"
#include "concurrentqueue.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <memory>
#include <numeric>
#include <stop_token>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

/// Truncate a duration to the given precision, then format as a human-readable string.
/// Wraps the common absl::FormatDuration(absl::Trunc(...)) two-step into a single call.
[[nodiscard]] inline auto FormatTruncatedDuration(absl::Duration duration, absl::Duration precision)
    -> std::string {
  return absl::FormatDuration(absl::Trunc(duration, precision));
}

}  // namespace

namespace lancet::core {

using base::EtaTimer;

namespace {

[[nodiscard]] inline auto InitWindowStats() -> PipelineExecutor::WindowStats {
  using enum VariantBuilder::StatusCode;
  return PipelineExecutor::WindowStats{{UNKNOWN, 0},
                                       {SKIPPED_NONLY_REF_BASES, 0},
                                       {SKIPPED_REF_REPEAT_SEEN, 0},
                                       {SKIPPED_INACTIVE_REGION, 0},
                                       {SKIPPED_NOASM_HAPLOTYPE, 0},
                                       {MISSING_NO_MSA_VARIANTS, 0},
                                       {FOUND_GENOTYPED_VARIANT, 0}};
}

}  // namespace

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
PipelineExecutor::PipelineExecutor(WindowBuilder builder,
                                   std::shared_ptr<VariantBuilder::Params const> params,
                                   usize num_threads, u32 window_length)
    : mWindowBuilder(std::move(builder)),
      mParams(std::move(params)),
      mNumThreads(num_threads),
      mWindowLength(window_length) {}

// ---------------------------------------------------------------------------
// LogWindowStats — static method for final window status breakdown
// ---------------------------------------------------------------------------
void PipelineExecutor::LogWindowStats(WindowStats const& stats) {
  using CodeCounts = std::pair<VariantBuilder::StatusCode const, u64>;
  static auto const SUMMER = [](u64 const sum, CodeCounts const& item) -> u64 {
    return sum + item.second;
  };
  auto const nwindows = std::accumulate(stats.cbegin(), stats.cend(), 0, SUMMER);

  std::ranges::for_each(stats, [&nwindows](CodeCounts const& item) -> void {
    auto const [status_code, count] = item;
    auto const pct_count = (100.0 * static_cast<f64>(count)) / static_cast<f64>(nwindows);
    auto const status_str = ToString(status_code);
    LOG_INFO("{} | {:>8.4f}% of total windows | {} windows", status_str, pct_count, count)
  });
}

// ---------------------------------------------------------------------------
// Execute — orchestrator: allocate → feed → launch → process → shutdown
// ---------------------------------------------------------------------------
auto PipelineExecutor::Execute(std::ostream& output) -> WindowStats {
  auto const num_total = mWindowBuilder.ExpectedTargetWindows();

  static thread_local auto const THREAD_ID = std::this_thread::get_id();
  LOG_INFO("Using main thread {:#x} to synchronize variant calling pipeline",
           absl::Hash<std::thread::id>()(THREAD_ID))
  LOG_INFO("Processing {} window(s) with {} VariantBuilder thread(s)", num_total, mNumThreads)

  mWindows.reserve(num_total);
  mSendQueue = std::make_shared<AsyncWorker::InputQueue>(num_total);
  mRecvQueue = std::make_shared<AsyncWorker::OutputQueue>(num_total);
  mVariantStore = std::make_shared<VariantStore>();

  // Reset per-execution state for potential reuse
  mRegionIdx = 0;
  mWindowStart = -1;
  mGlobalIdx = 0;
  mLastContiguousDone = 0;
  mIdxToFlush = 0;

  moodycamel::ProducerToken const producer_token(*mSendQueue);
  FeedInitialWindows(producer_token, num_total);
  LaunchWorkers();

  auto stats = ProcessAllResults(output, num_total, producer_token);

  ShutdownWorkers();
  mVariantStore->FlushAllVariantsInStore(output);

  return stats;
}

// ---------------------------------------------------------------------------
// FeedInitialWindows — seed the work queue before launching workers
//
// When total windows fit within 2×BATCH_SIZE, generate all windows upfront
// to avoid batched emission overhead for smaller runs.  For larger genomes
// (WGS), emit only the first batch and rely on FeedNextBatch() for the rest.
// ---------------------------------------------------------------------------
void PipelineExecutor::FeedInitialWindows(moodycamel::ProducerToken const& token, usize num_total) {
  static constexpr usize BATCH_THRESHOLD = 2 * WindowBuilder::BATCH_SIZE;
  bool const use_batching = num_total > BATCH_THRESHOLD;

  if (use_batching) {
    // Seed the queue with the first batch of windows
    FeedNextBatch(token);
  } else {
    // Small run: generate all windows upfront, no batching overhead
    mWindows = mWindowBuilder.BuildWindows();
    mGlobalIdx = num_total;  // Mark all as emitted
    mSendQueue->enqueue_bulk(token, mWindows.begin(), mWindows.size());
  }
}

// ---------------------------------------------------------------------------
// FeedNextBatch — emit the next batch from the window builder
//
// Generates the next batch of windows from the builder, enqueues them for
// workers, and appends them to the tracking vector. Callers are responsible
// for checking whether more windows are needed before invoking.
// ---------------------------------------------------------------------------
void PipelineExecutor::FeedNextBatch(moodycamel::ProducerToken const& token) {
  auto next_batch = mWindowBuilder.BuildWindowsBatch(mRegionIdx, mWindowStart, mGlobalIdx);
  if (!next_batch.empty()) {
    mSendQueue->enqueue_bulk(token, next_batch.begin(), next_batch.size());
    mWindows.insert(mWindows.end(), next_batch.begin(), next_batch.end());
  }
}

// ---------------------------------------------------------------------------
// LaunchWorkers — spawn jthreads with cooperative cancellation
//
// Each worker owns a VariantBuilder and processes windows from the shared
// input queue until stop is requested.  Absorbs the former PipelineWorker
// anonymous function.
// ---------------------------------------------------------------------------
void PipelineExecutor::LaunchWorkers() {
  mWorkerThreads.reserve(mNumThreads);
  for (usize idx = 0; idx < mNumThreads; ++idx) {
    mWorkerThreads.emplace_back([send_queue = mSendQueue, recv_queue = mRecvQueue,
                                 variant_store = mVariantStore, var_bldr_params = mParams,
                                 window_length = mWindowLength](std::stop_token stop_token) {
#ifdef LANCET_PROFILE_MODE
      if (ProfilingIsEnabledForAllThreads() != 0) ProfilerRegisterThread();
#endif
      auto worker = std::make_unique<AsyncWorker>(send_queue, recv_queue, variant_store,
                                                  var_bldr_params, window_length);
      worker->Process(std::move(stop_token));
    });
  }
}

// ---------------------------------------------------------------------------
// ShutdownWorkers — cooperative stop + join + profiler cleanup
// ---------------------------------------------------------------------------
void PipelineExecutor::ShutdownWorkers() {
  std::ranges::for_each(mWorkerThreads, std::mem_fn(&std::jthread::request_stop));
  std::ranges::for_each(mWorkerThreads, std::mem_fn(&std::jthread::join));

#ifdef LANCET_PROFILE_MODE
  ProfilerStop();
  ProfilerFlush();
#endif
}

// ---------------------------------------------------------------------------
// FlushCompletedVariants — coordinate-sorted VCF output synchronization
// ---------------------------------------------------------------------------
void PipelineExecutor::FlushCompletedVariants(std::ostream& output, usize num_total,
                                              absl::FixedArray<bool> const& done_windows) {
  // -------------------------------------------------------------------------
  // VCF Output Synchronization & Bulk Flushing
  // -------------------------------------------------------------------------
  // Worker threads finish windows out-of-order. We use `done_windows` to track
  // all completions, and `last_contiguous_done` traces the furthest unbroken
  // chain of sequentially finished windows starting from the beginning.
  while (mLastContiguousDone < num_total && done_windows[mLastContiguousDone]) {
    mLastContiguousDone++;
  }

  // Rather than writing to disk immediately, we lag behind by `NUM_BUFFER_WINDOWS`.
  // This safety gap allows active upstream threads to resolve large structural
  // variants that might overlap across sequential window boundaries.
  static constexpr usize NUM_BUFFER_WINDOWS = 100;
  usize target_flush_idx = 0;
  if (mLastContiguousDone > NUM_BUFFER_WINDOWS) {
    target_flush_idx = mLastContiguousDone - NUM_BUFFER_WINDOWS;
  }

  // varstore->FlushVariantsBeforeWindow takes a target window and dumps ALL
  // buffered variants chronologically up to that genomic threshold in one shot.
  //
  // Example Scenario (If NUM_BUFFER_WINDOWS = 2):
  //  - `last_contiguous_done` evaluates to 3.
  //  - Thread A finishes window #5 (done_windows[5] = true, chain unbroken, cursor stays 3).
  //  - Thread B finishes window #4 (done_windows[4] = true, chain completes!).
  //  - `last_contiguous_done` instantly slides from 3 up to 5.
  //  - `target_flush_idx` evaluates to (5 - 2) = 3.
  //  - `idx_to_flush` jumps from its old state directly to 3.
  //  - A single FlushVariantsBeforeWindow(*windows[3]) call fires, sweeping
  //    the store and safely dumping all variants prior to window #3 to disk.
  if (mIdxToFlush < target_flush_idx) {
    mIdxToFlush = target_flush_idx;
    mVariantStore->FlushVariantsBeforeWindow(*mWindows[mIdxToFlush], output);
  }
}

// ---------------------------------------------------------------------------
// ProcessAllResults — main event loop
//
// Dequeues worker results, tracks progress via ETA timer, dynamically feeds
// new window batches when queue runs low, and flushes completed variants in
// coordinate order.
// ---------------------------------------------------------------------------
auto PipelineExecutor::ProcessAllResults(std::ostream& output, usize num_total,
                                         moodycamel::ProducerToken const& token) -> WindowStats {
  auto stats = InitWindowStats();

  usize num_completed = 0;
  absl::FixedArray<bool> done_windows(num_total);
  done_windows.fill(false);

  AsyncWorker::Result result;
  moodycamel::ConsumerToken consumer_token(*mRecvQueue);

  base::Timer timer;
  EtaTimer eta_timer(num_total);

  // ---------------------------------------------------------------------------
  // Main pipeline loop: process results and dynamically feed new window batches
  // ---------------------------------------------------------------------------
  // NOTE: The feed_next_batch check operates unconditionally at the top of the loop.
  // This guarantees that worker queues are consistently supplied regardless of
  // whether the subsequent queue read succeeds or times out.
  //
  // wait_dequeue_timed serves as a blocking futex call preventing hardware thread-spin.
  // If the timeout triggers without a payload, we `continue` the loop,
  // allowing the top-level feed evaluation to re-poll state.
  constexpr auto QUEUE_TIMEOUT = std::chrono::milliseconds(10);
  static constexpr auto ELAPSED_PRECISION = absl::Seconds(1);
  static constexpr auto WINDOW_RT_PRECISION = absl::Microseconds(100);
  while (num_completed != num_total) {
    if (mGlobalIdx < num_total && mSendQueue->size_approx() < WindowBuilder::BATCH_SIZE) {
      FeedNextBatch(token);
    }

    // NOTE: Sleep is handled by the wait_dequeue_timed futex block.
    if (!mRecvQueue->wait_dequeue_timed(consumer_token, result, QUEUE_TIMEOUT)) continue;

    num_completed++;
    stats.at(result.mStatus) += 1;
    done_windows[result.mGenomeIdx] = true;

    eta_timer.Increment();
    auto const percent = 100.0 * static_cast<f64>(num_completed) / static_cast<f64>(num_total);
    auto const window_name = mWindows[result.mGenomeIdx]->ToSamtoolsRegion();
    auto const elapsed = FormatTruncatedDuration(timer.Runtime(), ELAPSED_PRECISION);
    auto const remaining = FormatTruncatedDuration(eta_timer.EstimatedEta(), ELAPSED_PRECISION);
    auto const window_runtime = FormatTruncatedDuration(result.mRuntime, WINDOW_RT_PRECISION);
    LOG_INFO("Progress: {:>8.4f}% | Elapsed: {} | ETA: {} @ {:.2f}/s | {} done with {} in {}",
             percent, elapsed, remaining, eta_timer.RatePerSecond(), window_name,
             ToString(result.mStatus), window_runtime)

    FlushCompletedVariants(output, num_total, done_windows);
  }

  return stats;
}

}  // namespace lancet::core
