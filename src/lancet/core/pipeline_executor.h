#ifndef SRC_LANCET_CORE_PIPELINE_EXECUTOR_H_
#define SRC_LANCET_CORE_PIPELINE_EXECUTOR_H_

#include "lancet/base/types.h"
#include "lancet/core/async_worker.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/variant_store.h"
#include "lancet/core/window.h"
#include "lancet/core/window_builder.h"

#include "absl/container/btree_map.h"
#include "absl/container/fixed_array.h"
#include "concurrentqueue.h"

#include <iosfwd>
#include <memory>
#include <thread>
#include <vector>

namespace lancet::core {

/// Owns the batch-fed, multi-threaded window processing pipeline.
/// Extracted from PipelineRunner::Run() to isolate the execution engine from
/// CLI concerns (param validation, I/O setup, VCF header construction).
class PipelineExecutor {
 public:
  using WindowStats = absl::btree_map<VariantBuilder::StatusCode, u64>;

  PipelineExecutor() = delete;
  PipelineExecutor(WindowBuilder builder, std::shared_ptr<VariantBuilder::Params const> params,
                   usize num_threads, u32 window_length);

  /// Run the full pipeline execution lifecycle:
  /// feed windows → launch workers → process results → flush variants → shutdown.
  /// Returns per-status-code window counts for downstream logging.
  [[nodiscard]] auto Execute(std::ostream& output) -> WindowStats;

  /// Log final window status breakdown to the application logger.
  static void LogWindowStats(WindowStats const& stats);

 private:
  // ── 8B Align (construction-time state) ──────────────────────────────────
  WindowBuilder mWindowBuilder;                           // 8B+ — owns region → window partitioning
  std::shared_ptr<VariantBuilder::Params const> mParams;  // 8B  — shared immutable build params
  usize mNumThreads;                                      // 8B  — number of worker jthreads

  // ── 8B Align (per-execution state, initialized in Execute) ──────────────
  std::vector<WindowPtr> mWindows;                       // 8B+ — all windows emitted so far
  std::shared_ptr<AsyncWorker::InputQueue> mSendQueue;   // 8B  — lock-free producer → consumer
  std::shared_ptr<AsyncWorker::OutputQueue> mRecvQueue;  // 8B  — lock-free consumer → producer
  std::shared_ptr<VariantStore> mVariantStore;           // 8B  — thread-safe variant dedup store
  std::vector<std::jthread> mWorkerThreads;              // 8B+ — C++20 cooperative cancellation

  // ── 8B Align (batch feeding state) ──────────────────────────────────────
  usize mRegionIdx = 0;   // 8B  — current region in builder
  i64 mWindowStart = -1;  // 8B  — current window start offset
  usize mGlobalIdx = 0;   // 8B  — next global window index

  // ── 8B Align (contiguous flush tracking) ────────────────────────────────
  usize mLastContiguousDone = 0;  // 8B  — watermark of contiguously done
  usize mIdxToFlush = 0;          // 8B  — last flushed window index

  // ── 4B Align ────────────────────────────────────────────────────────────
  u32 mWindowLength;  // 4B  — cached from params for workers

  /// Enqueue the next batch of windows from the window builder.
  /// Called when the send queue drops below BATCH_SIZE to keep workers fed.
  void FeedNextBatch(moodycamel::ProducerToken const& token);

  /// Feed all windows upfront (small genomes) or the initial batch (WGS).
  /// Uses BuildWindows() for small runs, BuildWindowsBatch() for large runs.
  void FeedInitialWindows(moodycamel::ProducerToken const& token, usize num_total);

  /// Launch num_threads jthreads, each owning a VariantBuilder instance.
  /// Uses C++20 cooperative cancellation via std::stop_token.
  void LaunchWorkers();

  /// Request stop on all worker threads and join them.
  /// Also stops gperftools profiler in LANCET_PROFILE_MODE.
  void ShutdownWorkers();

  /// Advance the contiguous-done watermark and flush completed variants
  /// to the output stream, maintaining coordinate-sorted VCF output.
  void FlushCompletedVariants(std::ostream& output, usize num_total,
                              absl::FixedArray<bool> const& done_windows);

  /// The main event loop: dequeue results, track progress, flush variants.
  /// Returns accumulated per-status-code window counts.
  [[nodiscard]] auto ProcessAllResults(std::ostream& output, usize num_total,
                                       moodycamel::ProducerToken const& token) -> WindowStats;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_PIPELINE_EXECUTOR_H_
