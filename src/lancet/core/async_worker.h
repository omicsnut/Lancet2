#ifndef SRC_LANCET_CORE_ASYNC_WORKER_H_
#define SRC_LANCET_CORE_ASYNC_WORKER_H_

#include "lancet/base/types.h"
#include "lancet/core/variant_builder.h"
#include "lancet/core/variant_store.h"
#include "lancet/core/window.h"

#include "absl/time/time.h"
#include "blockingconcurrentqueue.h"

#include <memory>
#include <stop_token>
#include <utility>

namespace lancet::core {

class AsyncWorker {
 public:
  struct Result {
    // ── 8B Align ────────────────────────────────────────────────────────────
    usize mGenomeIdx = 0;
    absl::Duration mRuntime = absl::ZeroDuration();
    // ── 1B Align ────────────────────────────────────────────────────────────
    VariantBuilder::StatusCode mStatus = VariantBuilder::StatusCode::UNKNOWN;
  };

  using InputQueue = moodycamel::BlockingConcurrentQueue<WindowPtr>;
  using OutputQueue = moodycamel::BlockingConcurrentQueue<Result>;

  using InQueuePtr = std::shared_ptr<InputQueue>;
  using OutQueuePtr = std::shared_ptr<OutputQueue>;
  using VariantStorePtr = std::shared_ptr<VariantStore>;
  using VariantBuilderPtr = std::unique_ptr<VariantBuilder>;
  using BuilderParamsPtr = std::shared_ptr<VariantBuilder::Params const>;

  /// `worker_index` identifies this worker within the pool; range is `[0, num_threads)`
  /// and PipelineExecutor assigns indices in construction order. Forwarded to VariantBuilder
  /// so each worker's per-thread graph shard gets a deterministic filename.
  AsyncWorker(InQueuePtr in_queue_ptr, OutQueuePtr out_queue_ptr, VariantStorePtr variant_store_ptr,
              BuilderParamsPtr params, u32 window_len, u32 worker_id)
      : mInPtr(std::move(in_queue_ptr)),
        mOutPtr(std::move(out_queue_ptr)),
        mStorePtr(std::move(variant_store_ptr)),
        mBuilderPtr(std::make_unique<VariantBuilder>(std::move(params), window_len, worker_id)) {}

  void Process(std::stop_token stop_token);

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  InQueuePtr mInPtr;
  OutQueuePtr mOutPtr;
  VariantStorePtr mStorePtr;
  VariantBuilderPtr mBuilderPtr;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_ASYNC_WORKER_H_
