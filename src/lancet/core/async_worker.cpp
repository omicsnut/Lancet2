#include "lancet/core/async_worker.h"

#include "lancet/base/logging.h"
#include "lancet/base/timer.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"

#include "concurrentqueue.h"

#include <absl/hash/hash.h>
#include <blockingconcurrentqueue.h>
#include <chrono>
#include <cxxabi.h>
#include <exception>
#include <memory>
#include <stop_token>
#include <thread>
#include <typeinfo>
#include <utility>
#include <vector>

#include <cstdlib>

namespace lancet::core {

// ============================================================================
// AsyncWorker::Process — thread pool worker loop
//
// Each AsyncWorker runs in its own std::jthread, pulling windows from a
// lock-free MPMC queue (moodycamel::BlockingConcurrentQueue).  The loop:
//   1. Check stop_token — cooperative cancellation from the main thread
//   2. Dequeue a window (10ms timeout — prevents busy-spinning)
//   3. Run VariantBuilder::ProcessWindow → assemble, call, genotype
//   4. Push result back to the output queue for progress reporting
// ============================================================================
// stop_token is a lightweight handle designed for by-value pass per C++20 jthread API
// NOLINTNEXTLINE(performance-unnecessary-value-param)
void AsyncWorker::Process(std::stop_token stop_token) {
  static thread_local auto const THREAD_ID =
      absl::Hash<std::thread::id>()(std::this_thread::get_id());
  LOG_DEBUG("Starting AsyncWorker thread {:#x}", THREAD_ID)

  lancet::base::Timer timer;
  usize num_done = 0;
  auto window_ptr = std::make_shared<Window>();
  moodycamel::ProducerToken const out_token(*mOutPtr);
  constexpr auto QUEUE_TIMEOUT = std::chrono::milliseconds(10);

  while (true) {
    if (stop_token.stop_requested()) break;

    // Blocking dequeue with timeout — prevents busy-spinning while allowing
    // periodic re-check of the stop_token.
    if (!mInPtr->wait_dequeue_timed(window_ptr, QUEUE_TIMEOUT)) continue;

    timer.Reset();
    try {
      auto variants = mBuilderPtr->ProcessWindow(std::const_pointer_cast<Window const>(window_ptr));
      mStorePtr->AddVariants(std::move(variants));
    } catch (std::exception const& exc) {
      LOG_CRITICAL("AsyncWorker thread {:#x} CRASHED on window idx={}: {}", THREAD_ID,
                   window_ptr->GenomeIndex(), exc.what())
      std::abort();
    } catch (...) {
      // abi::__cxa_current_exception_type() gives the mangled type name of
      // whatever was thrown (int, char*, custom class, etc.). A rethrow into
      // catch(std::exception&) is redundant — the catch above already covers
      // all std::exception subclasses.
      char const* type_name = "unknown";
#if defined(__GNUC__) || defined(__clang__)
      if (auto const* type_info = abi::__cxa_current_exception_type()) {
        type_name = type_info->name();
      }
#endif
      LOG_CRITICAL("AsyncWorker thread {:#x} CRASHED on window idx={}: non-std exception type={}",
                   THREAD_ID, window_ptr->GenomeIndex(), type_name)
      std::abort();
    }

    auto const status_code = mBuilderPtr->CurrentStatus();
    mOutPtr->enqueue(out_token, Result{.mGenomeIdx = window_ptr->GenomeIndex(),
                                       .mRuntime = timer.Runtime(),
                                       .mStatus = status_code});
    num_done++;
  }

  LOG_DEBUG("Quitting AsyncWorker thread {:#x} after processing {} windows", THREAD_ID, num_done)
}

}  // namespace lancet::core
