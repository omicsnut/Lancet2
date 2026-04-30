---
description: Lancet2 core/ layer rules — pipeline orchestration with lock-free queues, 256-shard mutex-protected VariantStore, the chunked-sorted VCF flush logic that produces tail-readable coordinate-sorted output, jthread cooperative cancellation, per-thread VariantBuilder reuse. Load when editing src/lancet/core/.
paths:
  - "src/lancet/core/**"
---

# core/ layer rules

`core/` is the pipeline orchestrator. It owns the worker pool, the
input/output queues, the variant store, and — most subtly — the
synchronization logic that produces a coordinate-sorted, incrementally-
flushed BGZF VCF from out-of-order parallel window completions. This
is where most concurrency bugs hide; the design is conservative on
purpose.

## VariantStore is sharded (256 buckets, alignas(64)) for write contention

`VariantStore` carries an `std::array<VariantBucket, 256>` where each
`VariantBucket` is `alignas(64)` (cache-line aligned, no false
sharing) and contains its own `absl::Mutex` plus a
`flat_hash_map<Key, Value>` annotated with `ABSL_GUARDED_BY(mMutex)`.
Shard selection is `absl::Hash<Key>()(identifier) & (NUM_SHARDS - 1)`
— a power-of-two mask, not modulo.

Workers contend at most on the bucket containing their variant; with
256 shards and ~30 worker threads, contention is statistically
negligible. **Do not** lower `NUM_SHARDS` below 256 without profile
data — the cost is one mutex and one hash map per shard (cheap), and
the benefit is essentially lock-free behavior at high thread counts.

The `alignas(64)` is load-bearing: without it, multiple buckets
share a cache line, and updating one bucket's mutex/map invalidates
the cache line for adjacent buckets — cache-line ping-pong silently
serializes the workers.

## Duplicate-variant policy: keep the higher-coverage version

`VariantStore::AddVariants` handles duplicates (same `CHROM+POS+REF`
identifier) by keeping the one with higher `TotalCoverage()`. The
rationale: when two windows overlap and both call the same locus,
the better-covered window assembled a more complete multi-allelic
picture. Replacing this with "first wins" or "last wins" silently
loses ALTs at the overlap boundary.

## The chunked sorted VCF flush is a three-component contract

The output VCF must be:
1. **Coordinate-sorted** by `(CHROM, POS)` — required by tabix
   indexing and every downstream consumer.
2. **Incrementally flushed** to disk during the run (not buffered
   until the end) — so the user can `tail -f` the output and so
   memory doesn't blow up on chr1.
3. **Correct under out-of-order parallel completion** — workers
   finish windows out of order, so naive "write as you finish"
   produces unsorted output.

The logic that satisfies all three lives in
`PipelineExecutor::FlushCompletedVariants` plus
`VariantStore::FlushVariantsBeforeWindow`:

```
mLastContiguousDone  ← tracks the longest unbroken prefix of
                       finished windows (slides forward only when
                       prefix is contiguous)
target_flush_idx     ← mLastContiguousDone − NUM_BUFFER_WINDOWS
                       (deliberately lags by 100 windows)
FlushVariantsBeforeWindow(*windows[target_flush_idx])
                     ← sweeps all 256 shards, extracts variants
                       whose (ChromIndex, StartPos1) is before the
                       target window, sorts in-place by VariantCall::
                       operator<, fmt::print's the VCF records,
                       calls out.flush() on the BgzfOstream
```

**`NUM_BUFFER_WINDOWS = 100` is load-bearing**: structural variants
can span window boundaries, and the 100-window buffer ensures the
next window's results have been written into the store before we
flush variants from earlier windows. Lowering it produces missing
SVs at window boundaries; raising it pushes more memory into the
store without correctness benefit.

The `mLastContiguousDone` cursor "slides instantly" when the missing
predecessor finishes — example: workers complete 5, 4, then the
cursor jumps from 3 to 5 in one update. A single
`FlushVariantsBeforeWindow` call sweeps the store. Don't try to
flush per-window; the batch is what makes the in-place sort cheap.

`FlushAllVariantsInStore` runs once at pipeline shutdown to drain
remaining variants. Without this, the last `NUM_BUFFER_WINDOWS`
variants are never written.

## VariantCall::operator< is the sort key — don't reimplement

The flush sort uses `*lhs < *rhs` from `caller::VariantCall::
operator<`. This compares `(ChromIndex, StartPos1, REF, ALT)` —
match this convention exactly when comparing variants elsewhere.
Code that sorts by `StartPos1` alone produces correctly chromosome-
chunked but wrong intra-chromosome order.

## HAS_NO_SUPPORT filter excludes spurious calls

Both flush paths apply `HAS_NO_SUPPORT`: skip variants where
`!HasAltSupport()` or every category is `AlleleType::REF`. These
appear as flushable entries because the genotyper records them
during processing, but they have no informative content. Don't
remove the filter — the VCF gets bloated with no-op records.

## Workers communicate via lock-free queues, not direct calls

`PipelineExecutor` owns:
- `mSendQueue` (`moodycamel::BlockingConcurrentQueue<WindowPtr>`) —
  main thread → workers
- `mRecvQueue` (`moodycamel::BlockingConcurrentQueue<Result>`) —
  workers → main thread

`wait_dequeue_timed(consumer_token, result, 10ms)` is the futex
block — no hardware spin. The 10ms timeout on the main loop's
`mRecvQueue->wait_dequeue_timed` is intentional: it lets the loop
top re-check whether to feed the next batch even when no result is
ready.

Workers (`AsyncWorker::Process`) use the same pattern with their own
`mInPtr->wait_dequeue_timed`, plus a `stop_token` check so cooperative
cancellation works (`std::jthread::request_stop` from the executor's
shutdown).

## Per-thread VariantBuilder, with deterministic worker_index

Each worker thread owns one `VariantBuilder` instance — never
shared. `worker_index_for_thread = static_cast<u32>(thread_idx)` is
forwarded to the builder so per-thread artifact filenames (graph DOT
shards) are deterministic across runs. Without deterministic
indexing, debugging artifact ordering across runs is impossible.

## Crash context registration is per-thread, scoped per window

Each worker registers a thread-local crash slot at startup
(`RegisterThreadSlot()`), populates it before each `ProcessWindow`
with `(genome_index, region_str)`, clears it after, and unregisters
at exit. The crash handler reads these slots to print which window
each worker was processing when the crash hit. **The
`ClearSlotWindowInfo` in the catch handler before `std::terminate()`
is required** — without it, a downstream signal handler reads stale
context.

## WindowBuilder constants encode workload assumptions

`DEFAULT_WINDOW_LENGTH = 1000`, `DEFAULT_PCT_OVERLAP = 20`,
`DEFAULT_REGION_PADDING = 500`. The bounds enforce sanity:
windows in `[1000, 2500]`, overlap in `[10, 90]%`, padding ≤ 1000.
Lowering window size below 1000 makes the per-window overhead
dominate; raising it past 2500 makes the assembly memory blow up on
repetitive regions. The defaults reflect chr1/chr4 profiling.

`BATCH_SIZE` (in `WindowBuilder.h`) controls how many windows are
generated up-front before the executor starts pulling more.
`BATCH_THRESHOLD = 2 * BATCH_SIZE` is the cutoff at which
`PipelineExecutor::FeedInitialWindows` switches from upfront-
generate-all to streaming-batches. For small runs (single region),
upfront is cheaper; for WGS, batched is the only feasible path.

## ReadCollector caps coverage to bound memory

`ReadCollector::DEFAULT_MAX_WINDOW_COVERAGE = 1000.0`. Windows in
ultra-deep regions (PCR duplicates, centromeres) get downsampled to
1000× before assembly. Without this cap, the graph for a 100,000×
window has millions of nodes and the worker times out. The
downsampling preserves paired-end mate pairs (uniform paired
downsampling) — don't replace with single-end downsampling, which
breaks the genotyper's mate-distance features.

## ProcessWindow exceptions terminate the run, not the worker

When `VariantBuilder::ProcessWindow` throws, `AsyncWorker::Process`
logs `LOG_CRITICAL` with the window context and calls
`std::terminate()`. **One bad window kills the pipeline.** This is
intentional: silently skipping a window produces a VCF with a
mysterious gap, and partial output is worse than no output. Don't
"helpfully" wrap the worker in a try/catch that continues — every
window that fails should be visible to the user.
