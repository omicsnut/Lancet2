---
name: sanitizer-triage
description: Use when investigating a crash, intermittent test failure, suspected data race, memory error, or any sanitizer-detectable bug, AND the goal is to land a fix with a regression test. Trigger on "this crashes", "race condition", "use-after-free", "ASan says", "TSan reports", "MSan", "UBSan", "intermittent flake in tests/". Drives pixi-managed sanitizer build trees (cmake-build-asan, cmake-build-msan, cmake-build-tsan, cmake-build-ubsan), reproduces the failure, analyzes the trace, and lands the minimum fix. Per-sanitizer detail (what each catches, runtime options, common findings, the mimalloc-coverage caveat, Clang doc URLs) lives in references/sanitizer_matrix.md, loaded on demand. For analysis-only of an existing report, use the sanitizer-triage subagent instead.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Sanitizer triage on Lancet2

This skill formalizes the end-to-end procedure for taking a sanitizer-detectable bug from "I have a reproduction" to "merged with regression test." It composes with the `sanitizer-triage` subagent for analysis; this skill is the workflow the calling session executes.

The workflow uses pixi-managed sanitizer build tasks (`pixi run configure-asan`, `build-asan`, `test-asan`, and the parallel families for MSan, TSan, UBSan). These tasks live in `pixi.toml`. Pipeline-level sanitizer runs additionally require the `LANCET_SANITIZE_BUILD` source guard in `src/lancet/main.cpp`, `benchmarks/main.cpp`, and `CMakeLists.txt` — without that guard, sanitizer builds either fail to link or run silently because Lancet2's mimalloc allocator override bypasses the sanitizer interceptors. See "Setup verification" below to confirm both pieces are in place before starting.

The per-sanitizer detail — which sanitizer catches what, runtime option semantics, common findings, the mimalloc coverage caveat, and the Clang documentation URLs — lives in `references/sanitizer_matrix.md`. Load it on demand at Step 1 and Step 4.

## Setup verification (one-time, per repository)

Before this skill works, two things must be in place — both are part of the project's standard configuration:

1. **Sanitizer pixi tasks present.** Verify with `pixi task list | grep -E '(asan|msan|tsan|ubsan)'`. Twelve task names should appear: `configure-{asan,msan,tsan,ubsan}`, `build-{asan,msan,tsan,ubsan}`, `test-{asan,msan,tsan,ubsan}`. If any are missing, the project's `pixi.toml` is out of date for sanitizer support — surface this to the user and stop; don't try to hand-roll cmake invocations.

2. **`LANCET_SANITIZE_BUILD` source guard present (only required for pipeline-level runs).** Verify with `grep -l LANCET_SANITIZE_BUILD CMakeLists.txt src/lancet/main.cpp benchmarks/main.cpp`. All three files should match. The guard gates Lancet2's mimalloc allocator override behind a CMake option that the sanitizer pixi tasks set; without it, sanitizer builds of the `Lancet2` pipeline binary either fail to link (`multiple definition of malloc`) or run but report nothing (mimalloc bypasses the sanitizer interceptors).

If you only need to sanitize unit tests (`pixi run test-asan`, etc.), the source guard is NOT required — `TestLancet2` does not link mimalloc. The guard is required when sanitizer-running the pipeline binary itself.

## Step 1 — Pick the right sanitizer

Match the symptom to the sanitizer. The matrix in `references/sanitizer_matrix.md` covers each in detail; the short version:

- **ASan** (the default first try): crashes with `SEGV`, `heap-buffer-overflow`, `use-after-free`, `stack-buffer-overflow`. Linux ASan also includes LSan for leaks. Cheap (2× CPU); always combined with UBSan in the project's `configure-asan` task.
- **TSan**: intermittent test failures, hangs, output that varies between runs even with deterministic input. Roughly 5-15× CPU; needs multiple threads to actually surface races (use `--num-threads $(nproc)` even for normally-2-thread bugs).
- **MSan**: suspected uninitialized-memory reads in pure Lancet2 code paths. **Use sparingly** — MSan produces false positives on every call into vendored deps (htslib, abseil, spdlog, SPOA, minimap2) because those aren't MSan-instrumented. See `references/sanitizer_matrix.md` for the full caveat.
- **UBSan-standalone**: inconsistent arithmetic, "shift exponent too large", "signed integer overflow", alignment-dependent failures. Already bundled into the `configure-asan` task; reach for `configure-ubsan` only when ASan's allocator interposition is interfering with the bug.

If the symptom doesn't clearly match one, start with ASan. It catches the broadest class and has the best output quality.

## Step 2 — Build under the chosen sanitizer

```bash
pixi run build-asan   # or build-msan, build-tsan, build-ubsan
```

The `build-X` task depends on `configure-X`, so a fresh tree configures and builds in one command. Subsequent runs incrementally rebuild — pixi's dependency tracker skips configure if `cmake-build-X/CMakeCache.txt` exists.

The `cmake-build-asan/`, `cmake-build-msan/`, `cmake-build-tsan/`, `cmake-build-ubsan/` directories are matched by the protected-paths hook's `cmake-build-*` glob — agent writes are blocked but build-system writes pass through (the hook only fires on `Edit|Write|MultiEdit` tool calls). Read access from the agent into these trees is unrestricted, so reading vendored-dep source under `cmake-build-asan/_deps/` for triage works as expected.

## Step 3 — Reproduce under the sanitizer

For unit-test reproduction (works on a stock checkout):

```bash
pixi run test-asan -- "[failing test name]"   # or test-msan/test-tsan/test-ubsan
```

The `--` separates the pixi task args from arguments forwarded to `TestLancet2`. The trailing `[tag]` filter limits the run to one Catch2 test, which is what you want during triage.

For pipeline-level reproduction (requires the `LANCET_SANITIZE_BUILD` source guard per Setup verification §2):

```bash
pixi run ./cmake-build-asan/Lancet2 pipeline \
  --tumor $LANCET_TEST_SOMATIC_TUMOR \
  --normal $LANCET_TEST_SOMATIC_NORMAL \
  --reference $LANCET_TEST_SOMATIC_REFERENCE \
  --region $LANCET_TEST_SOMATIC_REGION_SMALL \
  --num-threads $(nproc) \
  --out-vcfgz /tmp/triage.vcf.gz \
  2> /tmp/sanitizer.log
```

The `pixi run` prefix is required for sanitizer pipeline runs because sanitizer builds are dynamically linked (`-DLANCET_BUILD_STATIC=OFF`); the sanitizer runtime libraries (`libasan.so`, `libtsan.so`, etc.) and `libc++` resolve from the pixi env's library path. Without `pixi run`, the binary fails to start with a missing-library error. This differs from production and profiling builds, which are statically linked and run fine outside the pixi env.

Use `_REGION_SMALL` (1 Mb) regions to keep iteration fast — sanitizer-instrumented binaries are 2-15× slower depending on which sanitizer (the per-sanitizer overhead estimates are in `references/sanitizer_matrix.md`). Save the sanitizer output to a file (`2> /tmp/sanitizer.log` above) — output is verbose and easy to lose to terminal scrollback.

If the bug doesn't reproduce on first try:
- TSan: increase thread count, run multiple times (10×) to catch intermittents.
- ASan/UBSan: try `ASAN_OPTIONS=detect_stack_use_after_return=1:check_initialization_order=1` for additional checks (see `references/sanitizer_matrix.md`).
- Any sanitizer: try a different region, or a different fixture (germline NA12878 vs somatic HCC1395).

## Step 4 — Analyze with the sanitizer-triage subagent

Hand the captured sanitizer output to the subagent:

```
Use the sanitizer-triage subagent to analyze this output:
[paste the sanitizer report]

The reproduction is: [your command]
```

The subagent reads the cited frames in `src/lancet/`, identifies the root cause, and proposes the minimum fix. It has its own memory file at `.claude/agent-memory/sanitizer-triage.md` recording known-benign warnings and patterns of new warnings that turned out to be real bugs — consult that before reasoning from scratch.

For unfamiliar diagnostic categories ("what does this UBSan report mean?"), `references/sanitizer_matrix.md` has per-sanitizer notes including option semantics and Lancet2-specific gotchas.

## Step 5 — Write a regression test FIRST

Before implementing the fix, add a Catch2 test in `tests/<layer>/` that reproduces the bug deterministically. Use the `add-cpp-test` skill for the canonical procedure. The test must fail under the sanitizer build before the fix is applied.

For race conditions, deterministic reproduction is hard but not impossible. Patterns that work: constraining thread count to 2-4 (small enough for the test to be fast, large enough for contention), injecting `std::this_thread::yield()` calls at the suspected race window, using `absl::Mutex::AssertHeld` to verify locking invariants.

## Step 6 — Apply the minimum fix

In `src/lancet/<layer>/`, make the smallest change the subagent proposed. The fix must address the root cause, not the symptom. Common fix patterns:
- Per-thread instances of a shared resource (most HTSlib iterator races).
- Explicit synchronization on a previously-implicit shared variable (most variant-store races).
- Bounds-checking on a previously-trusted index (most ASan findings).
- Proper RAII for an HTSlib resource (most leak findings).

## Step 7 — Verify under sanitizer and under release

Run the regression test under the sanitizer build to confirm the fix:

```bash
pixi run test-asan -- "[regression test name]"
```

Then verify the fix doesn't regress the release build:

```bash
pixi run build-release
pixi run test
pixi run lint-check
```

A fix that passes the sanitizer but fails clang-tidy needs further work.

## Step 8 — Commit

Use the conventional-commit prefix `fix:` and reference the sanitizer in the body:

```
fix: per-thread BAM iterator in read_collector

Reproduced under TSan (cmake-build-tsan, --num-threads 8) on
chr1:1000000-2000000. Race was on shared bam1_t* between worker
threads — HTSlib iterators are not thread-safe. Fix: per-worker
iterator instance with RAII cleanup.

Regression test: tests/hts/read_collector_test.cpp
```

Lancet2's `.chglog/config.yml` does not support scopes; mention the layer in the subject text instead.

## Step 9 — Optional: keep the sanitizer tree

The `cmake-build-{asan,msan,tsan,ubsan}/` trees are excluded by `.claudeignore` and `.gitignore` (the at-write hook protects `cmake-build-*`). You can keep them around for future triage — incremental rebuilds via `pixi run build-asan` are much faster than re-configuring from scratch.

## When NOT to use this skill

Do not use this skill for bugs that aren't sanitizer-detectable (logic errors, wrong VCF schema, race-free correctness gaps); use the standard `add-cpp-test` workflow. Do not use it as a routine pre-merge check — the project's CI already runs ASan on the test suite. Use this skill specifically when you have a symptom and need to find the root cause.

## When the bug is in a vendored dependency

If the subagent identifies the root cause inside HTSlib, minimap2, SPOA, longdust, or another vendored library under `cmake-build-*/_deps/`, the fix path changes. Do NOT patch the vendored code (which is write-blocked anyway by protected-paths and would be lost on next reconfigure). Instead, work around it at the Lancet2 callsite. Common workarounds:
- **Per-thread instances** (for thread-unsafe libraries — most HTSlib uses).
- **Defensive copying** (for libraries with unstable iterator semantics).
- **Pre-validation** (for libraries with insufficient input checking).

Document the workaround with a comment citing the upstream behavior. The dependency's source remains freely readable (read-only protection — see AGENTS.md "What you must not do") so you can confirm the workaround is necessary, not just convenient.

## A note on mimalloc coverage

Sanitizer builds use libc malloc, not Lancet2's production mimalloc (the `LANCET_SANITIZE_BUILD` source guard makes this unconditional in sanitizer trees). A clean ASan run does NOT prove production-with-mimalloc is clean. Bugs that depend on mimalloc's allocation pattern, thread-local heap layout, or page-recycling timing will not surface here. The trade-off is intentional — mimalloc's preprocessor-level override of `malloc/free/new/delete` bypasses sanitizer interceptors entirely, so the choice is "sanitize without mimalloc" or "don't sanitize at all." The guard's rationale and exact scope are documented inline in `src/lancet/main.cpp` and `CMakeLists.txt`; read them when adjusting the sanitizer build configuration.

## References

- `references/sanitizer_matrix.md` — Per-sanitizer detail for ASan, MSan, TSan, UBSan, and LSan: what each catches, build flag rationale, runtime options, common findings, Lancet2-specific notes (per-thread BAM iterators, variant-store sharding, vendored-dep coverage, the mimalloc situation), and the canonical Clang documentation URL for each. Load this at Step 1 (picking the sanitizer) and Step 4 (interpreting unfamiliar output).
