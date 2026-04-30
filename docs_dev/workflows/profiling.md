# Lancet2 CPU Profiling Guide

This guide covers the end-to-end workflow for profiling Lancet2: building with profiling enabled, collecting CPU profiles, analyzing them with the profiling toolkit, diffing between runs, and tracking performance trends across commits. All profiling tasks run in an isolated pixi environment (`-e profiling`) to avoid polluting the default build environment with Python dependencies.

---

## Prerequisites

1. **pixi** installed and working (`pixi --version`)
2. **Profiling environment** installed: `pixi install -e profiling`
3. **Lancet2** built in Release mode with profiling enabled (see [Building with Profiling](#building-with-profiling))

---

## Building with Profiling

Lancet2 uses [gperftools](https://github.com/gperftools/gperftools) for CPU profiling. The build system provides dedicated profiling tasks that automatically configure the correct build type:

```bash
pixi run build-profile
```

This is equivalent to the expanded form:

```bash
pixi run configure-profile    # sets RelWithDebInfo + LANCET_PROFILE_MODE=ON
pixi run build-profile        # builds with profiling enabled
```

The profiling build uses `CMAKE_BUILD_TYPE=RelWithDebInfo`, which applies the same optimization flags as Release (`-O3`, architecture targeting, `-ffast-math`) plus `-g` for DWARF debug symbols. This preserves full-speed codegen while enabling pprof's `--list` source-line annotation — a critical capability for sub-function bottleneck resolution.

If `LANCET_PROFILE_MODE=ON` is passed with any other build type (e.g., Release), `platform_checks.cmake` automatically overrides to `RelWithDebInfo` and prints a status message.

The build links gperftools `libprofiler.a` (CPU profiler only — Lancet2 uses mimalloc as its allocator, not tcmalloc) into the `lancet_core` and `lancet_cli` targets. Profiling overhead is negligible (~2–5%).

> [!WARNING]
> **Do NOT profile Debug builds.** Debug builds have optimizations disabled (`-O0`), inlined functions are not inlined, and timing data is meaningless for performance analysis. Always use `pixi run build-profile` for profiling.

---

## Collecting a Profile

When built with `LANCET_PROFILE_MODE=ON`, profiling is fully automatic — no environment variables needed. The `PipelineRunner` constructor:

1. Sets `CPUPROFILE_PER_THREAD_TIMERS=1` (profile all worker threads, not just main)
2. Sets `CPUPROFILE_FREQUENCY=250` (250 samples/second — 2.5× the gperftools default of 100 Hz; finer enough to resolve short hot functions, coarse enough that overhead stays negligible)
3. Calls `ProfilerStart()` with an auto-generated timestamped filename

Each worker thread registers itself via `ProfilerRegisterThread()` in `PipelineExecutor::LaunchWorkers()`. When the pipeline completes, `PipelineExecutor::ShutdownWorkers()` calls `ProfilerStop()` and `ProfilerFlush()`.

Run Lancet2 normally — the profile is written automatically:

```bash
./cmake-build-relwithdebinfo/Lancet2 pipeline \
    --normal normal.bam --tumor tumor.bam --reference ref.fa \
    --region chr1:1000000-2000000 --out-vcfgz test.vcf.gz --num-threads 8
```

The output file is written to the **current working directory** with an auto-generated name:

```
Lancet.cpu_profile.<YYYYMMDD>T<HHMMSS>.bin
```

The timestamp uses the local timezone and is formatted by `absl::FormatTime("%Y%m%d%ET%H%M%S")`. For example, a run started at 1:46:21 PM on April 16, 2026 produces:

```
Lancet.cpu_profile.20260416T134621.bin
```

Typical file sizes are **2–50 MB** depending on runtime.

> [!TIP]
> **Run for ≥15 minutes.** At 250 Hz sampling, a 15-minute run collects ~225,000 samples — enough for statistically stable per-function attribution. Runs of **15–30 minutes** on a representative workload (e.g., a full chromosome or a panel of ~50 targeted regions) produce the most actionable profiles. Short runs (<5 min) have high variance and may miss functions that only dominate under sustained load.

---

## Analyzing a Profile

All analysis commands use `pixi run -e profiling analyze-profile`. The `-e profiling` flag is **required** — the script depends on `rich` and `jinja2`, which are only available in the profiling environment.

### Quick Start

```bash
# Full analysis (all views except source-line)
pixi run -e profiling analyze-profile <profile.bin>

# Single view
pixi run -e profiling analyze-profile <profile.bin> -- --view overview
```

### Available Views

| View | Flag | In `all`? | Description |
|:-----|:-----|:---------:|:------------|
| `overview` | `--view overview` | ✓ | Total CPU time, function counts, top-5 summary |
| `top` | `--view top` | ✓ | Ranked function tables by self-time and cumulative time |
| `modules` | `--view modules` | ✓ | Functions grouped by Lancet2 layer and external dependency |
| `components` | `--view components` | ✓ | High-level component attribution with internal-vs-external split |
| `tree` | `--view tree` | ✓ | Caller/callee tree from pprof |
| `hotpaths` | `--view hotpaths` | ✓ | Drill-down into top Lancet2 functions by cumulative time |
| `lines` | `--view lines` | ✗ | Source-line attribution (high-detail, excluded from `all`) |
| `all` | `--view all` | — | All of the above except `lines` (default) |

### Additional Flags

| Flag | Description |
|:-----|:------------|
| `--top N` | Number of top functions to show (default: **30**) |
| `--focus REGEX` | Focus the tree view on functions matching a regex |
| `--nodecount N` | Max nodes for pprof to consider (default: **200**) |
| `--binary PATH` | Path to the Lancet2 binary for symbolization (auto-detected) |

---

## Output Modes

### Terminal Output (CLI)

The default mode. All views render in the terminal using `rich` — tables, panels, trees, and color-coded severity indicators. Output is TTY-aware: colors are suppressed when piping to a file.

Every view listed above has CLI output.

### Static HTML Report

Generate a self-contained HTML report with `--html`:

```bash
pixi run -e profiling analyze-profile <profile.bin> -- --html report.html
```

The HTML report includes:
- **Overview** stat cards (total time, binary, functions shown)
- **Component breakdown** table with gradient bar charts
- **Top 50 functions** with sortable columns (click headers)
- **Zoomable SVG call graph** with ➕/➖/↺/⊞ controls and Ctrl+scroll

The report is a single self-contained HTML file (~60–100 KB) with no external dependencies. It uses a dark theme and can be shared directly or committed to a PR for review.

> [!NOTE]
> **HTML vs CLI coverage.** The HTML report contains overview, components, and top functions. Tree, hotpaths, modules, lines, and annotated source views are CLI-only. This is intentional — tree and hotpath views are interactive exploration tools best suited to the terminal.

### Annotated Source View

Inspect which source lines consume the most CPU for a specific function:

```bash
pixi run -e profiling analyze-profile <profile.bin> -- --list "BuildGraph"
```

This invokes `pprof -list=<regex>` and renders the annotated C++ source with syntax highlighting via `rich.syntax.Syntax`. Useful for identifying hot loops within a function.

---

## Diffing Two Profiles

### Same-Binary Diff (`--diff-base`)

Compare two profiles collected from the **same binary** (e.g., different inputs, different thread counts):

```bash
pixi run -e profiling analyze-profile <new_profile.bin> -- --diff-base <old_profile.bin>
```

The diff report shows:
1. **Overview panel** — base total, new total, net delta (seconds and percentage), and a verdict (REGRESSION ▲ or IMPROVEMENT ▼)
2. **Top changes table** — functions sorted by absolute self-time delta, with 🔴 regressed / 🟢 improved status
3. **Regression warnings** — auto-detected functions with >10% regression AND >1s absolute increase

> [!WARNING]
> **`--diff-base` only works when both profiles were collected with the same binary.**
> pprof symbolizes both profiles using one binary's symbol table. If the code changed between profiles (different inlining decisions, renamed functions, different function layout), the diff produces garbage — functions appear as "new" or "disappeared" when they actually exist in both.

### Cross-Binary Diff (`--diff-tag`)

Compare profiles collected from **different binaries** (e.g., before and after a code change). This is the primary mechanism for measuring optimization impact:

```bash
# Step 1: save each profile with its own binary at collection time
pixi run -e profiling analyze-profile <old.bin> -- --binary <old_binary> --save-summary before-change
pixi run -e profiling analyze-profile <new.bin> -- --binary <new_binary> --save-summary after-change

# Step 2: diff by tag — no binary needed, uses pre-symbolized data
pixi run -e profiling analyze-profile -- --diff-tag before-change after-change
```

The diff report shows:
1. **Overview panel** — base tag, new tag, version info, total time delta, and verdict
2. **Component-level changes** — per-component delta table (minimap2, graph, HTSlib, etc.)
3. **Function-level changes** — top 40 functions sorted by absolute self-time delta
4. **Summary panels** — top improvements (>1s reduction) and regressions (>1s increase)

> [!IMPORTANT]
> **`--save-summary` captures all symbolized function data permanently.** The binary is only needed once — at save time. After that, `--diff-tag` works by function name from the pre-symbolized snapshots, producing accurate diffs regardless of binary changes.

> [!NOTE]
> Old history entries saved before `--diff-tag` support only contain `top_10` functions. The diff will work but with limited coverage. Re-run `--save-summary` with the original binary to capture full function data.

---

## Performance History Tracking

The profiling toolkit maintains a JSONL-based history file at `profiling/history.jsonl` (git-tracked). Each entry captures a snapshot of a profile run with full metadata.

### Saving a Summary

After analyzing a profile, save a tagged summary:

```bash
pixi run -e profiling analyze-profile <profile.bin> -- --view overview --save-summary <tag>
```

The `<tag>` is a human-readable label for this data point (e.g., `first-baseline`, `post-ksw2-opt`, `v2.10.0-rc1`).

### What Gets Saved

Each history entry contains:

| Field | Description |
|:------|:------------|
| `tag` | User-provided label |
| `version.full` | Git version string: `VERSION_TAG-BRANCH-COMMIT[-dirty]` |
| `version.version_tag` | Semver tag (e.g., `v2.9.0`) |
| `version.branch` | Current git branch |
| `version.commit` | Short commit hash (10 chars) |
| `version.dirty` | Whether the working tree has uncommitted changes |
| `timestamp` | UTC ISO 8601 timestamp |
| `total_sec` | Total CPU time in seconds |
| `top_10` | Top 10 functions by self-time (backward compat) |
| `functions` | **All** profiled functions — the pre-symbolized snapshot used by `--diff-tag` |
| `components` | Per-component breakdown (flat_sec, flat_pct) |
| `modules` | Per-module breakdown for modules ≥0.1s |
| `internal_vs_external` | Lancet2 own code vs minimap2 vs HTSlib vs SPOA vs other |

### Viewing Trends

```bash
pixi run -e profiling analyze-profile -- --history
```

This renders a trend table showing total time and key component percentages across all saved profiles:

```
                      Profile History — Performance Trends
┏━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━┓
┃ Tag             ┃ Version                     ┃ Date       ┃   Total ┃minimap2 ┃   graph ┃genotype ┃    SPOA ┃
┡━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━┩
│ first-baseline  │ v2.9.0-main-358a2fcd15      │ 2026-04-16 │  89.5m  │  36.1%  │  19.2%  │   0.2%  │   0.8%  │
│ post-ksw2-opt   │ v2.9.1-main-7b3c4d5e6f      │ 2026-04-18 │  72.3m  │  28.4%  │  22.1%  │   0.3%  │   1.0%  │
└─────────────────┴─────────────────────────────┴────────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```

---

## Typical Workflow

A data-driven performance optimization loop:

1. **Baseline.** Build with profiling, run on a representative workload, save the summary:
   ```bash
   pixi run -e profiling analyze-profile <profile.bin> -- --save-summary first-baseline
   ```

2. **Identify bottleneck.** Use `--view components` to find which component dominates, then `--view top` and `--view modules` to find specific functions. Use `--list "FunctionName"` for source-line detail.

3. **Optimize.** Make your code change.

4. **Measure.** Rebuild with profiling, re-run the same workload, save with a new tag:
   ```bash
   pixi run -e profiling analyze-profile <new.bin> -- --binary <new_binary> --save-summary post-my-optimization
   ```

5. **Diff.** Compare against the baseline by tag:
   ```bash
   pixi run -e profiling analyze-profile -- --diff-tag first-baseline post-my-optimization
   ```

6. **Repeat.** Check `--history` to verify the trend is going in the right direction.

---

## Component Classification

The profiling toolkit classifies every profiled function into a **module** (fine-grained) and a **component** (coarse-grained). Classification is rule-based — the first matching prefix wins.

### Modules

| Prefix | Module |
|:-------|:-------|
| `lancet::base::` | `lancet/base` |
| `lancet::hts::` | `lancet/hts` |
| `lancet::cbdg::` | `lancet/cbdg` |
| `lancet::caller::` | `lancet/caller` |
| `lancet::core::` | `lancet/core` |
| `lancet::cli::` | `lancet/cli` |
| `spoa::` | `spoa` |
| `ksw_*` | `minimap2/ksw2` |
| `mm_*`, `mg_*` | `minimap2` |
| `bgzf_*`, `hts_*`, `sam_*` | `htslib` |
| `deflate_*`, `inflate_*` | `htslib/zlib` |
| `absl::` | `abseil` |
| `mi_*`, `_mi_*` | `mimalloc` |
| `std::` | `libstdc++` |
| `__memcpy`, `__read`, etc. | `sys/*` |

### Components

| Component | Modules |
|:----------|:--------|
| **Lancet2 (graph assembly)** | `lancet/base`, `lancet/cbdg` |
| **Lancet2 (genotyping)** | `lancet/caller` |
| **Lancet2 (pipeline)** | `lancet/core`, `lancet/cli` |
| **Lancet2 (HTS I/O)** | `lancet/hts` |
| **minimap2** | `minimap2`, `minimap2/ksw2` |
| **SPOA (MSA)** | `spoa` |
| **HTSlib** | `htslib`, `htslib/zlib` |
| **Memory management** | `mimalloc`, `allocator` |
| **System** | all `sys/*` modules |

To add a new library or Lancet2 module, edit the `_MODULE_RULES` and `_COMPONENT_MAP` tables in `scripts/analyze_profile.py`.

---

## File Reference

| File | Purpose |
|:-----|:--------|
| `scripts/analyze_profile.py` | Main profiling analysis script (all modes) |
| `scripts/profile_report.html.j2` | Jinja2 template for static HTML reports |
| `profiling/history.jsonl` | Git-tracked performance history (JSONL) |
| `pixi.toml` (`[feature.profiling]`) | Profiling environment definition (`rich`, `jinja2`) |

---

## Command Reference

```bash
# ── Single profile analysis ──────────────────────────────────────────────
pixi run -e profiling analyze-profile <profile.bin>                              # all views
pixi run -e profiling analyze-profile <profile.bin> -- --view overview           # single view
pixi run -e profiling analyze-profile <profile.bin> -- --view top --top 50       # top 50
pixi run -e profiling analyze-profile <profile.bin> -- --view lines             # source-line
pixi run -e profiling analyze-profile <profile.bin> -- --list "BuildGraph"       # annotated source
pixi run -e profiling analyze-profile <profile.bin> -- --html report.html        # HTML report

# ── Diff mode ────────────────────────────────────────────────────────────
pixi run -e profiling analyze-profile <new.bin> -- --diff-base <old.bin>          # same-binary diff
pixi run -e profiling analyze-profile <new.bin> -- --diff-base <old.bin> --html d.html  # HTML diff
pixi run -e profiling analyze-profile -- --diff-tag <base-tag> <new-tag>           # cross-binary diff

# ── History tracking ─────────────────────────────────────────────────────
pixi run -e profiling analyze-profile <profile.bin> -- --save-summary <tag>       # save snapshot
pixi run -e profiling analyze-profile -- --history                                # view trends
```
