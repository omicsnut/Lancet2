# Probe Tracking Infrastructure

Developer documentation for the Lancet2 probe variant forensic pipeline.
This system answers: **for every missed truth variant, exactly where in the
pipeline did it get lost?**

The architecture has two halves:

1. **C++ data emitter** — records raw diagnostic facts per
   `(probe_id, window, comp_id, k)` attempt with no attribution logic.
2. **Python attribution engine** (`scripts/analyze_probe_results.py`) —
   derives `lost_at_stage` from those facts using a bottom-up cascade.

---

## Phase 0: CLI Setup

**File:** `src/lancet/cli/cli_interface.cpp`

```
--probe-variants <missed_variants.txt>   (input: truth variants to track)
--probe-results  <probe_results.tsv>     (output: raw fact table)
```

Both flags are **bidirectionally required** (`probe_variants_opt->needs(probe_results_opt)` and vice versa). If either is omitted, the entire probe system is inactive — zero overhead in production.

---

## Phase 1: Global Initialization (once, main thread)

**File:** `src/lancet/cli/pipeline_runner.cpp` → `SetupProbeTracking()`

```
SetupProbeTracking():
  1. Load variants from missed_variants.txt → vector<ProbeVariant>
  2. Build ProbeIndex (precomputed k-mer hashes for all k in [min_k, max_k, step])
  3. Create shared ProbeResultsWriter (thread-safe, owns the variant list)
  4. Store both in VariantBuilder::Params (shared_ptr, immutable after this point)
```

**ProbeIndex**: for each truth variant × each k-value, precomputes ALT-unique
k-mer hashes. This avoids recomputing them per-window. Shared read-only across
all worker threads.

**ProbeResultsWriter**: owns the output file and the `mWrittenProbeIds` set
(tracks which probe_ids have been written). Thread-safe via `absl::Mutex`.

---

## Phase 2: Per-Thread Wiring (one VariantBuilder per worker)

**File:** `src/lancet/core/variant_builder.cpp` → constructor

```
VariantBuilder(params):
  mProbeDiagnostics.Initialize(
      params->mProbeVariantsPath,      // loads variants into ProbeTracker
      params->mProbeResultsWriter,     // shared writer (thread-safe)
      params->mProbeIndex              // shared index (read-only)
  )
  mDebruijnGraph.SetProbeTracker(mProbeDiagnostics.Tracker())
```

### Ownership chain

```
VariantBuilder (per-thread)
  └── ProbeDiagnostics (owns ProbeTracker)
        └── ProbeTracker
              ├── mRecords: vector<ProbeKRecord>  (per-window accumulator)
              ├── mNodeTags: NodeTagMap            (k-mer → graph node mapping)
              ├── mResultsWriter: shared_ptr       (→ global writer)
              └── mProbeIndex: shared_ptr           (→ global index)

Graph (per-thread)
  └── mProbeTrackerPtr: ProbeTracker*  (non-owning, nullable)
```

When `--probe-variants` is not set: `Tracker()` returns `nullptr`,
`mProbeTrackerPtr` is null, every `Probe*` method in Graph early-returns.
**Zero overhead**.

---

## Phase 3: Per-Window Graph Processing

**File:** `src/lancet/cbdg/graph.cpp` → `BuildComponentResults()`

### 3a. Context Construction

```cpp
Context probe_ctx{
    .mChrom = region_chrom,
    .mRefSeq = region_seq,
    .mRegStr = region_str,        // "chr1:1000-2000" (the window)
    .mRegionStart = region_start0
};
```

This context flows through **every** probe operation in the window, ensuring
all records carry the window identity.

### 3b. K-value Loop

For each `k` in `[min_k, min_k+step, ..., max_k]`:

```
probe_ctx.mKmerSize = k
probe_ctx.mCompId = 0  (reset before component loop)

┌─ BuildGraph(mate_mers)
│
├─ ProbeGenerateAndTag(probe_ctx)     ← tags ALT-unique k-mers in graph nodes
├─ ProbeCountInReads(probe_ctx)       ← counts how many reads carry each ALT k-mer
├─ ProbeLogStatus(BUILD, probe_ctx)   ← snapshot survival after graph build
│
├─ RemoveLowCovNodes(0)
├─ ProbeLogStatus(LOWCOV1, probe_ctx) ← snapshot after first low-cov removal
│
├─ MarkConnectedComponents
│
└─ For each component:
     probe_ctx.mCompId = component_index
     │
     ├─ FindSource/FindSink
     │   └─ ProbeSetNoAnchor(probe_ctx)     ← if no anchor found
     │   └─ ProbeSetShortAnchor(probe_ctx)  ← if anchor too short
     │
     ├─ ProbeCheckAnchorOverlap(source, sink, probe_ctx)  ← variant_in_anchor
     │
     ├─ PruneComponent(component_index)
     │   ├─ ProbeLogStatus(COMPRESS1, prune_ctx)  ← 4 more survival snapshots
     │   ├─ ProbeLogStatus(LOWCOV2, prune_ctx)
     │   ├─ ProbeLogStatus(COMPRESS2, prune_ctx)
     │   └─ ProbeLogStatus(TIPS, prune_ctx)
     │
     ├─ HasCycle?
     │   └─ ProbeSetGraphCycle(probe_ctx) → break (retry at higher k)
     │
     ├─ IsComplex?
     │   └─ ProbeSetGraphComplex(probe_ctx) → break (retry at higher k)
     │
     ├─ BuildHaplotypes(comp_id, ..., probe_ctx)
     │   └─ ProbeSetTraversalLimit(probe_ctx)  ← if BFS budget exhausted
     │
     └─ ProbeCheckPaths(haplotypes, probe_ctx) ← records which hap indices hit
```

### 3c. Key Invariant: 4-Key Records

Every `Probe*` method calls `FindOrCreateRecord` with the full
`(probe_id, window, comp_id, k)` key. If a record doesn't exist yet, it's
created. If it does, the existing record is updated.

**Pre-component stages** (BUILD, LOWCOV1): `comp_id = 0`, survival counts span
**all** nodes. One record per `(probe, window, k=0)`.

**Post-component stages** (COMPRESS1 through TIPS): `comp_id = N`, survival
counts filtered to **only** nodes belonging to component N via
`CountSurvivingKmers(nodes, comp_id)`.

### 3d. Structural Flag Filtering

`SetNoAnchor` and `SetShortAnchor` check
`HasTaggedNodesInComponent(probe_id, comp_id, nodes)` before setting the flag.
Without this, a structural failure in component 3 would falsely flag probes
that only exist in component 1.

---

## Phase 4: MSA + Genotyper Annotation

**File:** `src/lancet/core/variant_builder.cpp`

```
For each component:
  extracted = ExtractVariants(component, component_idx, window)
  if (extracted.IsEmpty()) continue;
  // Probes with paths in skipped components keep all MSA flags false.
  // Python classifies these as msa_not_extracted.

  CheckMsaExtraction(extracted, window)   ← annotates ALL records with paths
  Genotyper::Genotype(...)
  CheckGenotyperResult(result, extracted) ← annotates ALL records with paths

SubmitCompleted()  ← flushes records to shared writer
```

### CheckMsaExtraction

**File:** `src/lancet/core/probe_diagnostics.cpp`

Iterates `FindRecordsWithPaths(probe_id)` — returns **every** record that has
non-empty `mHapIndices`. Each record gets independent classification:

| Tier | Condition | Field Set |
|------|-----------|-----------|
| 1 | Same position + same alleles | `mIsMsaExactMatch = true` |
| 2 | Same alleles, different position | `mIsMsaShifted = true`, `mMsaShiftBp` |
| 3 | Position matches, allele absorbed into MNV | `mIsMsaSubsumed = true` |
| — | No match at any tier | All flags remain false |

### CheckGenotyperResult

Same iteration over `FindRecordsWithPaths`. Sets:

- `mGenoTrueAltReads`, `mGenoTotalRefReads`
- `mGenoStolenToRef`, `mGenoStolenToWrongAlt`
- `mGenoNonOverlapping`
- `mIsGenoHasAltSupport` — true if ALT reads > 0
- `mIsGenoNoOverlap` — true if zero read alignments overlapped the variant

---

## Phase 5: TSV Emission

### Per-Window (SubmitCompleted)

**File:** `src/lancet/cbdg/probe_tracker.cpp`

```cpp
void ProbeTracker::SubmitCompleted() {
  if (mRecords.empty() || !mResultsWriter) return;
  mResultsWriter->Append(mRecords);  // thread-safe (mutex inside)
  mRecords.clear();
}
```

### Post-Pipeline (EmitUnprocessedProbes)

**File:** `src/lancet/cli/pipeline_runner.cpp`

Called after `executor.Execute()` completes. Emits a synthetic
`ProbeKRecord{.mProbeId = idx, .mIsNotProcessed = true}` for every variant
that was never activated (fell in N-only, repeat, or inactive windows).
Ensures **every input variant has at least one output row**.

### TSV Schema (40 columns)

**File:** `src/lancet/cbdg/probe_results_writer.cpp`

Pure fact table. No `lost_at_stage` column. No attribution. Every boolean is
emitted as `0`/`1`.

---

## Phase 6: Python Attribution Engine

**File:** `scripts/analyze_probe_results.py`

### Step 1: Load raw TSV

```python
probes = load_probe_results(args.probe_results)  # multi-row per probe
```

### Step 2: Derive lost_at_stage — bottom-up cascade per row

```python
probes = derive_lost_at(probes)
```

Checks deepest evidence first. First match wins:

| Priority | Condition | Attribution |
|----------|-----------|-------------|
| 1 (deepest) | `is_geno_has_alt_support == 1` | `survived` |
| 2 | `is_msa_exact_match` and `is_geno_no_overlap` | `geno_no_overlap` |
| 3 | `is_msa_exact_match` and stolen > 0 | `geno_stolen` |
| 4 | `is_msa_exact_match` and no ALT support | `geno_zero_alt_reads` |
| 5 | `hap_indices != "."` and `is_msa_shifted` | `msa_shifted` |
| 6 | `hap_indices != "."` and `is_msa_subsumed` | `msa_subsumed` |
| 7 | `hap_indices != "."` (no MSA match) | `msa_not_extracted` |
| 8 | `n_surviving_tips > 0` and `traversal_limited` | `bfs_exhausted` |
| 9 | `n_surviving_tips > 0` | `no_path` |
| 10–15 | Pruning stages (reverse pipeline order) | `pruned_at_*` |
| 16 | `is_graph_cycle` | `graph_has_cycle` |
| 17 | `is_graph_complex` | `graph_too_complex` |
| 18 | `is_no_anchor` | `no_anchor` |
| 19 | `is_short_anchor` | `short_anchor` |
| 20 | `is_variant_in_anchor` | `variant_in_anchor` |
| 21 | fallback | `not_processed` |

### Step 3: Select best record per probe

```python
probes = compute_depth(probes)               # integer score 0–20
attribution = select_best_per_probe(probes)  # deepest stage, then highest k
```

### Step 4: Report sections

| Section | Description |
|---------|-------------|
| §1 Scorecard | Coverage validation, vital signs |
| §2 Funnel | Stage attribution distribution (21-level cascade) |
| §3 Survival | K-mer attrition through 6 pruning stages |
| §4 Breakdown | Type × Size × Stage cross-tabulation |
| §5 Genotyper | MSA/Genotyper read assignment forensics |
| §6 Targets | Top 30 variants closest to success |
| §7 Deep Dive | Detailed analysis of top 2 loss stages |
| §8 Windows | Cross-window attribution, boundary sensitivity |

---

## Gap Analysis

### Every pipeline exit point is covered

| Exit Point | Recording Mechanism | Python Attribution |
|------------|--------------------|--------------------|
| Window skipped (N-only, repeat, inactive) | `EmitUnprocessedProbes` | `not_processed` |
| No anchor found | `ProbeSetNoAnchor` | `no_anchor` |
| Anchor too short | `ProbeSetShortAnchor` | `short_anchor` |
| Graph has cycle | `ProbeSetGraphCycle` | `graph_has_cycle` |
| Graph too complex | `ProbeSetGraphComplex` | `graph_too_complex` |
| K-mers pruned at any stage | `ProbeLogStatus` (6 stages) | `pruned_at_*` |
| K-mers survive but no path carries ALT | `ProbeCheckPaths` (empty hap_indices) | `no_path` |
| Paths found, MSA extracts nothing | `extracted.IsEmpty()` continue | `msa_not_extracted` |
| MSA extracts but variant shifted | `CheckMsaExtraction` tier 2 | `msa_shifted` |
| MSA extracts but variant subsumed | `CheckMsaExtraction` tier 3 | `msa_subsumed` |
| MSA exact match, zero read overlap | `CheckGenotyperResult` | `geno_no_overlap` |
| MSA exact match, reads stolen | `CheckGenotyperResult` | `geno_stolen` |
| MSA exact match, zero ALT reads | `CheckGenotyperResult` | `geno_zero_alt_reads` |
| MSA exact match, ALT reads present | `CheckGenotyperResult` | `survived` |

### Diagnostic annotations (not exit points)

These flags are set but do **not** terminate processing — the component
continues through pruning, paths, MSA, and genotyping normally:

| Annotation | Recording Mechanism | Effect in Python |
|------------|--------------------|--------------------|
| Variant falls inside anchor sequence | `ProbeCheckAnchorOverlap` | `variant_in_anchor` (standalone attribution at priority 20) |
| BFS walk-tree budget exhausted | `ProbeSetTraversalLimit` | Refines `no_path` → `bfs_exhausted` |

**`variant_in_anchor`**: the truth variant's genomic position overlaps the
invariant source/sink anchor nodes. The graph cannot distinguish this variant
from reference by construction, but the pipeline still attempts full assembly.
In Python's cascade it sits at priority 20 (lowest), so it only becomes the
final attribution when nothing deeper succeeded.

**`is_traversal_limited`**: set inside `BuildHaplotypes` when the
BFS budget is hit, but the function still returns whatever (partial) haplotypes
were found. Processing continues into `ProbeCheckPaths`, MSA, and genotyping
with those partial results. In Python's cascade, it distinguishes two
`no_path` scenarios:

- `bfs_exhausted`: k-mers survived pruning, no ALT path found, **but** BFS
  budget was hit — the ALT path may exist but was not enumerated.
- `no_path`: k-mers survived pruning, no ALT path found, BFS completed
  fully — the ALT path genuinely does not exist in the graph.

### No diagnostic data is lost

- Every `(probe_id, window, comp_id, k)` attempt creates its own record.
- `FindRecordsWithPaths` returns **all** records (not just "best"), so every
  attempt gets MSA/Genotyper annotation.
- Component-filtered survival counting prevents cross-component noise.
- Structural flags are only set on probes with tagged nodes **in that
  component**.

### Thread safety

- `ProbeTracker` is per-thread (owned by each `VariantBuilder`).
- `ProbeResultsWriter::Append` acquires `absl::Mutex` before writing.
- `ProbeIndex` is `shared_ptr<const>` — immutable, lock-free reads.
