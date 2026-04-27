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

## Phase 4: Variant Extraction + Genotyper Annotation

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
- `mGenoReassignedToRef`, `mGenoReassignedToWrongAlt` — ALT-carrying reads
  that the genotyper's alignment scoring assigned to REF or a different ALT
  haplotype instead of the truth ALT allele
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

Checks deepest evidence first. First match wins. The 27-level cascade spans
7 categories:

| Category | Stages | Description |
|----------|--------|-------------|
| Genotyper | `survived`, `geno_reads_reassigned`, `geno_zero_alt_reads`, `geno_no_overlap` | Deepest — variant reached genotyping |
| Variant extraction | `msa_not_extracted`, `msa_shifted`, `msa_subsumed` | Haplotype paths found, MSA outcome |
| Path finding | `bfs_exhausted`, `no_path` | K-mers survived pruning, path enumeration outcome |
| Pruning | `pruned_at_tips`, `pruned_at_compress2`, `pruned_at_lowcov2`, `pruned_at_compress1`, `pruned_at_lowcov1`, `pruned_at_build` | 6 graph simplification stages |
| Graph construction | `graph_has_cycle`, `graph_too_complex`, `no_anchor`, `short_anchor`, `variant_in_anchor` | Structural graph failures |
| Not processed | `not_processed` (+ 6 sub-stages with `--log`) | Variant never entered the graph pipeline |
| Survived | `survived` | Variant reached genotyper with ALT support |

### Step 3: Select best record per probe

```python
probes = compute_depth(probes)               # integer score 0–26
attribution = select_best_per_probe(probes)  # deepest stage, then highest k
```

### Step 4: Sub-classify not_processed (optional, requires --log)

When `--log` points to the Lancet2 debug log, probes attributed as
`not_processed` are sub-classified by their window's completion status:

| Sub-stage | Window Status | Meaning |
|-----------|--------------|---------|
| `not_processed:ref_all_n` | `SKIPPED_NONLY_REF_BASES` | Window reference is all N bases |
| `not_processed:ref_repeat` | `SKIPPED_REF_REPEAT_SEEN` | Window has k-mer repeats (guaranteed graph cycle) |
| `not_processed:inactive` | `SKIPPED_INACTIVE_REGION` | Window had no mutation evidence |
| `not_processed:low_coverage` | `SKIPPED_LOW_COVERAGE` | Window coverage below MinAnchorCov (5×) |
| `not_processed:no_alt_haplotype` | `SKIPPED_NOASM_HAPLOTYPE` | Assembly ran but found no variant haplotypes |
| `not_processed:other_variant_called` | `FOUND_GENOTYPED_VARIANT` | Window called other variants, not this one |

When multiple windows overlap a probe's position (250bp overlap between
adjacent 500bp windows), the window with the deepest pipeline stage wins.
On ties, the window where the variant is most centered wins.

### Step 5: Report sections

| Section | View name | Description |
|---------|-----------|-------------|
| §1 Scorecard | `scorecard` | Coverage validation, vital signs |
| §2 Funnel | `funnel` | Stage attribution distribution (27-level cascade) |
| §3 Survival | `survival` | K-mer attrition through 6 pruning stages |
| §4 Breakdown | `breakdown` | Type × Size × Stage cross-tabulation |
| §5 Genotyper | `genotyper` | Variant extraction and read assignment forensics |
| §6 Targets | `targets` | Top 30 variants closest to success |
| §7 Deep Dive | `deepdive` | Detailed analysis of top 2 loss stages |
| §8 Windows | `windows` | Cross-window attribution, boundary effects |

---

## Gap Analysis

### Every pipeline exit point is covered

| Exit Point | Recording Mechanism | Python Attribution |
|------------|--------------------|--------------------|
| Window skipped (N-only, repeat, inactive) | `EmitUnprocessedProbes` | `not_processed` (sub-classified with `--log`) |
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
| MSA exact match, reads reassigned | `CountReassignedReads` | `geno_reads_reassigned` |
| MSA exact match, zero ALT reads | `CheckGenotyperResult` | `geno_zero_alt_reads` |
| MSA exact match, ALT reads present | `CheckGenotyperResult` | `survived` |

### Diagnostic annotations (not exit points)

These flags are set but do **not** terminate processing — the component
continues through pruning, paths, MSA, and genotyping normally:

| Annotation | Recording Mechanism | Effect in Python |
|------------|--------------------|--------------------|
| Variant falls inside anchor sequence | `ProbeCheckAnchorOverlap` | `variant_in_anchor` (lowest priority — only attributed when nothing deeper succeeded) |
| BFS walk-tree budget exhausted | `ProbeSetTraversalLimit` | Refines `no_path` → `bfs_exhausted` |

**`variant_in_anchor`**: the truth variant's genomic position overlaps the
invariant source/sink anchor nodes. The graph cannot distinguish this variant
from reference by construction, but the pipeline still attempts full assembly.
In Python's cascade it sits at the lowest priority, so it only becomes the
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

---

## Running the Analysis

Both scripts require the `hts-tools` pixi environment (provides Polars,
Rich, pysam, edlib, tqdm).

The end-to-end workflow has three steps:

1. **Truth concordance** — identify which truth variants Lancet2 missed
2. **Lancet2 probe run** — re-run Lancet2 with probe tracking on the missed set
3. **Probe analysis** — attribute each missed variant to a pipeline stage

### Step 1: Truth Concordance

`scripts/truth_concordance.py` compares a truth VCF against a Lancet2
output VCF. For each truth variant, it assigns a concordance level
(L0, LD, L1–L3, MISS) and collects per-filter read evidence for missed
variants.

```bash
pixi run -e hts-tools python3 scripts/truth_concordance.py \
    --truth-small expected_snvs_indels.vcf.gz \
    --lancet lancet_output.vcf.gz \
    --ref reference.fa \
    --samples normal.bam tumor.bam \
    --output-dir data/ \
    --mode small
```

| Flag | Required | Description |
|------|----------|-------------|
| `--truth-small` | one of small/large | Small variant truth VCF (GIAB SNVs + indels) |
| `--truth-large` | one of small/large | Large variant truth VCF (Manta PASS INS/DEL) |
| `--lancet` | yes | Lancet2 output VCF |
| `--ref` | yes | Reference FASTA (required for CRAM and edit distance) |
| `--samples` | yes | One or more BAM/CRAM alignment files |
| `--mode` | no | `small`, `large`, or `all` (default: `all`) |
| `--workers` | no | Parallel workers for read evidence (default: 16) |
| `--output-dir` | no | Output directory (default: `data/`) |

Key outputs:

- `missed_variants.txt` — missed variants with per-filter read support counts
  (**required input for step 2**)
- `concordance_details.txt` — all truth variants with match levels
  (optional enrichment for probe analysis §1)
- `truth_concordance_report.txt` — rich diagnostic report (§1–§4)

### Step 2: Lancet2 Probe Run

Re-run Lancet2 with the missed variants file to collect forensic data:

```bash
Lancet2 pipeline \
    -n normal.bam -t tumor.bam -r reference.fa -o out.vcf.gz \
    --probe-variants data/missed_variants.txt \
    --probe-results data/probe_results.tsv
```

Both `--probe-variants` and `--probe-results` must be provided together.
Enable `--verbose` logging if you plan to use the `--log` flag in step 3
for not_processed sub-classification.

### Step 3: Probe Analysis

```bash
pixi run -e hts-tools python3 scripts/analyze_probe_results.py \
    --probe-results data/probe_results.tsv \
    --missed-variants data/missed_variants.txt \
    --concordance-details data/concordance_details.txt \
    --log data/lancet2_debug.log \
    --output-dir data/ \
    --view all
```

| Flag | Required | Description |
|------|----------|-------------|
| `--probe-results` | yes | Raw TSV from Lancet2 `--probe-results` |
| `--missed-variants` | yes | From step 1 (`missed_variants.txt`) |
| `--concordance-details` | no | From step 1 (enriches §1 scorecard) |
| `--log` | no | Lancet2 debug log (enables not_processed sub-classification) |
| `--output-dir` | no | Output directory (default: `data/`) |
| `--view` | no | Single section or `all` (default: `all`) |

Output files:

| File | Format | Content |
|------|--------|---------|
| `probe_analysis_report.txt` | Rich text | Full diagnostic report (§1–§8) |
| `probe_stage_attribution.txt` | TSV | One row per probe: final `lost_at_stage`, type, size |
| `probe_survival_matrix.txt` | TSV | One row per `(probe, window, comp, k)`: 6 survival counts |

### Interpreting Results

Start with **§1 Scorecard** — verify that the coverage gap is zero (every
input variant produced at least one output row). A non-zero gap indicates
probes in regions Lancet2 never processed.

**§2 Funnel** shows the stage attribution distribution. The dominant loss
stage identifies the primary pipeline bottleneck. With `--log`, the
`not_processed` category expands into 6 sub-stages showing exactly why
each window was skipped.

**§6 Targets** lists the 30 variants closest to success (highest depth
score) — the highest-value debugging targets.

Use `--view` to render a single section during iterative debugging:

```bash
pixi run -e hts-tools python3 scripts/analyze_probe_results.py \
    --probe-results data/probe_results.tsv \
    --missed-variants data/missed_variants.txt \
    --view funnel
```
