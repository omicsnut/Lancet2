# Read Filtering & Downsampling

This guide documents which reads Lancet2 includes in analysis, how per-sample downsampling works, and how the `--extract-pairs` flag modifies read collection behavior.

---

## Read Exclusion Criteria

The following reads are silently excluded during both active region detection and read collection:

| Flag | SAM Bit | Effect |
|:-----|:--------|:-------|
| QC-failed | `0x200` | Excluded — sequencer-flagged quality failures |
| Duplicate | `0x400` | Excluded — optical/PCR duplicate reads |
| Unmapped | `0x4` | Excluded — reads with no alignment |
| MAPQ = 0 | — | Excluded — multi-mapped reads with ambiguous placement |

Additionally, during the **Profile pass** (Pass 1), reads with MAPQ < 20 are counted but do not contribute to the downsampling candidate pool.

---

## Two-Pass Read Collection

Read collection uses a memory-efficient two-pass strategy per sample:

### Pass 1: Profile & Downsample Math

Iterates the region using **zero-copy alignment proxies** — no sequence or quality data is deep-copied. This pass:

1. Counts total passing reads and bases.
2. Collects unique qname hashes for all passing reads.
3. Computes a coverage-based downsampling threshold from `--max-sample-cov`.
4. If `--extract-pairs` is enabled, tracks out-of-region mate locations.

The downsampling threshold is computed as:

```
max_reads = ceil(max_sample_cov × region_length / mean_read_length)
```

If the number of passing reads exceeds `max_reads`, the qname hash list is shuffled with a **deterministic seed** and truncated. Both mates of a pair are symmetrically accepted or rejected — if one mate is in the keep set, its partner is always included.

!!! note "Deterministic Downsampling"
    The downsampling shuffle uses a fixed seed (`std::default_random_engine(0)`) to guarantee **identical variant calls** between runs on the same input. This is a deliberate design choice for reproducibility in clinical and production pipelines.

### Pass 2: Deep Copy & Object Construction

Re-iterates the region. Only reads whose qname hash is in the keep set trigger expensive operations — `BuildSequence()` and `BuildQualities()` via the `Read` constructor. Reads not in the keep set are skipped entirely.

### Pass 3: Mate Recapture

If `--extract-pairs` is enabled, mates that fall **outside** the current window region are retrieved via targeted single-position queries using `htslib`. Mates whose qnames were not kept during downsampling are also excluded.

---

## `--extract-pairs`

When this flag is enabled, the read collector changes behavior:

1. **SA tag extraction**: The supplementary alignment tag (`SA`) is additionally requested from the BAM alongside `AS` and `XS`. Without this flag, only `AS` and `XS` are extracted.

2. **Mate tracking**: For each read that is **not** in a proper pair (`FLAG & 0x2 == 0`) **or** has a supplementary alignment (`SA` tag present), the mate's location is recorded for out-of-region retrieval in Pass 3. Reads where both mates are already in-region are excluded from mate retrieval since both copies are already captured.

3. **Runtime impact**: Mate retrieval issues additional targeted `htslib` queries (one per out-of-region mate location), which can slow down read collection significantly — particularly in regions with many discordant pairs or structural variant breakpoints.

**When to enable:** Use `--extract-pairs` when you expect structural variant breakpoints at window edges where one mate maps inside and the other maps far away. For standard SNV/indel calling, this flag is not needed.

* **CLI reference:** [`--extract-pairs`](../reference.md#flags)

---

## `--max-sample-cov`

Controls the maximum per-sample read coverage retained per window.

- **Default:** 1000×
- **Effect:** Windows with >1000× per-sample coverage are downsampled to ~1000× using the deterministic paired downsampling described above.
- **Trade-off:** Lower values = faster runtime and lower memory, but reduced sensitivity for subclonal variants. Higher values = more reads kept for assembly, but diminishing returns above ~500× for most variant types.

* **CLI reference:** [`--max-sample-cov`](../reference.md#parameters)
