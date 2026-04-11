# Active Region Detection

Before assembling the de Bruijn graph for a window, Lancet2 runs a fast pre-scan to check whether the region contains any evidence of variation. Windows with no mutation signal are skipped entirely, reducing WGS runtime by ~80%.

## How It Works

The active region detector iterates all reads overlapping the window using zero-copy alignment proxies (no sequence or quality deep-copy). For each read, it inspects three evidence sources:

1. **MD tag mismatches** — If the BAM contains `MD` tags, each mismatch position (with base quality ≥ Q20) is tracked. If **≥2 reads** show a mismatch at the **same** genomic position, the window is active.

2. **CIGAR-based events** — Insertions, deletions, and explicit sequence mismatches (`X` ops) are counted at each genomic position. Again, if **≥2 independent reads** report the same event type at the same position, the window is active.

3. **Soft clips** — If a read contains soft-clipped bases, the boundary positions are tracked. **≥2 soft clips** at the same position trigger an active call.

The detector terminates immediately on the first evidence hit (`≥2` at any position), so the cost is **O(R)** in the best case and **O(R × C)** in the worst case, where R = number of overlapping reads and C = mean CIGAR length.

## The MD Tag Requirement

!!! warning "BAMs without MD tags cause ~5–10× slower runtime"
    If the input BAM/CRAM files lack the `MD` auxiliary tag, Lancet2 **automatically disables active region detection** and assembles every window. This is because the MD tag is the only source of per-base mismatch information when the CIGAR uses generic `M` operations (which do not distinguish matches from mismatches).

    You will see a runtime log message:

    ```
    MD tag not found in input BAM. Disabling active region detection.
    ```

    To restore fast-skip behavior, generate MD tags before running Lancet2:

    ```bash
    samtools calmd -b input.bam reference.fasta > input_with_md.bam
    samtools index input_with_md.bam
    ```

## `--no-active-region`

This flag forces Lancet2 to assemble **every** window regardless of mutation evidence. Use cases:

- **Completeness auditing** — ensure no variant is missed by the heuristic at the cost of runtime.
- **Known variant regions** — if you know a window contains a variant that the heuristic might miss (e.g., a mosaic variant with only 1 supporting read per position), forcing assembly guarantees it will be evaluated.

Enabling this flag has the same runtime impact as BAMs lacking MD tags — expect **5–10× longer** wall-clock time on WGS runs.

* **CLI reference:** [`--no-active-region`](../reference.md#flags)
