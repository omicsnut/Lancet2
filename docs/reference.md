---
hide:
  - navigation
---

# Reference

## pipeline
The pipeline subcommand runs the entire Lancet variant calling pipeline on one (or) more region(s) of interest.
The full help text for the subcommand can be generated using the following command line

```bash
Lancet2 pipeline --help
```

### Required

#### `-r`,`--reference`
> [PATH]

Path to the reference FASTA file. Supports local paths and cloud URIs (`s3://`, `gs://`, `http(s)://`, `ftp(s)://`) when built with cloud I/O support.
See [Native Cloud Streaming](guides/cloud_streaming.md) for setup instructions.

#### `-o`,`--out-vcfgz`
> [PATH]

Output path to the compressed VCF file (`.vcf.gz`). Supports cloud URIs for direct streaming uploads.
When writing to a cloud bucket, Lancet2 performs an [upfront authentication check](guides/cloud_streaming.md#cloud-authentication-pre-validation) before processing any windows.
See [VCF Output Reference](guides/vcf_output.md) for the full output format specification.

### Datasets

#### `-n`,`--normal`
> [PATH...]

Path to one (or) more normal BAM/CRAM file(s). Required.
Multiple paths enable multi-sample mode where all normal samples share the `NORMAL` graph color. Sample names are read from BAM `SM` read group tags. See [Multi-Sample & Germline Mode](guides/architecture.md#multi-sample-germline-mode) for caveats.

#### `-t`,`--tumor`
> [PATH...]

Path to one (or) more tumor BAM/CRAM file(s). Optional — when omitted, Lancet2 runs in germline-only mode (no `SHARED`/`NORMAL`/`TUMOR` INFO tags in the VCF).
Multiple paths enable multi-sample somatic mode. See [Multi-Sample & Germline Mode](guides/architecture.md#multi-sample-germline-mode).

### Regions

#### `-R`,`--region`
> [REF:[:START[-END]]...]

One (or) more genomic regions (1-based, both inclusive). When neither `--region` nor `--bed-file` is specified, Lancet2 processes all non-decoy, non-mitochondrial reference contigs.

#### `-b`,`--bed-file`
> [PATH]

Path to BED file with regions to process. Supports cloud URIs (`s3://`, `gs://`, `http(s)://`) via `htslib`'s remote file streaming. The BED must have exactly 3 tab-separated columns (chrom, start, end). Comment lines starting with `#` are ignored.

#### `-P`,`--padding`
> [0-1000]. Default value --> 500

Extends each input region by N bases on both sides before windowing. Captures indel breakpoints that straddle region boundaries — a 50 bp deletion centered at a BED endpoint would be missed without sufficient padding.
See [Windowing & Overlap](guides/architecture.md#8-windowing-overlap) for how padding interacts with window size and overlap.

#### `-p`,`--pct-overlap`
> [10-90]. Default value --> 20

Percent overlap between consecutive windows. At defaults (1000 bp windows, 20% overlap), the step size is 800 bp, creating 200 bp overlaps. Higher overlap → fewer edge-effect blind spots but more total windows to process.
See [Windowing & Overlap](guides/architecture.md#8-windowing-overlap) for the step size formula and trade-offs.

#### `-w`,`--window-size`
> [1000-2500]. Default value --> 1000

Width of each micro-assembly window in base pairs. Larger windows provide more flanking context for large indels but consume more memory per thread.
See [Windowing & Overlap](guides/architecture.md#8-windowing-overlap) for detailed trade-offs.

### Parameters

#### `-T`,`--num-threads`
Number of async worker threads for parallel window processing. Default value --> 2.
Each thread owns an independent `VariantBuilder` instance with no shared mutable state. See [Performance & Parallelism](guides/architecture.md#7-performance-parallelism) for the threading architecture.

#### `-k`,`--min-kmer`
Minimum k-mer length to try for micro-assembly graph nodes. Default value --> 13. Allowed range: [13–253].
The graph construction starts at this k-mer size and increments by `--kmer-step` on retry. Smaller values increase sensitivity for short variants but produce more complex (slower) graphs.
See [K-mer Retry Cascade](guides/architecture.md#k-mer-retry-cascade) for the retry logic.

#### `-K`,`--max-kmer`
Maximum k-mer length to try for micro-assembly graph nodes. Default value --> 127. Allowed range: [15–255].
If no cycle-free graph is produced by this k-mer size, the window yields no assembled haplotypes. Higher values can resolve longer repeat motifs at the cost of requiring longer exact matches in reads.
See [K-mer Retry Cascade](guides/architecture.md#k-mer-retry-cascade).

#### `-s`,`--kmer-step`
K-mer step size for the retry cascade. Must be one of {2, 4, 6, 8, 10}. Default value --> 6.
Smaller steps try more intermediate k values before exhausting the search space — more chances to find a clean graph, but more rebuild iterations.
See [K-mer Retry Cascade](guides/architecture.md#k-mer-retry-cascade).

#### `--min-anchor-cov`
Minimum coverage for source/sink anchor nodes in the De Bruijn graph. Default value --> 5.
Anchors are reference k-mers at the start and end of the assembled component — they must be well-supported to define reliable walk boundaries. Lower values may enable assembly in low-coverage regions but risk anchoring on noise.

#### `--min-node-cov`
Minimum coverage for non-anchor nodes in the De Bruijn graph. Default value --> 2.
Nodes below this threshold are pruned during the [Graph Pruning Pipeline](guides/architecture.md#graph-pruning-pipeline). Higher values aggressively filter noise (faster runtime) but reduce sensitivity for subclonal variants with low allele frequency.

#### `--max-sample-cov`
Maximum per-sample coverage before downsampling. Default value --> 1000.
Windows exceeding this threshold per sample are downsampled using a deterministic paired strategy (fixed seed for reproducibility). Both mates of a pair are symmetrically accepted or rejected.
See [Read Filtering & Downsampling](guides/read_filtering.md) for the full downsampling algorithm.

### Flags

#### `--verbose`
Turn on verbose logging. Emits per-window status messages (skipped/assembled/genotyped) and graph complexity metrics to stderr.

#### `--extract-pairs`
Extract out-of-region mate reads for discordant and supplementary-aligned pairs.
When enabled, reads that are **not** in a proper pair or have a supplementary alignment (`SA` tag) trigger targeted `htslib` queries to retrieve their out-of-region mates. This captures split reads spanning structural variant breakpoints at window edges, at the cost of additional I/O per discordant pair.
See [Read Filtering & Downsampling](guides/read_filtering.md) for details.

#### `--no-active-region`
Force assembly of all windows regardless of mutation evidence.
By default, Lancet2 skips windows where no read shows variation (the [Active Region Detection](guides/active_region.md) heuristic). This flag disables that fast-skip, forcing assembly of every window. Useful for completeness auditing but causes **5–10× slower** runtime on WGS.
See [Active Region Detection](guides/active_region.md) for the heuristic algorithm and the MD tag requirement.

#### `--no-contig-check`
Skip contig name validation between the reference FASTA and BAM/CRAM headers.
Use when contig naming conventions differ across files (e.g., `chr1` vs `1`). Without this flag, mismatched contig names cause Lancet2 to exit with an error.

### Optional

#### `--graphs-dir`
Output directory to write per-window graphs in DOT and GFA format.
Must be a non-existing directory path that will be created. Produces two subdirectories: DOT files for Graphviz visualization (one per pruning stage) and GFA files for the SPOA POA graph.
See [Custom Visualization](guides/custom_visualization.md) for rendering instructions and interpretation.

#### `--enable-graph-complexity-features`
Emit `GRAPH_CX` INFO tag with per-variant graph complexity metrics (GEI, TipToPathCovRatio, MaxSingleDirDegree).
Topology-derived features that are coverage-stable above 20×.
See [Graph Complexity](guides/graph_complexity.md) and [VCF Output Reference](guides/vcf_output.md#optional-complexity-annotations) for details.

#### `--enable-sequence-complexity-features`
Emit `SEQ_CX` INFO tag with 11 multi-scale sequence complexity metrics (HRun, entropy, TR motifs, LongdustQ).
These features are perfectly coverage-invariant (sequence-only) and designed for ML-based variant filtering.
See [Sequence Complexity](guides/sequence_complexity.md) and [VCF Output Reference](guides/vcf_output.md#optional-complexity-annotations) for details.

#### `--genome-gc-bias`
> [0.0-1.0]. Default value --> 0.41

Global genome GC fraction for LongdustQ score correction.
Set to 0.5 to disable GC correction (uniform model).
Default value (0.41) is the human genome-wide GC average. Adjust for non-human genomes or targeted panels with significantly different GC content.

```bash
Lancet2 pipeline \
    --enable-graph-complexity-features \
    --enable-sequence-complexity-features \
    --genome-gc-bias 0.41 \
    --normal normal.bam --tumor tumor.bam \
    --reference ref.fasta --region "chr22" \
    --out-vcfgz output.vcf.gz
```

## VCF Output

See the [VCF Output Format](guides/vcf_output.md) guide for complete documentation
of all INFO and FORMAT fields.
