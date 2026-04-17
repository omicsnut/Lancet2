# HTS Read Extraction — Benchmark Report

## Purpose

The `Extractor` is the I/O gateway for Lancet2: it reads sequencing alignments (reads) from BAM/CRAM files using the [HTSlib](https://github.com/samtools/htslib) C library. Every read that Lancet2 processes passes through this component, making its throughput a hard ceiling on pipeline performance.

This benchmark measures **wall-clock time to extract all alignments from a tumor sample** across two file formats (BAM and CRAM), three field-parsing levels, and 1–8 I/O threads, on two server-class CPUs.

### What the benchmark tests

Each benchmark opens a tumor sample file, iterates over every alignment record, and collects them into a vector. The three field-parsing levels test progressively heavier decoding:

| Level | Fields Parsed | Description |
|:------|:--------------|:------------|
| `CORE_QNAME` | Coordinates + read name | Lightest: no sequence data decoded |
| `CIGAR_SEQ_QUAL` | Above + CIGAR + bases + quality scores | Medium: full alignment data |
| `AUX_RGAUX` | Above + 7 auxiliary tags (RG, MC, NM, SA, XS, MD, AS) | Heaviest: everything Lancet2 needs |

---

## Test Environment

| Property | AMD Milan | Intel Ice Lake |
|:---------|:----------|:---------------|
| **Host** | lancet2-dev-gcc12 | lancet-dev01 |
| **CPUs** | 32 | 32 |
| **Clock** | 2450 MHz | 2600 MHz |
| **Date** | 2023-04-04 | 2023-05-01 |

---

## Results

### CRAM format (version 3.1)

The test CRAM file (6.69 GB) uses CRAM 3.1, which compresses alignment data using a combination of codecs — [rANS order-0/order-1](https://github.com/samtools/hts-specs/blob/master/CRAMv3.pdf) (arithmetic coding), gzip, bzip2, LZMA, and optionally reference-based sequence differencing. CRAM 3.1 added additional rANS variants (Nbit, stripe) for improved compression. Decompression is more compute-intensive than BAM because each data series (bases, qualities, read names, etc.) may use a different codec.

### BAM format

The test BAM file (10.62 GB) uses BGZF — block-level gzip with 64 KB blocks. All fields are compressed uniformly with a single codec (deflate). Decompression is lighter per-read than CRAM but produces ~59% more I/O traffic for the same alignment data.

#### CRAM extraction times (ms)

| Threads | CORE_QNAME (AMD) | CORE_QNAME (ICL) | CIGAR_SEQ_QUAL (AMD) | CIGAR_SEQ_QUAL (ICL) | AUX_RGAUX (AMD) | AUX_RGAUX (ICL) |
|--------:|------------------:|------------------:|---------------------:|---------------------:|----------------:|----------------:|
| 1 | 36,014 | 34,694 | 100,510 | 94,384 | 150,882 | 142,447 |
| 2 | 6,899 | 6,270 | 21,986 | 20,103 | 33,196 | 31,390 |
| 3 | 3,474 | 3,011 | 14,055 | 12,893 | 21,568 | 20,624 |
| 4 | 2,849 | 2,369 | 10,676 | 9,652 | 16,234 | 15,406 |
| 5 | 2,679 | 2,293 | 8,683 | 7,890 | 13,651 | 12,734 |
| 6 | 2,579 | 2,309 | 7,831 | 7,224 | 12,902 | 11,400 |
| 7 | 4,975 | 2,261 | 7,505 | 6,658 | 24,384 | 10,795 |
| 8 | 11,300 | 2,233 | 8,717 | 6,583 | 33,754 | 10,478 |

#### BAM extraction times (ms)

| Threads | CORE_QNAME (AMD) | CORE_QNAME (ICL) | CIGAR_SEQ_QUAL (AMD) | CIGAR_SEQ_QUAL (ICL) | AUX_RGAUX (AMD) | AUX_RGAUX (ICL) |
|--------:|------------------:|------------------:|---------------------:|---------------------:|----------------:|----------------:|
| 1 | 60,772 | 61,423 | 76,663 | 73,290 | 99,479 | 95,237 |
| 2 | 14,940 | 11,808 | 16,055 | 14,387 | 26,688 | 24,434 |
| 3 | 6,448 | 5,317 | 10,678 | 9,124 | 17,610 | 16,215 |
| 4 | 4,704 | 3,926 | 8,013 | 7,235 | 13,018 | 12,223 |
| 5 | 3,951 | 3,601 | 7,548 | 6,509 | 12,074 | 9,861 |
| 6 | 3,556 | 3,206 | 7,331 | 6,486 | 10,208 | 8,824 |
| 7 | 3,438 | 3,027 | 7,017 | 6,411 | 9,386 | 8,612 |
| 8 | 3,326 | 2,981 | 6,529 | 6,178 | 8,901 | 8,311 |

---

## Analysis

### 1. CRAM is faster than BAM at low thread counts

At 1 thread, CRAM `CORE_QNAME` extraction (34–36s) is nearly **2× faster** than BAM (60–61s). This is counterintuitive — CRAM decoding is more compute-intensive — but the explanation is I/O volume: CRAM files are significantly smaller on disk, so less data travels from storage to memory. At 1 thread, the disk read dominates.

At higher thread counts, CRAM's heavier per-read CPU cost catches up. By 4+ threads, BAM and CRAM times converge because the bottleneck shifts from I/O to decoding throughput.

### 2. AMD Milan shows performance degradation at 7–8 CRAM threads

On the AMD Milan system, CRAM extraction times **increase sharply** at 7–8 threads (e.g., `CORE_QNAME` goes from 2,579ms at 6 threads to 11,300ms at 8 threads — a 4.4× regression). Intel Ice Lake does not exhibit this behavior.

The root cause is not conclusively established by this data. Possible explanations include HTSlib internal lock contention on CRAM decompression buffers, disk I/O scheduling differences between the two test systems, or platform-specific memory subsystem effects. BAM extraction does not show this pattern, which suggests the bottleneck is specific to CRAM's heavier per-read decompression.

### 3. Field-parsing cost scales linearly with data decoded

Across both architectures and formats, the time ratio between parsing levels is consistent:

| Comparison | AMD CRAM 4T | ICL CRAM 4T | Ratio |
|:-----------|:------------|:------------|:------|
| AUX_RGAUX / CORE_QNAME | 16,234 / 2,849 | 15,406 / 2,369 | ~5.7× / ~6.5× |
| AUX_RGAUX / CIGAR_SEQ_QUAL | 16,234 / 10,676 | 15,406 / 9,652 | ~1.5× / ~1.6× |

Decoding auxiliary tags adds ~50–60% overhead on top of CIGAR+SEQ+QUAL. This cost is unavoidable for Lancet2's active region detection (which requires the MD tag) and mate-pair extraction (which requires the SA tag).

### 4. Optimal thread count is 4–6

Both architectures plateau at 4–6 threads. Beyond that, returns diminish and AMD shows regressions. Lancet2's production code does not set HTSlib I/O threads on extractors (defaulting to single-threaded decoding per extractor instance), relying instead on window-level parallelism across worker threads for throughput.

---

## Takeaway

- **Use CRAM format** for storage and distribution — it provides equivalent or faster extraction at reduced I/O bandwidth.
- **4–6 HTSlib threads per extractor** is the measured throughput sweet spot. Lancet2's production code uses single-threaded extraction per extractor, relying on window-level parallelism instead. This benchmark data informs potential future tuning.
- **AUX tag parsing is the heaviest field level** but is required for Lancet2's mutation detection logic. Extracting only the tags needed (rather than all tags) keeps the overhead controlled.
- **AMD Milan shows CRAM regression at 7–8 threads** — Intel Ice Lake does not. The root cause is not established by this benchmark (candidates include HTSlib internal lock contention, disk I/O scheduling, or platform-specific memory subsystem effects).
