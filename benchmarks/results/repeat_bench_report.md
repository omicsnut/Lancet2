# Repeat Detection — Benchmark Report

## Purpose

Lancet2 is a variant caller that performs localized colored de Bruijn graph micro-assembly. It slices the genome into ~1 kbp windows, re-assembles all reads into a colored de Bruijn graph, and extracts variant haplotypes from the graph topology. Graph construction iterates k-mer sizes from 13 to 127 in steps of 6, retrying at larger k when the current k produces a tangled graph. Before each assembly attempt, Lancet2 checks whether the reference sequence contains repeated k-mers at the current k — repeats create cycles by construction, making assembly at that k pointless. Two functions implement this check:

- **`HammingDist`** — counts byte-level mismatches between two equal-length strings. Called inside the inner loop of repeat detection, making it the single hottest function in the pipeline (10.1% of total CPU in profiling).
- **`HasRepeat`** — detects whether any two k-mers in a sliding window are identical (exact repeat) or nearly identical (approximate repeat, within a mismatch threshold).

This report documents the empirical comparison of four `HammingDist` implementations and the `HasRepeat` behavior across entropy levels and mismatch thresholds. The results justify the production implementation choices in [`repeat.cpp`](../../src/lancet/base/repeat.cpp).

---

## Test Environment

| Property | Value |
|:---------|:------|
| **Host** | ni1dm5-001.nygenome.org |
| **CPU** | Intel Skylake, 144 cores, 2300 MHz |
| **Compiler** | GCC 15.2 with `-O3 -march=x86-64-v3` (AVX2) |
| **Date** | 2026-04-17 |
| **Seeds** | Deterministic (42, 137, 271) for reproducibility |

---

## Part 1: Hamming Distance — Accumulator Strategy Comparison

### The problem

A comparison like `first[i] != second[i]` produces a 1-byte boolean result (0 or 1). When accumulated into a wider type, the compiler must insert instructions to widen each byte to the accumulator width before adding. This "widening cascade" bloats the SIMD instruction count and wastes vector lane capacity.

### Approaches tested

| Variant | Accumulator | Strategy |
|:--------|:------------|:---------|
| **Production (u8 batch)** | `u8`, batched in groups of 255 | Proves to the compiler's value-range analysis that no byte lane overflows. Eliminates all widening. |
| **Usize accumulator** | `usize` (8 bytes) | Forces `vpmovzxbq` — expands each byte to 8 bytes. 46 widening instructions per loop iteration. |
| **U32 accumulator** | `uint32_t` (4 bytes) | Forces `vpmovzxbd` — expands to 4 bytes. Halves the widening cost but still leaves 18 instructions. |
| **SWAR + popcount** | `usize` with manual bit tricks | Processes 8 bytes per scalar iteration using XOR + nibble-fold + `std::popcount`. Defeats auto-vectorization entirely. Also has undefined behavior (strict aliasing, alignment). |

### Results — k-mer sized inputs (11–121 bp), medium entropy

These lengths match Lancet2's operating range ([`graph_params.h`](../../src/lancet/cbdg/graph_params.h): min=13, max=127).

| k (bp) | Production (ns) | Usize (ns) | U32 (ns) | SWAR (ns) |
|-------:|---------:|----------:|---------:|---------:|
| 11 | 4.3 | 4.3 | 4.9 | 5.3 |
| 21 | 5.2 | 4.4 | 5.2 | 7.9 |
| 31 | 5.3 | 6.1 | 5.4 | 9.2 |
| 41 | 5.7 | 9.1 | 6.2 | 9.2 |
| 51 | 10.1 | 12.6 | 8.4 | 12.9 |
| 61 | 10.3 | 12.6 | 8.7 | 12.3 |
| 71 | 10.5 | 12.6 | 10.1 | 14.7 |
| 81 | 10.6 | 17.3 | 10.1 | 16.2 |
| 91 | 10.6 | 17.6 | 10.2 | 18.2 |
| 101 | 13.2 | 23.2 | 12.3 | 20.4 |
| 111 | 12.6 | 22.1 | 12.3 | 22.1 |
| 121 | 12.3 | 22.1 | 12.4 | 24.8 |

At short lengths (k < 31), all four approaches perform similarly because the data fits in a single SIMD register and the widening overhead is negligible. The gap widens as k increases. At k=101 (a typical Lancet2 operating point):

| Variant | Time (ns) | Relative to Production |
|:--------|----------:|-----------------------:|
| Production (u8 batch) | 13.2 | 1.0× |
| U32 accumulator | 12.3 | 0.93× |
| Usize accumulator | 23.2 | **1.76×** slower |
| SWAR + popcount | 20.4 | **1.55×** slower |

The U32 approach is within noise of production at these lengths. The real separation appears at longer inputs.

### Results — power-of-two stress test (8–2048 bytes), medium entropy

These lengths exceed Lancet2's k-mer range but stress-test the batching logic and reveal scaling behavior.

| Bytes | Production (ns) | Usize (ns) | U32 (ns) | SWAR (ns) |
|------:|---------:|---------:|---------:|---------:|
| 8 | 4.3 | 4.1 | 4.8 | 4.9 |
| 16 | 4.3 | 4.4 | 4.8 | 6.4 |
| 32 | 5.5 | 5.7 | 5.2 | 8.9 |
| 64 | 8.3 | 12.7 | 8.3 | 13.4 |
| 128 | 10.9 | 21.2 | 12.2 | 23.6 |
| 256 | 17.7 | 41.7 | 21.1 | 54.4 |
| 512 | 30.7 | 94.4 | 40.8 | 88.2 |
| 1024 | 54.5 | 169.2 | 74.1 | 180.1 |
| 2048 | 105.5 | 349.1 | 145.4 | 357.1 |

At 2048 bytes, the separation is clear:

| Variant | Time (ns) | Relative to Production |
|:--------|----------:|-----------------------:|
| Production (u8 batch) | 105.5 | 1.0× |
| U32 accumulator | 145.4 | 1.38× slower |
| Usize accumulator | 349.1 | **3.31×** slower |
| SWAR + popcount | 357.1 | **3.38×** slower |

### Entropy independence

HammingDist processes every byte unconditionally — it does not short-circuit on match/mismatch. The data confirms this: across low, medium, and high entropy sequences, all four variants produce identical timings (within measurement noise). The entropy generators affect only `HasRepeat` behavior, not `HammingDist`.

### Why the production implementation wins

The u8 batch accumulation eliminates the widening cascade by proving to the compiler's value-range analysis that each byte lane stays in [0, 255]. On GCC 15.2 with AVX2, this reduces the SIMD instruction count from 168 (usize) to 36 (u8 batch) per vectorized iteration — a 4.7× reduction in instruction count that translates directly to the observed 3.3× wall-clock improvement at scale.

The SWAR approach defeats auto-vectorization entirely by using `reinterpret_cast<u64*>` loads and scalar `std::popcount`, forcing the compiler into a serial dependency chain. It also contains three instances of undefined behavior (buffer over-read, alignment violation, strict aliasing violation).

---

## Part 2: HasRepeat — Repeat Detection Across Entropy and Mismatch Thresholds

### The problem

Before building the de Bruijn graph at a given k-value, check whether any two k-mers in the reference window are identical or nearly identical. If so, the graph would contain cycles, making assembly at this k pointless. The function must handle both exact repeats (distance = 0) and approximate repeats (distance ≤ `max_mismatches`).

### Implementation strategy

`HasRepeat` uses a **dual strategy** based on the mismatch threshold. For exact repeats (`max_mismatches = 0`), it uses an O(n) hash-set duplicate check — each k-mer is inserted into an `absl::flat_hash_set`, and the function returns `true` on the first insertion failure (duplicate). For approximate repeats (`max_mismatches > 0`), it falls back to O(n(n-1)/2) pairwise Hamming distance comparisons using upper-triangle indexing, short-circuiting on the first qualifying pair.

### Benchmark parameters

- **Sequence length**: 600 bp (typical Lancet2 genomic window)
- **K-mer sizes**: 13, 55, 97, 127 (from [`graph_params.h`](../../src/lancet/cbdg/graph_params.h))
- **Mismatch thresholds**: 0 through 8
- **Entropy levels**:
  - **Low**: 75% A, ~8% each C/G/T — creates homopolymer runs and frequent repeats
  - **Medium**: uniform random DNA — realistic genomic sequence
  - **High**: uniformly sampled from a 16-symbol dinucleotide alphabet — few repeats

### Results — k=13 (588 k-mers per window)

At k=13, there are only 4^13 ≈ 67 million possible k-mers. With 588 drawn from a 600bp window, the probability of collisions depends heavily on entropy.

| Mismatches | Low (ns) | Medium (ns) | High (ns) |
|-----------:|---------:|------------:|----------:|
| 0 | 390 | 5,200 | 5,505 |
| 1 | 513 | 1,121,772 | 1,132,201 |
| 2 | 95 | 740,170 | 346,952 |
| 3 | 17 | 21,868 | 39,504 |
| 4 | 17 | 14,283 | 10,715 |
| 5 | 11 | 5,071 | 75 |
| 6 | 11 | 11 | 76 |
| 7 | 11 | 11 | 75 |
| 8 | 11 | 11 | 75 |

**Key observations:**

1. **`mm=0` (exact repeat)** uses the O(n) hash-set fast path. Low entropy finds a duplicate early (390 ns). Medium/high entropy require scanning more k-mers before finding a duplicate (5.2–5.5 µs) but are still **200× faster** than the O(n²) pairwise approach that would take ~1.1ms.

2. **`mm=1` is the worst case for medium/high entropy** — almost no two random 13-mers differ by exactly 0 or 1 positions, so the function scans nearly all 172K pairs before returning `false` (~1.1ms).

3. **Monotonic decrease with increasing mismatch threshold** — as the threshold rises, more pairs qualify, and the short-circuit triggers earlier. At `mm=6+` with medium entropy, the very first pair checked already qualifies (11 ns ≈ 1 HammingDist call).

4. **Low entropy short-circuits immediately at all thresholds** — the biased base distribution creates many near-identical k-mers, so even `mm=1` finds a match quickly. The 513 ns for `mm=1` (vs 390 ns for `mm=0`) reflects the cost difference between pairwise scanning (finds match in ~2 pairs) vs hash-set lookup.

### Results — k=55, 97, 127 (fewer k-mers per window)

At larger k, the sliding window produces fewer k-mers (600 - k + 1), and the k-mer space grows exponentially (4^k), making collisions extremely unlikely.

| k | K-mers | mm=0 (ns, all entropies) | mm=1–8 (ns, all entropies) |
|--:|-------:|-------------------------:|---------------------------:|
| 55 | 546 | 6,900–7,200 | ~1,430,000 (flat) |
| 97 | 504 | 8,700–8,800 | ~804,000 (flat) |
| 127 | 474 | 8,200–8,300 | ~1,326,000 (flat) |

At k ≥ 55, approximate repeats (`mm ≥ 1`) never short-circuit — the k-mer space is too vast for any two k-mers to differ in ≤ 8 positions by chance. The function exhausts all pairs every time, producing a flat timing profile regardless of mismatch threshold or entropy.

The `mm=0` hash-set path remains fast (7–9 µs) because it performs O(n) insertions regardless of k-mer size. The pairwise path at k=97 (504 k-mers, 126,756 pairs) costs ~804 µs, while k=127 (474 k-mers, 112,101 pairs) costs ~1,326 µs — the per-pair cost scales with k because each HammingDist call processes k bytes.

---

## Production Configuration

In [`graph.h`](../../src/lancet/cbdg/graph.h), the production call site uses `NUM_ALLOWED_MISMATCHES = 3`:

```cpp
static constexpr usize NUM_ALLOWED_MISMATCHES = 3;
return lancet::base::HasRepeat(absl::MakeConstSpan(klen_seqs), NUM_ALLOWED_MISMATCHES);
```

At `mm=3`, the benchmark shows:
- **Low entropy, k=13**: 17 ns (short-circuits immediately)
- **Medium entropy, k=13**: 21,868 ns (scans ~1,800 pairs before finding a match)
- **High entropy, k=13**: 39,504 ns (scans ~3,200 pairs)
- **k ≥ 55, all entropies**: 800–1,467 µs (exhausts all pairs, no match found)

For the typical Lancet2 pipeline, `HasRepeat` with `mm=3` serves as a fast filter: most real genomic windows at small k contain near-duplicate k-mers (especially in homopolymer-rich regions), so the short-circuit triggers early. At larger k-values where repeats are rare, the function runs to completion but the cost is acceptable because larger k-values are tried fewer times (the outer k-stepping loop uses increments of 6).

In [`variant_builder.cpp`](../../src/lancet/core/variant_builder.cpp), a separate call uses `HasExactRepeat` (equivalent to `HasRepeat(kmers, 0)`) with the O(n) hash-set path, which completes in 5–9 µs regardless of sequence content.

---

## Design Decisions

1. **u8 batch accumulation over usize/u32/SWAR for HammingDist** — the only approach that eliminates all widening instructions. 3.3× faster than usize at 2048 bytes. No intrinsics, no undefined behavior, portable across x86 AVX2 and ARM64 NEON.

2. **Dual-strategy HasRepeat** — O(n) hash-set for exact repeats, O(n²) pairwise scan for approximate repeats. The hash-set path is 200× faster than pairwise for the exact case. The pairwise path with short-circuit is optimal for approximate repeats where the first qualifying pair is found early.

3. **Upper-triangle indexing** — comparing only pairs (i, j) where i < j halves the comparison count from n² to n(n-1)/2. Combined with short-circuit return, the typical case at small k (where repeats exist) terminates early at a fraction of the full cost.

4. **Deterministic benchmark seeds** — fixed PRNG seeds (42, 137, 271) ensure reproducible results across runs, enabling meaningful A/B comparison when testing implementation changes.
