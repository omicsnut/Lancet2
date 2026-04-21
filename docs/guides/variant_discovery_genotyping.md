# Variant Discovery & Genotyping

This guide covers the two core algorithmic phases that transform assembled haplotypes into genotyped VCF records:

1. **Phase 1 — Variant Discovery**: Multiple Sequence Alignment (MSA) and multiallelic variant extraction from the Partial Order Alignment (POA) graph.
2. **Phase 2 — Genotyping**: Read-to-haplotype alignment via minimap2 and allele assignment via the combined scoring function.

---

## Phase 1: Multiple Sequence Alignment

After the de Bruijn graph walk enumerator produces a set of haplotype sequences (one reference + N assembled ALT paths per connected component), Lancet2 must determine **where** these haplotypes differ from the reference and extract those differences as VCF-ready variant records.

### SPOA Alignment Engine

Haplotype sequences are aligned into a Multiple Sequence Alignment using **SPOA** (SIMD Partial Order Alignment) with custom convex (dual-affine) gap scoring. The POA graph captures the alignment as a Directed Acyclic Graph (DAG) of nucleotide nodes, where shared sequences converge onto shared paths and divergent mutations open topological "bubbles."

### SPOA Scoring Parameters

The scoring parameters are **not** the standard SPOA defaults — they were specifically tuned for micro-assembly variant extraction:

| Parameter | Value | Standard Value | Rationale |
|:----------|:------|:---------------|:----------|
| Match | **0** | +2 | Keeps SPOA in the faster int16 SIMD path (see below) |
| Mismatch | −6 | −4 | Moderate penalty keeps MSA intact through dense MNVs |
| Gap Open 1 | −6 | — | Short-gap model: penalizes sequencer noise micro-indels |
| Gap Extend 1 | −2 | — | Short-gap model: moderate extension cost |
| Gap Open 2 | −26 | — | Long-gap model: high upfront cost, but... |
| Gap Extend 2 | −1 | — | ...cheap per-base extension for large biological indels |

#### Why Match = 0 (SIMD Lane Width)

SPOA 4.1.5 dynamically dispatches between **int16** and **int32** SIMD lanes based on `WorstCaseAlignmentScore()`. Setting Match=0 keeps all runtime scores non-positive, which keeps the worst-case score above the int16 threshold (−32768) — so the engine selects the faster int16 path (16 lanes per AVX2 register) instead of falling back to int32 (8 lanes, half throughput).

All other parameters are shifted downward by −2 to compensate: a standard +2/−4 match/mismatch scheme becomes 0/−6.

#### Dual-Affine Gap Model: The 20 bp Crossover

The two affine gap models target fundamentally different biological events:

- **Model 1** (`O₁=−6, E₁=−2`): Penalizes short (1–20 bp) gaps that typically arise from sequencer noise — homopolymer stutters, PCR slippage, and alignment artifacts. The moderate extension penalty applies friction that prevents false 1 bp gaps from opening in noisy regions.

- **Model 2** (`O₂=−26, E₂=−1`): Has a steep open cost but a cheap extension cost, allowing the aligner to extend through large biological indels (deletions >20 bp, large insertions) without truncating the alignment.

The two models **intersect at exactly 20 bp**: solving `6 + 2L = 26 + 1L` yields `L = 20`. For gaps shorter than 20 bp, Model 1 dominates (stricter). For gaps longer than 20 bp, Model 2 dominates (cheaper extension). The aligner takes the minimum cost at each gap length, automatically selecting the appropriate model.

This design enables detection of somatic deletions up to ~500 bp within the ~1 kbp micro-assembly window context.

**Complexity:** `O(N × L²)` where N = number of haplotypes and L = contig length.

---

## Multiallelic Variant Extraction

After the MSA, variants are extracted from the POA graph's topology using a **Greedy Sink** algorithm that sweeps pointer arrays across all haplotype paths simultaneously.

### The Sweep Algorithm

The extractor maintains an array of active pointers — one per haplotype — initialized to each path's first node in the topologically-sorted POA graph:

1. **Convergence check**: If all pointers reference the same graph node, the paths agree at this position. The node is recorded as the last "match" anchor and all pointers advance.

2. **Bubble detection**: If pointers disagree (they reference different nodes), a topological bubble has opened — this is a variant.

3. **Greedy sink**: The pointer with the **lowest topological rank** (leftmost in 5′→3′ order) advances first, appending its character to a per-path string buffer. This process repeats until all pointers converge onto the same node again — the bubble's "sink."

4. **Normalization**: The accumulated per-path strings form a raw REF/ALT alignment. VCF parsimony (right-trim, then left-trim) is applied **globally** across all alleles simultaneously. This ensures consistent minimal representation for multiallelic records.

5. **Classification**: Each ALT allele's "Sequence Core" (the mutation stripped of shared padding) is classified as SNV, INS, DEL, MNP, or CPX via `ClassifyVariant`.

### The Multiallelic Shielding Problem

When multiple ALT alleles share a locus, VCF normalization must operate across the **entire multiallelic cluster simultaneously**. Consider:

```
REF: ATGC    ALT1: ACGC    ALT2: ATCC
```

If normalized independently, ALT1 would trim to `T→C` at position 2 and ALT2 would trim to `G→C` at position 3. But in a VCF multiallelic record, they must share the same POS and REF. The `NormalizeVcfParsimony` method handles this by only trimming characters that **all** alleles share with the reference — it stops the moment any single allele's trimmed sequence would become empty or lose its VCF anchor base.

The Sequence Core is computed per-allele using an **O(N)** 5′/3′ squeeze to determine the true biological variant length independent of the VCF padding context.

---

## Phase 2: Read-to-Haplotype Genotyping

Once variants are discovered, Lancet2 must determine which allele each read supports. This is done by realigning every read to every assembled haplotype using **minimap2**, then computing a combined scoring function to assign each read to its best-matching allele.

### The Alignment Paradigm Shift

There is an intentional contrast between the MSA (Phase 1) and read-to-haplotype (Phase 2) alignment parameters:

| Property | Phase 1 (MSA) | Phase 2 (Genotyping) |
|:---------|:--------------|:---------------------|
| **Source of divergence** | True biological mutations | Sequencer noise and artifacts |
| **Strategy** | Forgiving (cheap gaps, generous mismatches) | Strict (heavy penalties for noise) |
| **Goal** | Force end-to-end alignment to discover all variants | Precisely assign reads to their parent haplotype |

In Phase 1, the assembled contigs are high-confidence and divergence is real biology — so we use cheap gap extensions to capture large indels. In Phase 2, the biological variants are **already baked into the haplotype sequences** — divergence between a read and its parent haplotype is sequencer error, adapter contamination, or chimeric artifacts. Heavy gap/mismatch penalties prevent noisy reads from "bleeding" to the wrong allele.

### minimap2 Parameter Overrides

| Parameter | Override | Default | Rationale |
|:----------|:---------|:--------|:----------|
| Flag | `MM_F_SR` | 0 | Activates the SR extension-region code path: full-query boundaries (`qs0=0, qe0=qlen`) and `end_bonus`-based reference expansion at read edges. Disables irrelevant long-read paths (inversion detection, `bw_long` re-chain). All 10 SR code paths are verified safe for single-segment haplotype alignment. |
| Z-Drop | 100,000 | 400 | Effectively disables DP truncation. A 300 bp somatic deletion incurs gap penalty O + 300·E = 912, exceeding the default zdrop=400. |
| Bandwidth (`bw`) | 10,000 | 500 | Envelopes insertions up to ~2 kbp. Prevents the banding boundary from terminating alignment within large assembled insertions. |
| Seed k/w | 11/5 | 15/10 | Increases sensitivity for highly mutated fragments. Maps reads through dense mutation clusters and STRs where 15 bp exact matches are rare. Per-haplotype indexing makes k=11's higher false-positive rate harmless. |
| Gap model | Single-affine (12/3) | Dual-affine | In Phase 2, gaps are noise (not biology). A single strict affine model penalizes all gaps uniformly without a cheap-extension path. |
| `end_bonus` | 10,000 | -1 | Forces KSW2's EXTZ_ONLY mode to always backtrack to the query end instead of the max-score cell. 10,000 is 66× the max theoretical alignment score (151), keeping `mqe + end_bonus` within int32 range. `INT_MAX` causes signed overflow (undefined behavior) at two call sites: `ksw2_extd2_sse.c:393` and `align.c:699-702`. |
| `max_gap` | 200 | 5,000 | Caps the backward search radius in `mg_lchain_dp` on the **query** dimension. No valid chain on a 151 bp read spans a 200 bp query gap. |
| `max_gap_ref` | 5,000 | −1 (defaults to `max_gap`) | Caps the chaining gap on the **reference (haplotype)** dimension. Must be set independently because minimap2 defaults it to `max_gap` when ≤0 (`map.c:271`). With `max_gap_ref=200`, seeds separated by >200 bp on the haplotype cannot chain, blocking alignments across insertions >200 bp. Since assembled haplotypes range 200–2000 bp and routinely contain large InDels, 5,000 (minimap2's own default) covers all cases. |
| `best_n` | 1 | — | Keeps only the single best hit per haplotype. |

Since alignment is restricted to the local contig window, these inflated parameters have minimal runtime impact compared to whole-genome alignment.

### Exhaustive Haplotype Alignment

Each read is aligned to **every** haplotype — REF and all ALTs — independently. There is no early-exit: even if a read perfectly matches one haplotype, it must still be scored against all others to:

- Compute the **reference edit distance** (NM against REF haplotype), used for the ASMD FORMAT field.
- Guarantee correct relative scoring so a read isn't incorrectly assigned due to incomplete cross-haplotype comparison.

### The Combined Scoring Function

For each read × variant × haplotype combination, the genotyper computes a combined score:

```
combined = (global_score − local_raw_score − sc_penalty) + (local_pbq_score × local_identity)
```

Each component captures a distinct signal:

| Component | Definition | Purpose |
|:----------|:-----------|:--------|
| `global_score` | minimap2 DP score of the full read→haplotype alignment | How well the entire read fits this haplotype, including flanking context |
| `local_raw_score` | Raw substitution-matrix score (M/=/X ops only) within the variant's physical boundaries on the haplotype. Gap penalties (I/D) are excluded. | **Subtracted** to remove the variant region's contribution from the global score so the PBQ-weighted replacement does not double-count it. Excluding gap penalties prevents the subtraction from refunding them: for a 150bp insertion, `global − (−450) = global + 450` would inflate the score by 450, making it indistinguishable from a perfect match. |
| `sc_penalty` | Soft-clip penalty (clipped bases × mismatch cost) | **Subtracted** to crush chimeric supplementary mappings that only partially align |
| `local_pbq_score` | PBQ-weighted DP score within the variant region: each base's substitution matrix contribution is scaled by Phred confidence `(1 − 10^(−Q/10))` | High-confidence bases contribute more; Q5 bases are down-weighted to near zero |
| `local_identity` | Fraction of exact matches within the variant region from the CIGAR trace | **Gates** the PBQ score — a high PBQ score from a fragmented, low-identity alignment (e.g., STR artifact) is discounted toward zero |

#### Why Subtract `local_raw_score`?

The `global_score` from minimap2 already includes the raw alignment cost spanning the variant region. If `local_pbq_score` were added on top, the variant region would be **double-counted**. Subtracting `local_raw_score` removes the variant region's contribution from the global score, allowing the higher-fidelity PBQ-weighted local score to replace it.

`local_raw_score` captures only substitution-matrix scores from M/=/X CIGAR operations — gap penalties (I/D) are excluded. If gap penalties were included, the subtraction would refund them: `global − (−gap_penalty) = global + gap_penalty`. For a 150bp insertion, this is a +450 refund that makes a read with a large gap score nearly as high as a perfect match, destroying allele discrimination for large InDels.

#### Why Subtract `sc_penalty`?

Soft-clipped read tails are unaligned sequence. Minimap2 exempts them from its DP score. Without an explicit penalty, a chimeric read that only partially aligns to a haplotype could achieve a competitive global score against a clean end-to-end alignment. The soft-clip penalty crushes these partial matches.

#### Why `local_pbq_score × local_identity`?

This product acts as a **confidence gate**. In low-complexity repeat regions (STRs, VNTRs), a noisy read might accumulate a high raw DP score simply by traversing a chaotic repeat locus. The identity fraction measures alignment "cleanliness" — the fraction of bases that are exact matches. A perfect alignment (identity = 1.0) retains the full PBQ weight; a fragmented, heavily-gapped alignment (identity < 0.7) is aggressively discounted, preventing repeat-region artifacts from winning allele assignments.

### Base Quality at Variant Position

The genotyper extracts a representative base quality for each read at each variant:

- **SNVs**: The single base quality at the variant position.
- **Indels**: The **minimum** base quality across all read positions spanning the variant region (weakest-link summary).

The minimum is used because the confidence in observing a complete indel is bounded by the least confident base in the region. This follows the same convention as `bcftools mpileup` and ensures that downstream PL and PBQ computations correctly treat each read as one independent observation.

### Folded Read Position

For the RPCD (Read Position Cohen's D) FORMAT field, the variant's query position on the read is converted to a **folded position**: `min(p, 1−p)` where `p = variant_query_pos / read_length`.

- 0.0 = variant at read edge
- 0.5 = variant at read center

Folding maps both read ends to the same low-value space. Without folding, artifacts clustering at positions 5 and 145 in a 150 bp read would average to ~75 — indistinguishable from a centered variant. With folding, both map to ~0.03, creating a clear signal that ALT alleles are systematically positioned near read edges — a hallmark of alignment artifacts.

### Genotype Likelihood Model

After allele assignment, genotype likelihoods are computed using a **Dirichlet-Multinomial (DM)** count-based model that replaces the per-read product structure common to both pileup-based and haplotype-aware variant callers.

#### Why Not the Per-Read Product Model?

Variant callers — whether pileup-based or haplotype-aware — compute the total genotype likelihood as a **product of independent per-read terms**:

```
P(Data | GT) = Π P(read_i | GT)
```

Under the simplest error model, the per-read likelihood for a diploid genotype `GT = (a₁, a₂)` expands to:

```
P(read | GT=(a₁,a₂)) = 0.5 × P(read | a₁) + 0.5 × P(read | a₂)
P(read | allele a) = (1−ε) if read matches a, else ε/(K−1)
```

Haplotype-aware callers replace this inner `P(read | allele)` with more sophisticated likelihoods (e.g., alignment-based HMM scores that integrate over gaps and mismatches), but the **outer product structure** is shared. Converting to Phred-scaled PLs via `PL = -10 × log₁₀(P)`, each read adds an independent contribution to the sum, so **PLs scale linearly with depth**. A true heterozygous variant at 30× might produce PL[0/0]=90; at 3000× (the same variant, just deeper), PL[0/0]=9,000. This creates two problems regardless of how sophisticated the per-read likelihood is:

1. **Runaway scaling.** Ultra-deep targeted panels and UMI-deduplicated libraries produce PL values in the tens of thousands. Downstream ML models trained on 30× WGS data see entirely different feature distributions at 2000×.
2. **Missing overdispersion.** The per-read product model assumes reads are independent observations. In reality, correlated errors (flow-cell lane effects, PCR duplication, systematic mismapping) cause allele count variance to exceed the independent-observation expectation at high depth. No per-read product model can capture this — it treats every read as equally informative.

#### The Dirichlet-Multinomial (DM) Model

The DM model replaces per-read iteration with a single count-based evaluation. For K alleles with observed counts `c = (c₀, c₁, ..., c_{K−1})`, expected allele fractions `μ_i` under a genotype hypothesis, and precision parameter `M = (1−ρ)/ρ` controlling overdispersion:

```
ln P(c | μ, M) = lnΓ(ΣMμ_i) − lnΓ(N + ΣMμ_i) + Σ[lnΓ(c_i + Mμ_i) − lnΓ(Mμ_i)]
```

For each diploid genotype `(a₁, a₂)`:

- **Homozygous** (`a₁ == a₂`): `μ[a₁] = 1 − ε`, `μ[other] = ε/(K−1)`
- **Heterozygous**: `μ[a₁] = μ[a₂] = (1−ε)/2`, `μ[other] = ε/(K−1)`

**In plain terms:** imagine drawing N marbles (reads) from a bag of colored marbles (alleles). Each genotype hypothesis predicts a mix of colors, and the DM formula asks: "how surprised am I by what I actually drew?" The `lnΓ` (log-gamma) terms account for all the ways the observed counts could arise from the expected mix. Unlike the per-read product model which treats each draw independently, the DM model says "draws from the same sequencing run are correlated" — like asking 100 people vs. 10,000 the same question: the extra 9,900 answers share the same biases (same methodology, same population) and don't make you 100× more confident. In sequencing, those shared biases are PCR duplicates, flow-cell artifacts, and systematic mismapping.

This produces two key improvements over the per-read product model:

1. **Natural asymptote.** As depth N → ∞, the lgamma ratio terms converge, causing PLs to plateau at a ceiling determined by ρ. No artificial `PL / SAMPLE_DP` normalization is needed.
2. **Correlated error absorption.** The overdispersion parameter ρ (default 0.01, precision M=99) widens the count variance beyond the independent-observation expectation. At 2000× depth, additional reads of the same evidence produce diminishing confidence returns — matching the empirical reality of sequencing data.

The default parameters `ε = 0.005` (background error) and `ρ = 0.01` (overdispersion) are tuned for Illumina short-read sequencing. See the in-code documentation in `variant_support.cpp` for guidance on adjusting these for ONT, HiFi, or ultra-deep targeted panels.

#### Continuous Mixture Log-Odds (CMLOD)

DM-based PLs solve the depth-scaling and overdispersion problems, but they evaluate only rigid diploid genotype states (0%, 50%, 100% ALT). For somatic and mosaic variants at continuous frequencies like 2% or 15% VAF, none of these states represents the true frequency. Each ALT allele therefore also receives a **CMLOD** (Continuous Mixture Log-Odds) score that operates over continuous allele frequency space and integrates exact per-read base qualities:

```
P(read called as s | f) = Σ_t f[t] × P(s | true origin = t)
CMLOD[alt] = max(0, LL(f_MLE) − LL(f_null))    (log₁₀ scale)
```

**In plain terms:** PLs answer "which genotype best explains the data?" — choosing between fixed states like 0% ALT, 50% ALT, or 100% ALT. CMLOD answers a different question: "is there *any* evidence this ALT allele is real, regardless of its frequency?"

It works by comparing two scenarios:

- **Scenario A (MLE):** the variant exists at whatever frequency the data suggests (e.g., 15% VAF)
- **Scenario B (null):** the variant doesn't exist at all — those reads are just errors

CMLOD is the log-odds between these two scenarios. Crucially, each read's contribution is weighted by its base quality: a high-confidence Q40 read counts for much more than a noisy Q10 read. This means CMLOD naturally separates true low-frequency variants from sequencing noise — a mosaic variant at 5% VAF supported by 10 high-quality reads will score much higher than 10 noisy reads at the same frequency.

A Q40 read supporting ALT contributes ~4.0 to the LOD; a Q10 read contributes ~0.5. This separates high-confidence signal from sequencing noise without requiring an external error model.

CMLOD is emitted as a `Number=A` FORMAT field (one value per ALT allele). It is particularly useful for low-VAF mosaic variant detection and somatic variant calling, where the discrete diploid PL model cannot represent continuous frequency states like 2% or 15% VAF.

* **Read more:** [Alignment-Derived Annotations](alignment_annotations.md), [VCF Output Reference](vcf_output.md)
