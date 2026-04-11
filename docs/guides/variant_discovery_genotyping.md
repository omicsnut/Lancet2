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
| Match | **0** | +2 | Prevents 8-bit AVX2 SIMD overflow (see below) |
| Mismatch | −6 | −4 | Moderate penalty keeps MSA intact through dense MNVs |
| Gap Open 1 | −6 | — | Short-gap model: penalizes sequencer noise micro-indels |
| Gap Extend 1 | −2 | — | Short-gap model: moderate extension cost |
| Gap Open 2 | −26 | — | Long-gap model: high upfront cost, but... |
| Gap Extend 2 | −1 | — | ...cheap per-base extension for large biological indels |

#### Why Match = 0 (SIMD Overflow Prevention)

SPOA uses **8-bit signed integer AVX2 SIMD lanes** for the DP matrix computation. With the classical match score of +2, a 1,000 bp window accumulating match bonuses would reach a score of 2,000 — severely overflowing the signed 8-bit range (max 127). By anchoring the match score at 0, all runtime scores accumulate negatively, staying within the SIMD lane boundaries.

All other parameters are mathematically shifted downward by −2 to compensate: a standard +2/−4 match/mismatch scheme becomes 0/−6.

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
| Z-Drop | 100,000 | 400 | Effectively disables DP truncation. Prevents alignment termination across large somatic deletions where the gap penalty exceeds the standard Z-drop threshold. |
| Bandwidth (`bw`) | 10,000 | 500 | Envelopes insertions up to ~2 kbp. Prevents the banding boundary from terminating alignment within large assembled insertions. |
| Seed k/w | 11/5 | 15/10 | Increases sensitivity for highly mutated fragments. Maps reads through dense mutation clusters and STRs where 15 bp exact matches are rare. |
| Gap model | Single-affine (12/3) | Dual-affine | In Phase 2, gaps are noise (not biology). A single strict affine model cleanly penalizes all gaps without the "cheap extension" loophole. |
| `best_n` | 1 | — | Keeps only the single best hit per haplotype. |

Since alignment is restricted to the local contig window (~1 kbp), these inflated parameters have minimal runtime impact compared to whole-genome alignment.

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
| `local_raw_score` | Raw substitution matrix score within the variant's physical boundaries on the haplotype | **Subtracted** to carve an algebraic "hole" and prevent double-counting when the PBQ score is added |
| `sc_penalty` | Soft-clip penalty (clipped bases × mismatch cost) | **Subtracted** to crush chimeric supplementary mappings that only partially align |
| `local_pbq_score` | PBQ-weighted DP score within the variant region: each base's substitution matrix contribution is scaled by Phred confidence `(1 − 10^(−Q/10))` | High-confidence bases contribute more; Q5 bases are aggressively downweighted |
| `local_identity` | Fraction of exact matches within the variant region from the CIGAR trace | **Gates** the PBQ score — a high PBQ score from a fragmented, low-identity alignment (e.g., STR artifact) is mathematically discounted |

#### Why Subtract `local_raw_score`?

The `global_score` from minimap2 already includes the raw alignment cost spanning the variant region. If `local_pbq_score` were naively added on top, the variant region would be **double-counted**. By subtracting `local_raw_score`, the formula carves an exact algebraic hole in the global alignment, allowing the higher-fidelity PBQ-weighted local score to replace it cleanly.

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

After allele assignment, genotype likelihoods are computed using the standard per-read model:

```
P(read | GT) = 0.5 × P(read | a₁) + 0.5 × P(read | a₂)
```

across all diploid genotypes, producing Phred-scaled likelihoods (PL). Genotype quality (GQ) is the second-lowest PL value, capped at 99.

* **Read more:** [Alignment-Derived Annotations](alignment_annotations.md), [VCF Output Reference](vcf_output.md)
