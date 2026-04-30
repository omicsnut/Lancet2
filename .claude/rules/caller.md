---
description: Lancet2 caller/ layer rules — SPOA convex scoring with int16 SIMD lane discipline, minimap2 read-to-haplotype realignment, Dirichlet-Multinomial genotype likelihoods, allele assignment scoring formula. Load when editing src/lancet/caller/.
paths:
  - "src/lancet/caller/**"
---

# caller/ layer rules

`caller/` turns assembled haplotypes into VCF records. The data flow
is: assembled paths from `cbdg/` → MSA via SPOA → variant extraction
→ per-read allele assignment via minimap2 → genotype likelihoods →
VCF record. Each step has tight constraints that aren't obvious from
the API surface.

## SPOA scoring constants are perf-correctness-coupled

`msa_builder.h` declares six scoring constants that **must not be
changed without understanding all their effects together**:

```cpp
MSA_MATCH_SCORE     =   0
MSA_MISMATCH_SCORE  =  -6
MSA_OPEN1_SCORE     =  -6   MSA_EXTEND1_SCORE  =  -2
MSA_OPEN2_SCORE     = -26   MSA_EXTEND2_SCORE  =  -1
```

Three independent constraints are encoded:

1. **Match=0 keeps SIMD in int16 lane width.** SPOA 4.1.5 dispatches
   between int16 and int32 paths via `WorstCaseAlignmentScore()`.
   Match=0 keeps all runtime scores non-positive, which keeps the
   worst-case score (~−2300 for typical alignments) above int16's
   −32768 minimum. int16 SIMD has 16 lanes per AVX2 register; int32
   has 8. Setting Match=+1 silently halves throughput.

2. **Convex (dual-affine) gap scoring is biological, not aesthetic.**
   Linear gaps overpenalize multi-bp events. Single-affine forces a
   tradeoff between small-variant strictness and large-indel
   tolerance. Convex takes `min(g1 + (i−1)·e1, g2 + (i−1)·e2)` —
   strict for short gaps (suppresses sequencer noise), cheap for
   large gaps (admits real structural variants). The two models
   intersect at exactly 20bp where `−6 + 2L = −26 + 1L`.

3. **Mismatch=−6, not −19 or −4.** The asm5 minimap preset uses −19
   which shatters dense MNV clusters; classic +2/−4 lets sequencer
   noise into the MSA. −6 globally forces alignments through complex
   variants without admitting noise.

Changing any constant requires re-justifying all three constraints.
Don't.

## Per-thread spoa::AlignmentEngine is required (and reused)

`spoa::AlignmentEngine` is not thread-safe across threads, and
constructing it per call is performance-fatal. The engine is created
once per `MsaBuilder` instance and reused across all of a window's
SPOA calls; `MsaBuilder` itself is owned per-thread by
`VariantBuilder`. New code that creates a fresh `AlignmentEngine`
per call is a correctness-safe but throughput-fatal regression.

## Genotype likelihoods: Dirichlet-Multinomial, lgamma-stable

`genotype_likelihood.cpp::ComputeGenotypePLs` uses a Dirichlet-
Multinomial (DM) model, not per-read binomial. Why DM:

1. Generalizes to K alleles (multi-allelic sites), not just K=2.
2. Models correlated errors at ultra-high depth via overdispersion ρ.
3. Numerically stable via `std::lgamma`, not `lgamma(...) ≈ ...`
   approximations.

The closed-form is `ln P(c | μ, M) = lnΓ(ΣMμᵢ) − lnΓ(N + ΣMμᵢ) +
Σ[lnΓ(cᵢ + Mμᵢ) − lnΓ(Mμᵢ)]`. The constants `DM_BACKGROUND_ERROR =
0.005`, `DM_OVERDISPERSION = 0.01`, `DM_ALPHA_FLOOR = 1e-6` are
tuned for Illumina short-read profiles. Each constant has a
documented effect on PL behavior at different coverage regimes (see
`genotype_likelihood.cpp` constants block); don't change in
isolation.

## PL ordering is VCF-standard unphased: (0,0)(0,1)(1,1)(0,2)(1,2)(2,2)...

For K alleles, `ComputeGenotypePLs` returns `K·(K+1)/2` PLs in
**VCF-standard unphased ordering** — the j-th allele pair is at
index `j·(j+1)/2 + i` for genotype `(i,j)` with `i ≤ j`. This is the
order downstream consumers (bcftools, GATK, hap.py) expect. Code
that emits `(0,0)(0,1)(0,2)(1,1)(1,2)(2,2)` (lex-by-allele) is
silently wrong — every consumer reads it as a different genotype.

Best PL is normalized to 0; others scale relative. **GQ = min2 −
min1** (second-smallest minus smallest), capped at 99. This is the
GATK convention; don't replace with raw min2 (which double-counts
the offset).

## Allele assignment scoring formula has subtractive structure

The combined score for a read-haplotype pair is:

```
combined = (global_score − local_raw_score − sc_penalty)
         + (local_pbq_score · local_identity)
```

Each term has a specific role:

- `global_score`: minimap2 DP score of the full read→haplotype
  alignment, including flanks.
- `local_raw_score`: substitution-matrix score within the variant
  overlap, **subtracted** to carve a hole for the PBQ-weighted
  re-scoring (no double-counting).
- `sc_penalty`: soft-clip penalty, suppresses chimeric supplementary
  mappings that minimap exempts from primary DP.
- `local_pbq_score · local_identity`: PBQ-weighted local score
  multiplied by exact-match fraction. The identity gate prevents
  "lucky" alignments in low-complexity regions from winning over
  clean alignments.

Removing any term silently degrades accuracy in a specific class of
reads. The closest existing test coverage is
`tests/caller/variant_support_metrics_test.cpp`, which exercises the
per-allele FORMAT-field computations downstream of this formula.
Direct unit tests for the combined-score terms themselves don't yet
exist; if you change this formula, add coverage there alongside the
change.

## Folded read position: min(p, 1−p), not raw p

`RPCD` (read position Cohen's D) uses **folded** position `min(p,
1−p)`. Artifacts cluster at BOTH read ends (3' quality decay, 5'
soft-clip misalignment). With raw positions, bimodal end-clusters
average to ~0.5, indistinguishable from a centered variant. Folding
maps both ends to the same low-value space, converting the bimodal
trap into a unidirectional "ALT closer to edges than REF?" signal.

## Edit distance (NM) is always vs the REFERENCE haplotype

The NM field used in ASMD computation is always the read's edit
distance against the **reference haplotype (hap_idx=0)**, regardless
of which allele the read is assigned to. This is so `ASMD = mean(ALT
NM) − mean(REF NM)` cancels the variant's own contribution and
isolates excess noise. Code that uses NM-vs-assigned-haplotype
silently breaks ASMD's interpretation.

The CIGAR-vs-encoded-REF computation excludes soft clips, hard
clips, and N-skips per SAM spec. Including soft clips here double-
counts the alignment endpoint flexibility.

## Representative base quality for indels: minimum, not mean

For SNVs, the read's representative base quality at the variant is
the single base quality at that position. **For indels, it's the
minimum base quality across all read positions spanning the variant
region** (weakest-link summary). bcftools mpileup uses the same
convention; GATK PairHMM integrates differently. The minimum is
intentional: a 10bp deletion where 9 bases are Q30 and 1 is Q5 is
not a high-confidence observation. Code that uses mean here lets
single low-Q bases hide in long indels.

## CMLOD is per-ALT, REF index = 0.0

`ComputeContinuousMixtureLods` returns one score per allele with REF
at index 0 and ALTs at 1..K−1. **REF is always 0.0** — CMLOD
measures evidence of an ALT existing at its observed frequency vs.
being absent; the REF "doesn't exist" hypothesis is meaningless.
Code that emits a non-zero CMLOD for REF is a schema bug.

The output is `Number=A` in VCF (per-ALT, REF excluded). See
`vcf-validator` for the schema invariants this depends on.
