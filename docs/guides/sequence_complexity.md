# Sequence Complexity (`SEQ_CX`)

This INFO field captures 11 sequence complexity features distilled from
four raw measurement modules (homopolymer runs, Shannon entropy,
LongdustQ k-mer concentration, tandem repeat motifs). Emitted as the `SEQ_CX`
INFO tag in every VCF record.

All 11 features are computed from assembled haplotype sequences and are
**completely coverage-invariant**: they measure properties of the genomic
sequence itself, not properties of the reads. A poly-A run is 12bp long
regardless of whether 20 or 2000 reads cover it.

**INFO tag**: `SEQ_CX` (11 comma-separated values)

---

## The Four Measurement Modules

The 11 distilled features are built from four raw measurement modules, each
capturing a different aspect of sequence complexity at different biological
scales. Understanding these building blocks makes the distilled features
intuitive.

### 1. Homopolymer Run (HRun)

The longest consecutive stretch of identical bases in a DNA window.

**What it measures**: The length of the single longest same-base run — e.g.,
`AATTTTTTTTGC` has an HRun of 8 (the T run).

**Why it matters**: Homopolymers are the #1 source of Illumina INDEL errors.
During sequencing, the polymerase can "stutter" on long same-base runs,
adding or dropping 1–2 bases. A 12bp poly-A run produces ~10× more false
1bp insertions than a random 12bp sequence. HRun directly quantifies this
stutter risk.

**Computation**: Single O(L) scan counting consecutive identical characters.
Returns 0 for empty sequences, 1 for single-base sequences. No parameters.

### 2. Shannon Entropy

Base frequency diversity in a DNA window, measured in bits (0.0–2.0).

**What it measures**: How evenly the four DNA bases (A, C, G, T) are
distributed. The formula `H = −Σ pᵢ log₂(pᵢ)` counts the information
content of the base composition.

**Interpretation**:

- **0.0** = all one base (e.g., `AAAAAAA`) — maximally repetitive
- **1.0** = two equally frequent bases (e.g., alternating `ACACAC`)
- **2.0** = perfectly balanced ACGT — maximally diverse

**Why it matters**: Low-entropy regions are harder to align uniquely because
many positions look identical to the aligner. Variants called in low-entropy
windows are more likely to be mismapping artifacts. Shannon entropy
captures this alignment ambiguity risk.

**Computation**: O(L) base counting, then the entropy formula. No parameters.

### 3. LongdustQ (k-mer Concentration)

A k-mer-based complexity score that measures how repetitive a sequence is
by comparing the observed k-mer count distribution against a random
(Poisson) null model.

**What it measures**: Whether any short DNA words (k-mers) appear more
often than expected by chance. In random DNA, each k-mer appears roughly
equally. In repetitive DNA, a few k-mers dominate — and LongdustQ
quantifies this concentration.

**Interpretation**:

| LongdustQ Score | Meaning |
|:----------------|:--------|
| ≈ 0 | Random / unique sequence |
| > 0.6 | Moderately repetitive (longdust's default LCR threshold) |
| > 1.0 | Highly repetitive (strong tandem repeat signal) |
| > 2.0 | Extremely repetitive (long homopolymer, satellite DNA) |

**Why it matters**: LongdustQ detects repetitive structures that HRun and
entropy miss — microsatellites, tandem repeats, and satellite DNA that cause
systematic assembly errors. A dinucleotide repeat (CACACACA) has high
entropy (1.0) and short HRun (1), but LongdustQ correctly flags it as
repetitive because two k-mers dominate the distribution.

**Provenance**: Adapted from Li, H. "Finding low-complexity filter with
longdust" (2025). Lancet2 extracts just the Q(x) scoring formula and applies
it directly to known subsequences around variants — no sliding-window
scanning is needed because variant locations are already known from graph
assembly. Two k-mer sizes are used: **k=4** for local ±50bp flanks (sensitive
to short repeats) and **k=7** for full haplotype scoring (sensitive to
macro-scale satellite structure).

**GC-bias correction**: LongdustQ supports an optional GC-bias correction
(default: human genome GC=0.41) to prevent AT-rich regions from being
falsely scored as repetitive. See [GC-Bias Correction](#gc-bias-correction).

### 4. Tandem Repeat (TR) Motif Detection

Finds exact and approximate short tandem repeats (period 1–6) in a flanking
window around the variant.

**What it measures**: The presence, proximity, purity, and period of the
nearest tandem repeat to the variant site. Searches for motifs that repeat
≥2.5× (exact) or ≥3.0× (approximate, allowing 1 edit per repeat unit).

**Why it matters**: Repeat proximity is the single strongest predictor of
Illumina INDEL artifacts. A 1bp insertion adjacent to a (CA)₁₀ dinucleotide
repeat is overwhelmingly likely to be polymerase stutter, not a real
mutation. TR detection directly identifies this pattern and flags the
canonical stutter signature.

**Key outputs**:

- **Distance**: bp distance from variant to nearest TR (0 = overlapping)
- **Purity**: fraction of the repeat span that is error-free (1.0 = perfect)
- **Period**: motif length (1=homopolymer, 2=dinucleotide, ..., 6=hexanucleotide)
- **Stutter flag**: 1 if INDEL length ≤ period AND adjacent to TR (≤1bp)

**Computation**: O(L × P) where L = window length and P = max period (6).
See [Motif Detection Parameters](#motif-detection-parameters) for defaults.

---

## Design: Context vs. Perturbation Paradigm

Raw sequence complexity metrics at multiple scales (±5/10/20/50/100bp, full
haplotype) are highly correlated — HRun at ±5bp, ±10bp, and ±20bp are ~0.95
correlated. Additive ML models — such as EBMs (Explainable Boosting Machines,
which learn one feature's effect at a time, then combine them) and GAMs
(Generalized Additive Models) — split importance weight across
correlated features, producing noisy shape functions instead of confident ones.

The distillation separates metrics into three groups that are independent
**between** groups (Context, Delta, TR Motif measure fundamentally different
things), while features **within** each group are deliberately chosen at
non-overlapping scales to minimize correlation:

1. **Context** (REF-only): "How brittle is the genome here, regardless of the variant?"
2. **Deltas** (ALT−REF): "How did the variant alter the local complexity?"
3. **TR Motif** (ALT-only): "What is the repeat environment of the final allele?"

This enables the model to learn rescue logic: even if context is dangerous
(high ContextHRun), a negative DeltaHRun (variant broke the homopolymer)
rescues the call.

---

## Feature Layout

### A. Context (4 features — strictly REF)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 0 | ContextHRun | Integer | [0, ∞) | Max homopolymer run in REF ±20bp. A poly-A of length 12 means any 1bp INDEL within 20bp is a stutter candidate. |
| 1 | ContextEntropy | Float | [0.0, 2.0] | Shannon entropy H = −Σ pᵢ log₂(pᵢ) in REF ±20bp. Measures sequence diversity: 0 = all one base (like AAAAAAA), 2.0 = perfectly balanced ACGT. Low entropy = harder to align uniquely. |
| 2 | ContextFlankLQ | Float | [0.0, ~1.6] | log₁p-squashed LongdustQ (k=4) in REF ±50bp. LongdustQ is a k-mer concentration score: it counts how many short DNA words (k-mers) repeat more often than expected by chance. Higher values = more repetitive DNA = more likely to cause sequencing and assembly errors. The log₁p transform compresses heavy tails (telomeric LQ > 4.0 → ~1.6). |
| 3 | ContextHaplotypeLQ | Float | [0.0, ~1.6] | log₁p-squashed LongdustQ (k=7) on full REF haplotype. Captures macro-scale repetitive structure of the assembled contig. Only REF because a 5bp INDEL changes < 0.5% of a 1000bp k-mer distribution. |

### B. Deltas (3 features — ALT minus REF)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 4 | DeltaHRun | Integer | [−∞, +∞) | ALT ±5bp HRun − REF ±5bp HRun. Positive = variant extended a homopolymer (artifact signal). Negative = variant broke a homopolymer (rescue signal). Computed at ±5bp to maximize sensitivity to immediate allelic change. |
| 5 | DeltaEntropy | Float | [−2.0, +2.0] | ALT ±10bp entropy − REF ±10bp entropy. Negative = variant reduced sequence diversity (gap mimicking deletion). Computed at ±10bp to balance between ultra-local and micro sensitivity. |
| 6 | DeltaFlankLQ | Float | [−1.6, +1.6] | log₁p(ALT ±50bp LQ) − log₁p(REF ±50bp LQ). Positive = variant exacerbated microsatellite repetitiveness. Computed in log-space to normalize the heavy-tailed LQ distribution. |

### C. TR Motif (4 features — strictly ALT ±50bp)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 7 | TrAffinity | Float | [0.0, 1.0] | 1/(1+dist) where dist = distance to nearest tandem repeat in ALT ±50bp. Sentinel-safe: dist < 0 (no TR found) → 0.0. dist = 0 (sitting on TR) → 1.0. Monotonically decaying with distance. |
| 8 | TrPurity | Float | [0.0, 1.0] | Purity (1 − errors/span) of nearest TR. dist < 0 → 0.0. A purity of 0.95 means near-perfect repeat fidelity, strongly predicting polymerase stutter. |
| 9 | TrPeriod | Integer | [0, 6] | Period of nearest TR. dist < 0 → 0. Period dictates the biophysics of slippage: period 1 (homopolymer) has highest stutter rate, period 6 (hexanucleotide) has lowest. |
| 10 | IsStutterIndel | Integer | {0, 1} | 1 if INDEL length ≤ repeat period AND variant is adjacent to a TR (dist ≤ 1bp). This is the canonical Illumina polymerase stutter signature. |

---

## Transform Details

### Sentinel-Safe TR Features

Raw motif detection uses `dist_to_nearest_tr = −1` as a sentinel for "no TR
found." This breaks EBM numerical binning because −1 sorts between −2 and 0,
making "no TR" appear closer to "sitting on the repeat" than a distance of 10.

The `TrAffinity` inverse transform (`1/(1+dist)`) maps the sentinel into the
continuous [0, 1] range with 0 = no TR and 1 = overlapping TR. All other TR
features are zeroed when no TR is found.

### Log₁p Squashing

LongdustQ reaches 4.0+ in telomeric and pericentromeric regions, creating
extreme heavy tails. The `log₁p(max(0, x))` transform compresses these tails
while preserving monotonicity and differentiability at zero.

### Multi-Haplotype and Multi-Allelic Merging

When a variant appears on multiple ALT haplotypes, the scorer produces one
`SequenceComplexity` per haplotype and merges via element-wise max (pessimistic
worst-case). Multi-allelic variants at the same locus are similarly merged.

---

## GC-Bias Correction

LongdustQ uses a null model for expected k-mer frequencies. By default, it
assumes GC fraction = 0.41 (human genome-wide average). Without correction,
AT-rich regions are scored as artificially "repetitive."

```bash
--genome-gc-bias 0.41    # Human genome (default)
--genome-gc-bias 0.5     # Uniform model (no correction)
--genome-gc-bias 0.42    # Mouse genome
```

---

## Coverage Stability

All 11 SEQ_CX features are **perfectly coverage-invariant** because they
are computed from assembled haplotype sequences, not from read-level metrics.
The scoring pipeline works as follows:

1. The de Bruijn graph assembler produces haplotype sequences (contigs)
   from the reads in each window.
2. The sequence complexity scorer operates on these haplotype strings.
3. Homopolymer runs, Shannon entropy, LongdustQ, and TR motif detection
   all operate on the sequence characters, not on read counts or qualities.

The only indirect coverage dependency is that the assembler itself requires
sufficient coverage to produce correct haplotypes. Below ~10×, the
assembler may fail to construct the ALT haplotype, in which case no variant
is called and no SEQ_CX is computed. Above ~15×, the assembled haplotypes
are identical regardless of depth, and all 11 features are deterministic:

| Feature | Coverage dependency |
|:--------|:--------------------|
| ContextHRun | None — counts bases in REF string |
| ContextEntropy | None — base frequency ratio in REF string |
| ContextFlankLQ | None — k-mer concentration in REF string |
| ContextHaplotypeLQ | None — k-mer concentration in full REF haplotype |
| DeltaHRun | None — ALT string vs REF string |
| DeltaEntropy | None — ALT string vs REF string |
| DeltaFlankLQ | None — ALT string vs REF string |
| TrAffinity | None — motif distance in ALT string |
| TrPurity | None — motif purity in ALT string |
| TrPeriod | None — motif period in ALT string |
| IsStutterIndel | None — binary flag from ALT string |

---

## Motif Detection Parameters

The tandem repeat motif detector scans for exact and approximate repeats.
Default configuration is optimized for somatic variant calling on Illumina
short reads:

| Parameter | Default | Rationale |
|:----------|:--------|:----------|
| `max_period` | 6 | Covers homopolymers through hexanucleotide repeats. Period 1–3 account for > 90% of Illumina INDEL errors. |
| `min_copies_exact` | 2.5 | Catches ATATAT (3× AT) but not ATAT (2× AT). |
| `min_copies_approx` | 3.0 | Higher threshold for approximate matches to compensate for relaxed matching. |
| `max_edits_per_unit` | 1 | Allows 1 mismatch/indel per repeat unit. |
| `min_purity` | 0.75 | Below 75% purity, the repeat is too degraded to cause polymerase stutter. |
