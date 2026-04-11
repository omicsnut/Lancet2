# Polar Coordinate Features

These per-sample FORMAT fields transform Cartesian allele depth coordinates
(`AD_Ref`, `AD_Alt`) into polar coordinates, separating biological
identity from sequencing depth so that neither depends on the other.

---

## Why Not DP and VAF?

Standard VCF metrics already encode depth (`DP = AD_Ref + AD_Alt`) and allele
fraction (`VAF = AD_Alt / DP`). At first glance, these seem sufficient. The
problem is geometric: DP and VAF have structural properties that degrade ML
classification, and the polar transform fixes both simultaneously.

### DP vs. PRAD: How We Measure "Weight of Evidence"

**DP** is a simple sum (the L₁ norm) of the allele depth vector.
**PRAD** is the log-compressed Euclidean distance (the L₂ norm) —
the magnitude of the evidence *vector*.

The difference matters because DP **couples** the two allele counts: for a
fixed DP, increasing one allele's count forces the other's count down. This
treats REF and ALT as interchangeable, which they are not — a somatic variant
with AD=(95, 5) and a germline het with AD=(50, 50) both have DP=100, but
represent fundamentally different signal geometries. PRAD treats REF and ALT
as **independent dimensions** of a vector, so the distance from the origin
(the "zero-information" state) reflects the true signal magnitude.

Crucially, sequencing noise scatters reads **equally in all directions**
around the true allele depth point — a read is just as likely to land 2 counts
above as 2 counts to the right. This means the zone of expected noise forms a
roughly circular region, not a diagonal band. The Euclidean radius (PRAD)
traces a circular arc that matches this noise shape naturally. A constant DP
traces a straight diagonal that cuts through the noise at an angle,
artificially coupling the two allele counts.

### VAF vs. PANG: The Unstable Precision Problem

**VAF** is a linear ratio: `Alt / (Ref + Alt)`. **PANG** is an angular
ratio: `atan2(Alt, Ref)`.

VAF has **depth-dependent precision** — its noise level changes with coverage.
A single-read fluctuation at 10× depth swings VAF by 10%; the same
fluctuation at 1000× barely moves it (0.1%). This means two VAF=0.5 data
points at different depths carry radically different *reliability*, but
the number itself gives no indication of this. An ML model using raw VAF
must learn the "fan" (wedge) shape where the meaning of a VAF value
*depends on DP* — a hidden dependency that forces the model to learn the
interaction between VAF and DP rather than treating each independently.

PANG computes the angle of the allele depth vector, which is a pure ratio
independent of magnitude. When paired with PRAD, the two make the noise spread
**uniform** across the feature space: the flaring Cartesian wedge (where
variance grows with depth) unrolls into a rectangular strip with constant-width
noise at every depth. PANG answers **"what is it?"** (somatic vs. germline),
PRAD answers **"how confident are we?"**, and these two questions become
**independent** — allowing ML models to learn each axis separately.

### The Diagonal Collapse

A germline het at 20× sits at AD=(10, 10) in Cartesian space; the same
variant at 2000× sits at AD=(1000, 1000). A distance-based classifier sees
these as enormously far apart despite being biologically identical.

- **With DP/VAF:** the model sees `VAF=0.5, DP=20` and `VAF=0.5, DP=2000` —
  same VAF but wildly different DP. It must learn that "VAF=0.5 means
  germline" regardless of DP, which requires seeing examples at every depth.
- **With PANG/PRAD:** both points have **PANG=0.785 (45°)**. The biological
  identity is a constant. Only the radius differs, and it does so on a
  bounded log scale (1.15 → 3.15). The entire diagonal collapses into a
  single PANG value, and the model only needs to learn that PANG ≈ π/4
  means "germline."

This is the fundamental advantage: the polar transform **algebraically
separates** what a variant *is* from how much evidence supports it, enabling
ML models trained at one coverage to generalize to any other.

The polar transform "unrolls" the wedge and "collapses" the diagonal into
a uniform rectangular strip where ML models can draw simple flat boundaries.

---

## Polar Radius (`PRAD`)

**In plain terms**: PRAD answers "how much total evidence do I have?" on a
compressed scale where 1.0 ≈ ~10 reads, 2.0 ≈ ~100 reads, 3.0 ≈ ~1000 reads.

**Computation**: `PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²))`, computed per-sample.

The log10-compressed Euclidean magnitude of the allele depth vector. Encodes
total signal strength (sequencing depth) independent of allele fraction, with
built-in normalization for cross-coverage ML generalization.

**Why log10?** The raw Euclidean radius `sqrt(AD_Ref² + AD_Alt²)` scales
linearly with coverage, spanning 0 to >2800 at 2000× depth — a 100× dynamic
range. An ML model trained at 30× WGS would see completely different PRAD
distributions when applied to 100× WES data. The `log10(1+r)` compression
reduces this to a compact [0, ~3.5] range while preserving monotonicity:

| Coverage | Het 50% VAF | PRAD (log10) | Interpretation |
|:---------|:------------|:-------------|:---------------|
| 0× | — | 0.0 | No reads — no evidence |
| 20× | AD=(10,10) | 1.15 | Low confidence |
| 60× | AD=(30,30) | 1.63 | Moderate confidence |
| 100× | AD=(50,50) | 1.85 | Good confidence |
| 1000× | AD=(500,500) | 2.85 | Very high confidence |
| 2000× | AD=(1000,1000) | 3.15 | Extreme confidence |

**Value range**: [0, ~3.5]

**Interpretation**: Higher PRAD means more sequencing evidence, regardless
of whether it supports REF or ALT. Acts as a confidence tie-breaker for
low-angle variants:

| PRAD | Meaning |
|:-----|:--------|
| < 1.0 + low PANG | Stochastic noise — insufficient evidence |
| < 1.0 + high PANG | Low-coverage heterozygous/homozygous call |
| > 2.0 + low PANG | True low-frequency variant (somatic/mosaic) with deep sequencing support |
| > 2.0 + high PANG | Well-supported heterozygous/homozygous variant |

**Coverage stability**: PRAD still varies with coverage (1.15 at 20× to 3.15
at 2000×), but the log10 compression reduces the dynamic range from 100× to
3× and preserves monotonic ordering. For ML models that use PRAD as a
confidence weight, this bounded range avoids feature-scale domination that
occurs with raw allele counts. Tree-based models (XGBoost, Random Forest)
are invariant to monotonic transforms and handle this naturally; linear
models benefit from the compressed range.

---

## Polar Angle (`PANG`)

**In plain terms**: PANG answers "what fraction of reads support the variant?"
expressed as an angle. The key property: this angle is the same whether you
have 10 reads or 10,000 — depth doesn't change it.

**Computation**: `PANG = atan2(AD_Alt, AD_Ref)`, computed per-sample.

The arctangent of the allele depth ratio, in radians. Encodes the biological
identity of the variant — the allele fraction — independent of total depth.

**Value range**: [0, π/2] radians (≈ [0, 1.571])

**Interpretation**:

| PANG (radians) | PANG (degrees) | Meaning |
|:---------------|:---------------|:--------|
| ≈ 0.0 | 0° | Homozygous REF — no ALT evidence |
| < 0.524 | < 30° | Low-frequency signal (somatic, mosaic, or artifact) |
| ≈ 0.785 | 45° | Heterozygous germline — balanced REF/ALT |
| ≈ 1.571 | 90° | Homozygous ALT — nearly all reads are ALT |

**Coverage stability**: Perfectly coverage-invariant. PANG is a pure ratio
(`atan2(AD_Alt, AD_Ref)`) — identical allele fractions produce identical
angles at any depth. A 50% VAF het produces PANG ≈ 0.785 (45°) whether
the sample is 20× or 2000×. This is the key advantage of the polar transform:
allele identity is algebraically separated from sequencing depth.

---

## Tumor-Normal "Four-Feature" Paradigm

Each sample independently computes its own PRAD and PANG from its own allele
depths. The downstream ML model receives a four-element feature vector:

`[PRAD_Normal, PANG_Normal, PRAD_Tumor, PANG_Tumor]`

Both PRAD (log10-compressed, range [0, ~3.5]) and PANG (ratio-based, range
[0, π/2]) are coverage-invariant, ensuring the four-feature paradigm produces
comparable values regardless of whether the sample is 20× or 2000×.

This avoids cross-sample relational metrics that would violate VCF FORMAT
semantics and eliminates sentinel/NaN imputation issues. The model learns
the relational delta itself:

| Pattern | PANG Normal | PANG Tumor | Classification |
|:--------|:------------|:-----------|:---------------|
| Germline | ≈ π/4 (45°) | ≈ π/4 (45°) | Inherited variant — same allele fraction in both |
| Somatic | ≈ 0 (0°) | ≈ π/6..π/4 | Tumor-specific variant — normal has no ALT support |
| LOH | ≈ π/4 (45°) | ≈ π/2 (90°) | Loss of heterozygosity — tumor lost REF allele |
