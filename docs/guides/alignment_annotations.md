# Alignment-Derived Annotations

These per-sample FORMAT fields detect alignment, sequencing, and library
preparation artifacts. All fields are designed to be **coverage-invariant**:
the same biological signal produces the same annotation value regardless of
sequencing depth, enabling ML models trained at one coverage to generalize
across 20×–2000× without retraining.

SCA, FLD, and MQCD use metadata from the original BAM/CRAM alignment.
RPCD, BQCD, and ASMD use metrics computed from the minimap2 re-alignment
during genotyping. SDFC uses window-level BAM coverage. FSSE and AHDD use
per-read alignment evidence from the genotyper. HSE uses SPOA haplotype
assignment data. PDCV uses k-mer coverage from the de Bruijn graph paths.

---

## Soft Clip Asymmetry (`SCA`)

**Purpose**: Detect unresolved larger structural variant events that
masquerade as smaller local variants in the assembly.

**Computation**:

1. During BAM/CRAM ingestion, each read is flagged as "soft-clipped" if
   the total soft-clip bases in its original whole-genome alignment CIGAR
   are ≥ 6% of the read length.
2. After genotyping assigns reads to REF or ALT alleles, the soft-clip
   fraction for each allele is computed:
   - `alt_frac = alt_soft_clip_count / alt_total_count` (0 if no ALT reads)
   - `ref_frac = ref_soft_clip_count / ref_total_count` (0 if no REF reads)
3. `SCA = alt_frac − ref_frac`, computed per-sample.

**Value range**: [−1.0, 1.0]

**Interpretation**:

| SCA Range | Meaning |
|:----------|:--------|
| ≈ 0.0 | Symmetric soft-clipping between ALT and REF — likely no hidden event |
| > 0.1 | ALT reads disproportionately soft-clipped — possible unresolved SV, translocation, or complex event |
| < −0.1 | REF reads more clipped — unusual, may indicate mapping artifact near the reference path |

**Coverage stability**: Inherently coverage-invariant — computed as a ratio
of fractions. A 30% soft-clip rate in ALT reads produces SCA ≈ 0.3 at any
depth. Using SCA together with SDFC helps distinguish real SVs (high SCA,
typical SDFC) from mapping artifacts (high SCA, elevated SDFC).

---

## Fragment Length Delta (`FLD`)

**Purpose**: Detect chimeric library artifacts (artificial bridging) or
somatic cfDNA fragment length shifts.

**Computation**:

1. During BAM/CRAM ingestion, each read's template length (`TLEN` /
   `bam1_t::core.isize`) is captured from the original alignment.
2. After genotyping, for each sample, properly-paired reads with non-zero
   insert size are grouped by their assigned allele.
3. Mean insert sizes are computed for REF-supporting and ALT-supporting
   reads separately.
4. `FLD = mean_alt_isize − mean_ref_isize`, computed per-sample.
   The sign is preserved: negative values indicate shorter ALT fragments
   (common in cfDNA tumor-derived fragments). If either group has no
   properly-paired reads, the metric is **untestable** and emitted as
   `.` (VCF 4.5 missing value).

**Value range**: (−∞, +∞), or `.` (untestable)

**Interpretation**:

| FLD Range | Meaning |
|:----------|:--------|
| < 20 bp | Expected variation — insert size distributions are consistent |
| 20–100 bp | Moderate discrepancy — inspect manually |
| > 100 bp | Large discrepancy — likely chimeric artifact, library prep issue, or cfDNA fragment size shift |

**Coverage stability**: The mean insert size converges rapidly (by ~20
reads per allele). FLD shows minor variation at very low coverage (N<5 per
allele) due to sampling noise, but is effectively stable above 20×.
Interpret FLD jointly with RPCD: a true structural variant often produces
both elevated FLD and edge-biased RPCD.

---

## Mapping Quality Cohen's D (`MQCD`)

**Purpose**: Detect paralogous mismapping — situations where ALT-supporting
reads originate from a different genomic locus (e.g., a segmental duplication)
and are assigned artificially low mapping quality by the aligner.

**Statistical method**: Coverage-normalized effect size from the Mann-Whitney
U test (Wilcoxon Rank-Sum test). The raw Z-score is divided by √N (where
N = total reads) to remove the √N power amplification from the Central Limit
Theorem, recovering a standardized effect size analogous to Cohen's d.

**In plain terms**: MQCD measures "are the variant-supporting reads less
confidently mapped than the reference-supporting reads?" If ALT reads
consistently have lower mapping quality, the variant may be a mismapping
artifact rather than a real mutation. The score is designed so that the same
biological bias produces the same number whether you have 20 reads or 20,000
— making it safe to use across different sequencing depths.

**Motivation for Z/√N normalization**: The raw Mann-Whitney Z-score scales
with √N: the same mild ALT MAPQ depression produces Z ≈ −1.5 at 20× but
Z ≈ −14.9 at 2000×. This makes raw Z-scores unusable for ML models that
must generalize across coverages. Dividing by √N produces a coverage-invariant
effect size that measures the *magnitude* of the MAPQ difference, not the
statistical significance.

**Computation**:

1. Collect original alignment MAPQ (`bam1_t::core.qual`) for all reads,
   grouped by their genotyper-assigned allele (REF vs ALT).
2. Pool all observations into a single ranked sequence. Tied values
   receive the mean (mid-rank) of the positions they span.
3. Compute the U statistic: `U = R_alt − n_alt·(n_alt+1)/2`.
4. Apply the tie-corrected variance formula (Lehmann, 2006):
   `Var(U) = (m·n/12) · [(N+1) − Σ(tₖ³−tₖ)/(N·(N−1))]`
5. `MQCD = Z / √N = [(U − E[U]) / √Var(U)] / √(m+n)`, computed per-sample.

**Value range**: [−2, +2] typically, or `.` (untestable — one or both groups empty)

**Coverage stability**: A constant mild bias (ALT MAPQ 2 units lower)
produces MQCD ≈ −0.34 at **every** depth from 20× to 2000× (< 3% variation).

**Interpretation**:

| MQCD Range | Meaning |
|:-----------|:--------|
| `.` | Untestable: one or both allele groups empty (no reads to compare). |
| ≈ 0.0 | No systematic MAPQ difference — strong evidence for true variant. |
| −0.2 to −0.5 | Moderate ALT MAPQ depression — possible repetitive region, inspect manually |
| < −0.5 | Strong ALT mismapping signal — likely paralogous or multi-mapping artifact |
| > 0 | ALT MAPQ higher than REF — unusual, may indicate REF mapping issues |

**Manual interpretation tip**: Use MQCD together with SDFC. Paralogous
mismapping typically produces both negative MQCD (low ALT MAPQ) and
elevated SDFC (excess depth from collapsed paralogs). A site with
MQCD < −0.3 and SDFC > 2.0 is very likely a false positive from
segmental duplication.

**References**:

- Mann, H.B. & Whitney, D.R. (1947). *Annals of Mathematical Statistics*, 18(1), 50–60.
- Lehmann, E.L. (2006). *Nonparametrics: Statistical Methods Based on Ranks*, Springer.
- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*, 2nd ed.

---

## Read Position Cohen's D (`RPCD`)

**Purpose**: Detect systematic read-edge bias in ALT-supporting reads.
Artifacts from 3' quality degradation and 5' soft-clip misalignment produce
false variant calls that cluster at read extremes.

**Key insight — folded read position**: The test uses the **folded** read
position rather than the raw position. Raw positions are bimodal for
edge-biased artifacts (clustered at both 5' and 3' ends), but the mean
of positions 5 and 145 in a 150bp read ≈ 75 — indistinguishable from a
truly centered variant. Folding maps both ends to the same low-value space:

`P_folded = min(p, 1 − p)` where `p = variant_query_pos / read_length`

**In plain terms**: folding maps both ends of the read to the same "edge
zone." A variant seen at position 5 (near the start) and position 145
(near the end of a 150bp read) both get folded to distance 5 from the
nearest edge. This way, edge-biased artifacts — which tend to cluster
at *either* end — all show up as "close to the edge" rather than
canceling out as "average position = middle."

This converts the bimodal trap into a unidirectional signal: "Are ALT alleles
systematically closer to read edges than REF alleles?"

**Statistical method**: Same Z/√N effect-size normalization as MQCD, applied
to folded read positions instead of mapping qualities.

**Computation**:

1. During minimap2 re-alignment in the genotyper, walk the best-allele CIGAR
   to find the query position corresponding to the variant's haplotype position.
2. Compute folded position: `min(qpos/read_length, 1 - qpos/read_length)`.
   Result: 0.0 = read edge, 0.5 = read center.
3. Group folded positions by allele (REF vs ALT).
4. Apply Mann-Whitney U test → Z/√N effect size, computed per-sample.

**Value range**: [−2, +2] typically, or `.` (untestable — one or both groups empty)

**Coverage stability**: The effect size is coverage-invariant. A consistent
edge bias produces the same RPCD at any depth from 20× to 2000×.

**Interpretation**:

| RPCD Range | Meaning |
|:-----------|:--------|
| `.` | Untestable: one or both allele groups empty (no reads to compare). |
| ≈ 0.0 | Uniform read position distribution — expected for true variants. |
| −0.2 to −0.5 | Moderate edge bias — ALT allele somewhat closer to read ends |
| < −0.5 | Strong edge bias — likely alignment artifact or 3' error cascade |

**Manual interpretation tip**: RPCD should be interpreted jointly with BQCD.
Many artifacts produce both read-edge clustering (negative RPCD) and low
ALT base quality (negative BQCD) because base quality degrades near read
ends. If RPCD < −0.3 and BQCD < −0.3, the variant is very likely an artifact.
True variants may show mild negative RPCD in isolation (e.g., near an indel
that shifts alignment) without the accompanying BQCD depression.

---

## Base Quality Cohen's D (`BQCD`)

**Purpose**: Detect chemistry-driven sequencing artifacts where the ALT allele
is systematically called with lower base confidence than the REF allele.

**Key artifact**: 8-oxoguanine (8-oxoG) oxidation is the dominant source of
G→T / C→A errors in Illumina sequencing. Oxidized guanine mispairs with
adenine, producing consistent low-quality G→T miscalls. The miscalled base
has characteristically low Phred scores that are detectable by this test.

**Statistical method**: Z/√N effect-size normalization on per-read base
qualities at the variant position (REF vs ALT groups).

**Computation**:

1. During genotyping, the representative base quality at the variant position
   is recorded per read (minimum across the variant region for indels).
2. Base qualities are grouped by allele (REF vs ALT), combining forward and
   reverse strand observations.
3. Apply Mann-Whitney U test → Z/√N effect size, computed per-sample.

**Value range**: [−2, +2] typically, or `.` (untestable — one or both groups empty)

**Coverage stability**: Coverage-invariant. A consistent 5-unit ALT quality
depression produces the same BQCD at any depth.

**Interpretation**:

| BQCD Range | Meaning |
|:-----------|:--------|
| `.` | Untestable: one or both allele groups empty (no reads to compare). |
| ≈ 0.0 | No systematic quality difference — expected for true variants. |
| −0.2 to −0.5 | Moderate ALT quality depression — inspect for oxidation or deamination |
| < −0.5 | Strong signal — likely chemistry artifact (8-oxoG, FFPE deamination) |

**Manual interpretation tip**: For targeted 8-oxoG detection, examine BQCD
jointly with the variant allele: G→T and C→A substitutions with BQCD < −0.3
are the classic oxidation signature. Other mutation types with negative BQCD
may indicate FFPE deamination (C→T/U) or other library damage.

---

## Allele-Specific Mismatch Delta (`ASMD`)

**Purpose**: Detect chimeric reads and paralogous mismapping. True variants
should be the only difference between a read and the reference. If ALT-supporting
reads carry excess random mismatches while REF reads are clean, the ALT allele
is likely a misaligned chimera or paralog.

**Computation**:

1. During minimap2 re-alignment, every read is aligned to the **REF haplotype**
   (haplotype index 0), regardless of which allele it is ultimately assigned to.
2. The edit distance (NM) is computed from the REF alignment CIGAR per SAM spec:
   mismatches (under M ops) + insertion bases + deletion bases. Soft clips,
   hard clips, and reference skips are excluded.
3. NM values are grouped by allele assignment (REF vs ALT).
4. The maximum variant length across all ALT alleles is computed as
   `max_var_len = max(|alt_length|)` (absolute value, since deletions are
   stored as negative lengths).
5. `ASMD = mean(ALT NM) − mean(REF NM) − max_var_len`, computed per-sample.

**Why variant length is subtracted**: ALT reads aligned back to the REF
haplotype carry the variant's own edit distance — a clean 50 bp deletion
contributes +50 NM against REF, but that is the variant itself, not noise.
Subtracting `max_var_len` removes this expected structural difference so ASMD
isolates only *excess noise* from chimeric or paralogous mismapping. For SNVs
(`max_var_len = 1`), the adjustment is minimal.

**Value range**: (−∞, +∞), typically [−5, 20], or `.` (untestable — one or both groups empty)

**Coverage stability**: Mean edit distance converges quickly. ASMD is stable
above 10× per allele. At extreme coverages (1000×+), the mean becomes very
precise, but the expected value remains the same.

**Interpretation**:

| ASMD Range | Meaning |
|:-----------|:--------|
| `.` | Untestable: one or both allele groups empty (no reads to compare). |
| ≈ 0 | ALT and REF reads have similar edit distance — true variant expected |
| 1–3 | Mild excess noise — possibly repetitive region or low-quality library |
| > 5 | Strong signal — ALT reads carry many extra mismatches, likely misaligned from a paralogous locus or chimeric junction |

**Manual interpretation tip**: ASMD should be interpreted together with MQCD.
Paralogous mismapping usually shows both elevated ASMD (excess mismatches in
mismapped ALT reads) and depressed MQCD (low ALT mapping quality). True
variants in repetitive regions may show mild ASMD elevation (1–2) without
MQCD depression.

---

## Site Depth Fold Change (`SDFC`)

**Purpose**: Detect collapsed paralogous mappings and abnormal depth at the
variant site relative to the local background.

**In plain terms**: SDFC answers "is this site getting more reads than its
neighbors?" A value of 1.0 means normal depth; 2.0 means twice the expected
depth — often a sign that reads from a duplicated region are all piling up
at one location.

**Computation**:

1. During `ProcessWindow`, a per-sample window coverage map is built:
   for each sample, `SampleWindowCov = sample_sampled_bases / window_length`.
2. For each variant and each sample, `SDFC = sample_DP / SampleWindowCov`, where
   sample_DP is this sample's total read depth at the variant site (sum across alleles).

**Why per-sample normalization**: Case and control samples have different
sequencing depths. A shared combined coverage would give both samples the same
SDFC value, masking depth spikes in the low-coverage sample and diluting them
in the high-coverage one. Per-sample normalization lets a 200× case and a 30×
control each detect their own paralogous mapping artifacts independently.

**Why window mean coverage**: The window (≥ 1000 bp, enforced minimum) averages
read depth across hundreds of positions, providing a stable local background
estimate that is immune to variant density non-uniformity and single-position
outliers. This is more robust than variant-DP-based approaches
(e.g., EMA of nearby variant depths) which are sparse and outlier-sensitive.

**Value range**: [0, ∞)

**Coverage stability**: Inherently coverage-invariant — SDFC is a ratio (site
depth / window depth). A 2× depth spike produces SDFC ≈ 2.0 at any overall
coverage level.

**Interpretation**:

| SDFC Range | Meaning |
|:-----------|:--------|
| 0.8–1.2 | Normal depth — variant site matches local background |
| > 2.0 | Elevated depth — possible collapsed paralog or segmental duplication mapped to one locus |
| > 5.0 | Extreme — strong paralogous collapse signal |
| < 0.5 | Depleted depth — possible allelic dropout, mapping hole, or deletion |

**Manual interpretation tip**: Use SDFC together with MQCD for a comprehensive
paralog detection strategy: collapsed paralogs produce both elevated SDFC (excess
reads mapped to one site) and depressed MQCD (the paralogous reads have lower
MAPQ). A variant at a site with SDFC > 2.0 and MQCD < −0.3 is very likely a
false positive.

---

## Fragment Start Shannon Entropy (`FSSE`)

**Purpose**: Detect residual PCR duplicates that survive upstream deduplication
tools (Picard MarkDuplicates, samtools markdup).

Traditional deduplication identifies duplicates by exact start position + strand.
Three mechanisms produce biological duplicates that share a PCR origin but
differ by 1–3 bp in mapped start position, evading detection:

1. **Exonuclease fraying** — proofreading 3'→5' exonucleases nibble 1–3 bp
   from fragment ends during library prep, shifting mapped starts.
2. **Alignment jitter** — soft-clip thresholds and indel penalties cause
   identical molecules to align to slightly different start coordinates.
3. **Representative read roulette** — when MarkDuplicates selects a
   representative from a duplicate group, different runs may pick different
   members, leaving residual PCR siblings.

FSSE quantifies the spatial diversity of ALT read start positions. True
variants produce reads from independent sampling events scattered across
many start positions (high entropy). PCR duplicates cluster at few start
positions (low entropy).

**Computation**:

1. Collect alignment start positions for all ALT-supporting reads.
2. Bin starts into 3 bp buckets (floor(start / 3)) to absorb exonuclease
   fraying. The 3 bp bin width corresponds to the typical exonuclease
   nibble range.
3. Compute Shannon entropy over the bin counts:
   `H = -Σ(p_i × log₂(p_i))`, where `p_i = count_in_bin_i / total_ALT_reads`.
4. Normalize: `FSSE = H / log₂(min(N, 20))`, where N = number of ALT reads.
   The cap at 20 prevents normalization from penalizing deeply sequenced
   sites with hundreds of unique start positions.
5. If fewer than 3 ALT reads, return `.` (untestable — entropy requires
   a meaningful distribution).

**Value range**: [0, 1], or `.` (untestable)

**Interpretation**:

| FSSE Range | Meaning |
|:-----------|:--------|
| > 0.7 | High start diversity — reads from independent sampling events |
| 0.3–0.7 | Moderate diversity — inspect manually |
| < 0.3 | Low diversity — reads cluster at few start positions, likely residual PCR duplicates |

**Coverage stability**: Perfectly coverage-invariant — the normalization by
log₂(min(N, 20)) ensures FSSE measures start position *diversity* independent
of total depth. A 50% duplicate contamination rate produces the same FSSE
value at 30× and 300×.

---

## ALT-Haplotype Discordance Delta (`AHDD`)

**Purpose**: Detect assembly hallucinations — variants where the graph
assembler constructed an ALT haplotype that does not actually match the
reads supporting it.

**Computation**:

1. During genotyping, each read is re-aligned (minimap2) against the SPOA
   consensus haplotype it was assigned to. The edit distance (NM) against
   its *own* haplotype is stored as `mOwnHapNm`.
2. Reads are grouped by allele assignment (REF vs ALT).
3. `AHDD = mean(ALT reads' NM vs ALT haplotype) − mean(REF reads' NM vs REF haplotype)`.
4. If either group is empty, return `.` (untestable).

**Rationale**: For a true variant, ALT reads should align well to the ALT
haplotype (low NM) and REF reads should align well to the REF haplotype
(low NM), producing AHDD ≈ 0. When the assembler hallucinates a spurious
path, the ALT reads do not actually match their assigned haplotype
(high NM), producing AHDD >> 0.

**Value range**: (−∞, +∞), or `.` (untestable)

**Interpretation**:

| AHDD Range | Meaning |
|:-----------|:--------|
| < 0.5 | Good haplotype-read concordance — reads match their assigned assembly path |
| 0.5–2.0 | Moderate discordance — possible noisy region or complex variant |
| > 2.0 | High discordance — ALT reads poorly match the assembled ALT haplotype, likely assembly hallucination |

**Coverage stability**: Near-invariant. The mean NM converges rapidly
(stabilizes above ~10 reads per allele). Minor variation at very low per-allele
coverage (< 5 reads) due to sampling noise.

---

## Haplotype Segregation Entropy (`HSE`)

**Purpose**: Measure how concentrated ALT reads are on a single SPOA
haplotype path vs. scattered across many.

**Rationale**: True variants arise from a single biological allele. All
ALT reads should consistently align to the same assembled haplotype path,
producing a low-entropy distribution (HSE near 0). Random
sequencing errors scatter across multiple haplotypes
because each error generates a slightly different assembly path, producing
a high-entropy distribution (HSE near 1).

**Computation**:

1. After genotyping assigns reads to alleles and SPOA paths, collect the
   haplotype ID for each ALT-supporting read.
2. Compute Shannon entropy over haplotype assignment counts:
   `H = -Σ(p_i × log₂(p_i))`, where `p_i = reads_on_haplotype_i / total_ALT_reads`.
3. Normalize: `HSE = H / log₂(total_haplotypes)`, where total_haplotypes is the
   total number of SPOA haplotype paths in this graph component (including REF).
4. If fewer than 3 ALT reads or only 1 haplotype, return `.` (untestable).

**Value range**: [0, 1], or `.` (untestable)

**Interpretation**:

| HSE Range | Meaning |
|:----------|:--------|
| < 0.2 | Concentrated — ALT reads consistently map to one haplotype (strong true variant signal) |
| 0.2–0.5 | Moderate concentration — possible complex variant or noisy assembly |
| > 0.5 | Scattered — ALT reads distribute across multiple haplotypes (likely noise or error) |

**Coverage stability**: Perfectly coverage-invariant — HSE measures
concentration shape (entropy), not magnitude. The same 80/20 read
split across two haplotypes produces the same HSE at any depth.

---

## Path Depth Coefficient of Variation (`PDCV`)

**Purpose**: Detect chimeric assembly artifacts from the de Bruijn graph
construction. Chimeric junctions (where unrelated sequences are stitched
together) produce paths with sharp, localized drops in k-mer coverage.

**Computation**:

1. During graph construction, each node in the de Bruijn graph tracks its
   k-mer coverage (the number of reads contributing to that k-mer).
2. For each assembled haplotype path, the coverage of each node along the
   path is collected, and the coefficient of variation is computed:
   `CV = σ(node_coverages) / μ(node_coverages)`.
3. PDCV takes the maximum CV across all ALT haplotype paths in the graph
   component. The worst-case path highlights the most suspicious junction.
4. If a path has fewer than 2 nodes, CV is undefined (variance requires
   ≥ 2 data points) and the metric returns `.`.

**Value range**: [0, ∞), or `.` (untestable)

**Interpretation**:

| PDCV Range | Meaning |
|:-----------|:--------|
| < 0.3 | Uniform path coverage — consistent k-mer support along the entire ALT path |
| 0.3–0.8 | Moderate variation — possible coverage fluctuation from repetitive sequences |
| > 0.8 | High variation — sharp coverage drops along the path, likely chimeric junction or assembly artifact |

**Coverage stability**: Near-invariant — CV is a ratio (σ/μ). A chimeric
junction producing a 10× drop at one node produces similar PDCV at 30×
and 300× total coverage. Minor noise at very low k-mer coverage (< 3×)
where integer quantization effects become significant.
