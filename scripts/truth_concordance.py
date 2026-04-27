#!/usr/bin/env python3
"""Truth Set Concordance — Deep Analysis for Lancet2.

Answers two questions independently:
1. Sensitivity: For every expected variant, is Lancet2 calling it? If not, why?
2. Specificity: For every Lancet2 call, do reads in the sample(s) support it?

Both questions use a single read evidence engine that mirrors Lancet2's
ReadCollector filters (duplicate, mapq<20, qcfail) and reports per-filter
ALT read counts — showing exactly why a variant appears unsupported.

Output files:
  truth_concordance_report.txt  Rich console report (sections 1-4, all tables)
  concordance_details.txt       Per-variant TSV: truth variants + match level
  missed_variants.txt           Per-variant TSV: missed variants + per-filter evidence
  lancet_unmatched_details.txt  Per-variant TSV: unmatched Lancet calls + evidence

Analysis pipeline:
  S1  Concordance funnel — 5 match levels + MISS:
        L0:   exact POS + REF + ALT
        LD:   truth SNV decomposed into a Lancet2 MNP spanning it
        L1:   exact POS, same type, size within 20% (indels only)
        L2:   +/-5bp POS, same type, size within 50% (indels only)
        L3:   +/-50bp POS, same type, size within 50% (indels only)
        MISS: no match found
  S2  Concordance details TSV
  S3  Sensitivity — classify missed truth variants by read support:
        Tier 1: CIGAR-based ALT allele detection at the variant position
        Tier 2: bwa-mem2 local realignment for all tier1=0 variants
  S4  Specificity — classify unmatched Lancet calls by read support

Usage:
    pixi run -e hts-tools python3 scripts/truth_concordance.py \
        --truth-small data/expected_small_variants_giab.chr1.vcf.gz \
        --truth-large data/expected_large_variants_manta.chr1.vcf.gz \
        --lancet data/post_weights.chr1.tmp.vcf.gz \
        --ref data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz \
        --samples data/NA12878.final.cram \
        --output-dir data/ --workers 64 --mode all
"""

from __future__ import annotations

import argparse
import bisect
import os
import random
import re
import shutil
import subprocess
import sys
import threading
import time

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import edlib
import pysam
from rich.console import Console
from rich.table import Table
from tqdm import tqdm


# ============================================================================
# Data Structures
# ============================================================================


@dataclass(frozen=True, slots=True)
class Variant:
    """Canonical variant representation.

    Attributes follow Lancet2's internal representation:
    - vtype: classified via Sequence Core (raw_variant.cpp:L40-L77)
    - variant_length: signed, from CalculateVariantLength (variant_bubble.cpp:L16-L49)
    - alt_index: 0-based position in multi-allelic ALT list (for tracing to VCF line)
    """

    chrom: str
    pos: int  # 1-based, VCF POS
    ref_seq: Optional[str]  # None for symbolic ALTs
    alt_seq: Optional[str]  # None for symbolic ALTs
    vtype: str  # SNV, INS, DEL, MNP, CPX, REF
    variant_length: int  # signed: positive=INS, negative=DEL, 1=SNV
    source: str  # truth_small, truth_large, lancet
    alt_index: int = 0  # 0-based index within multi-allelic ALT list


@dataclass(frozen=True, slots=True)
class ConcordanceResult:
    """Result of matching a truth variant against the Lancet2 call set.

    Attributes:
        truth: the truth variant being matched
        level: match level (L0, LD, L1, L2, L3, MISS)
        matched: the matching Lancet2 variant (None for MISS)
        edit_distance: Levenshtein distance between ALT haplotypes (None for L0/LD/MISS)
    """

    truth: Variant
    level: str
    matched: Optional[Variant] = None
    edit_distance: Optional[int] = None


# ============================================================================
# Output Helpers
# ============================================================================


class DualConsole:
    """Console wrapper that writes to both a file (clean) and terminal (rich colors).

    All rich renderables (tables, styled text) go to both outputs.
    The file gets clean ASCII; the terminal gets full color.
    """

    def __init__(self, output_path: Path) -> None:
        self._file = open(output_path, "w")
        self._file_console = Console(
            file=self._file, width=120, no_color=True, highlight=False,
        )
        self._term_console = Console(width=120, stderr=True)

    def print(self, *args, **kwargs) -> None:
        self._file_console.print(*args, **kwargs)
        self._term_console.print(*args, **kwargs)


def make_console(output_path: Optional[Path] = None) -> DualConsole | Console:
    """Create a console for report output.

    When output_path is provided, returns a DualConsole that writes
    clean ASCII to the file and colored output to the terminal.
    """
    if output_path is not None:
        return DualConsole(output_path)
    return Console(width=120)


# ============================================================================
# Core Algorithms — Sequence Core
# Exact port of RawVariant::ClassifyVariant (raw_variant.cpp:L40-L77)
# ============================================================================


def classify_variant(ref_seq: str, alt_seq: str) -> str:
    """Classify variant type using Lancet2's Sequence Core algorithm.

    Strips matching 5' prefix and 3' suffix between REF and ALT to isolate
    the biological core, then classifies based on core lengths. This decouples
    classification from VCF multi-allelic padding artifacts.

    Exact port of RawVariant::ClassifyVariant (raw_variant.cpp:L40-L77).
    """
    ref_len = len(ref_seq)
    alt_len = len(alt_seq)

    # 5' prefix match
    start_match = 0
    while (
        start_match < ref_len
        and start_match < alt_len
        and ref_seq[start_match] == alt_seq[start_match]
    ):
        start_match += 1

    if start_match == ref_len and start_match == alt_len:
        return "REF"

    # 3' suffix match (bounded by start_match to prevent overlap)
    end_match = 0
    while (
        end_match < (ref_len - start_match)
        and end_match < (alt_len - start_match)
        and ref_seq[-1 - end_match] == alt_seq[-1 - end_match]
    ):
        end_match += 1

    ref_core = ref_len - start_match - end_match
    alt_core = alt_len - start_match - end_match

    if ref_core == 0 and alt_core > 0:
        return "INS"
    if ref_core > 0 and alt_core == 0:
        return "DEL"
    if ref_core == 0 or alt_core == 0:
        return "REF"
    if ref_core != alt_core:
        return "CPX"
    return "SNV" if ref_core == 1 else "MNP"


# ============================================================================
# Variant Length Calculation
# Exact port of CalculateVariantLength (variant_bubble.cpp:L16-L49)
# ============================================================================


def calculate_variant_length(ref_seq: str, alt_seq: str, vtype: str) -> int:
    """Calculate variant length using Lancet2's exact logic.

    Returns signed length: positive for INS, negative for DEL, 1 for SNV.
    For MNP, returns the Sequence Core span (biological mutation width).

    Exact port of CalculateVariantLength (variant_bubble.cpp:L16-L49).
    """
    if vtype == "SNV":
        return 1

    diff = len(alt_seq) - len(ref_seq)

    # INS/DEL/CPX: net structural variance = string length difference
    if vtype in ("INS", "DEL", "CPX"):
        return diff

    # MNP: biological length is the Sequence Core span
    ref_len = len(ref_seq)
    alt_len = len(alt_seq)

    start_match = 0
    while (
        start_match < ref_len
        and start_match < alt_len
        and ref_seq[start_match] == alt_seq[start_match]
    ):
        start_match += 1

    end_match = 0
    while (
        end_match < (ref_len - start_match)
        and end_match < (alt_len - start_match)
        and ref_seq[-1 - end_match] == alt_seq[-1 - end_match]
    ):
        end_match += 1

    return alt_len - start_match - end_match


# ============================================================================
# Size Tolerance — Clamped Hybrid
# ============================================================================


def size_tolerance(abs_size: int, pct: float) -> int:
    """Clamped hybrid size tolerance: min(100, max(3, pct * abs_size)).

    - Floor of 3bp: prevents tiny indels from requiring exact size match
    - Cap of 100bp: prevents very large variants from matching unrelated events
    """
    return min(100, max(3, int(pct * abs_size)))


# ============================================================================
# VCF Loaders
# ============================================================================


def _is_symbolic(alt: str) -> bool:
    """Check if ALT allele is symbolic (e.g., <INS>, <DEL>)."""
    return alt.startswith("<") and alt.endswith(">")


def _extract_info_field(raw, alt_idx: int, fallback=None):
    """Extract per-allele value from pysam INFO field (may be scalar or tuple)."""
    if raw is None:
        return fallback
    if isinstance(raw, tuple):
        return raw[alt_idx] if alt_idx < len(raw) else fallback
    return raw


def _count_vcf_records(vcf_path: Path) -> int:
    """Fast record count via bcftools. Falls back to 0 (unknown total) on error."""
    try:
        result = subprocess.run(
            f'bcftools view -H "{vcf_path}" | wc -l',
            shell=True, capture_output=True, text=True, timeout=120,
        )
        return int(result.stdout.strip())
    except (subprocess.TimeoutExpired, FileNotFoundError, ValueError):
        return 0


def load_small_truth_variants(vcf_path: Path) -> list[Variant]:
    """Load small variant truth set (e.g., GIAB SNVs + indels).

    Applies Sequence Core classification to each ALT allele independently.
    Filters to PASS-only records.
    """
    variants = []
    total = _count_vcf_records(vcf_path)

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in tqdm(vcf, desc="truth-small", unit="rec", total=total or None):
            if rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue

            for alt_idx, alt in enumerate(rec.alts):
                if _is_symbolic(alt):
                    continue

                ref = rec.ref
                pos = rec.pos

                # Normalize: strip multi-allelic shielding padding
                if len(ref) > 1 and len(alt) > 1:
                    core_ref, core_alt, offset = _sequence_core(ref, alt)
                    if core_ref and core_alt:
                        ref = core_ref
                        alt = core_alt
                        pos = pos + offset

                vtype = classify_variant(ref, alt)
                if vtype == "REF":
                    continue

                vlen = calculate_variant_length(ref, alt, vtype)
                variants.append(
                    Variant(
                        chrom=rec.chrom,
                        pos=pos,
                        ref_seq=ref,
                        alt_seq=alt,
                        vtype=vtype,
                        variant_length=vlen,
                        source="truth_small",
                        alt_index=alt_idx,
                    )
                )

    return variants


def load_large_truth_variants(vcf_path: Path) -> list[Variant]:
    """Load large/structural variant truth set (e.g., Manta PASS INS/DEL).

    Uses INFO/SVTYPE when present, otherwise applies Sequence Core.
    Excludes DUP and BND. Filters to PASS-only.
    """
    variants = []
    total = _count_vcf_records(vcf_path)

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in tqdm(vcf, desc="truth-large", unit="rec", total=total or None):
            if rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue

            svtype = rec.info.get("SVTYPE", None)
            if svtype in ("DUP", "BND"):
                continue

            for alt_idx, alt in enumerate(rec.alts):
                if _is_symbolic(alt):
                    if svtype is None:
                        continue
                    svlen_raw = rec.info.get("SVLEN", None)
                    svlen = _extract_info_field(svlen_raw, alt_idx)
                    if svlen is None:
                        continue
                    variants.append(
                        Variant(
                            chrom=rec.chrom,
                            pos=rec.pos,
                            ref_seq=None,
                            alt_seq=None,
                            vtype=svtype,
                            variant_length=int(svlen),
                            source="truth_large",
                            alt_index=alt_idx,
                        )
                    )
                else:
                    ref = rec.ref
                    pos = rec.pos

                    # Normalize: strip multi-allelic shielding padding
                    if len(ref) > 1 and len(alt) > 1:
                        core_ref, core_alt, offset = _sequence_core(ref, alt)
                        if core_ref and core_alt:
                            ref = core_ref
                            alt = core_alt
                            pos = pos + offset

                    vtype = svtype if svtype else classify_variant(ref, alt)
                    if vtype == "REF":
                        continue
                    vlen = (
                        len(alt) - len(ref)
                        if vtype in ("INS", "DEL")
                        else calculate_variant_length(ref, alt, vtype)
                    )
                    variants.append(
                        Variant(
                            chrom=rec.chrom,
                            pos=pos,
                            ref_seq=ref,
                            alt_seq=alt,
                            vtype=vtype,
                            variant_length=vlen,
                            source="truth_large",
                            alt_index=alt_idx,
                        )
                    )

    return variants


def load_lancet_variants(vcf_path: Path) -> list[Variant]:
    """Load Lancet2 output VCF with Sequence Core normalization.

    Multi-allelic records are decomposed into individual REF/ALT pairs, then
    Sequence Core normalization strips matching 5'/3' bases (multi-allelic
    shielding padding) and adjusts POS. This ensures concordance matching
    compares the biological core, not VCF-padded representations.

    Returns list sorted by (chrom, pos) for binary search.
    """
    variants = []
    total = _count_vcf_records(vcf_path)

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in tqdm(vcf, desc="lancet-vcf", unit="rec", total=total or None):
            vtype_raw = rec.info.get("TYPE", None)
            vlen_raw = rec.info.get("LENGTH", None)

            for alt_idx, alt in enumerate(rec.alts):
                if _is_symbolic(alt):
                    continue

                ref = rec.ref
                pos = rec.pos

                # Normalize: strip multi-allelic shielding padding (matching 5'/3'
                # bases left by VCF parsimony when another ALT blocks trimming).
                # Adjusts pos to the core variant's true genomic position.
                if len(ref) > 1 and len(alt) > 1:
                    core_ref, core_alt, offset = _sequence_core(ref, alt)
                    if core_ref and core_alt:  # Guard: don't reduce to empty
                        ref = core_ref
                        alt = core_alt
                        pos = pos + offset

                # Extract per-allele type
                vtype_val = _extract_info_field(vtype_raw, alt_idx)
                # Re-classify after normalization (INFO/TYPE was pre-normalization)
                vtype = classify_variant(ref, alt)

                # Extract per-allele length
                vlen = calculate_variant_length(ref, alt, vtype)

                if vtype == "REF":
                    continue

                variants.append(
                    Variant(
                        chrom=rec.chrom,
                        pos=pos,
                        ref_seq=ref,
                        alt_seq=alt,
                        vtype=vtype,
                        variant_length=vlen,
                        source="lancet",
                        alt_index=alt_idx,
                    )
                )

    # Sort by position for O(log N) bisect lookup in concordance matching
    variants.sort(key=lambda v: (v.chrom, v.pos))
    return variants


# ============================================================================
# Concordance Matching Engine
# ============================================================================


def build_lancet_index(
    lancet_variants: list[Variant],
) -> dict[str, list[tuple[int, Variant]]]:
    """Build position index for O(log N) lookup via bisect.

    Returns dict[chrom] -> sorted list of (pos, Variant).
    """


    index: dict[str, list[tuple[int, Variant]]] = defaultdict(list)
    for v in lancet_variants:
        index[v.chrom].append((v.pos, v))

    # Already sorted from load_lancet_variants, but ensure
    for chrom in index:
        index[chrom].sort(key=lambda x: x[0])

    return dict(index)


def _search_range(
    chrom_index: list[tuple[int, Variant]], pos: int, margin: int
) -> list[Variant]:
    """Return all Lancet2 variants within [pos - margin, pos + margin]."""


    lo = bisect.bisect_left(chrom_index, (pos - margin,))
    hi = bisect.bisect_right(chrom_index, (pos + margin + 1,))
    return [v for _, v in chrom_index[lo:hi]]


def _sequence_core(ref: str, alt: str) -> tuple[str, str, int]:
    """Strip matching 5'/3' bases from REF/ALT, return (core_ref, core_alt, 5prime_offset).

    Port of Lancet2's ClassifyVariant algorithm (raw_variant.cpp). Used to
    normalize multi-allelic shielded variants at load time, stripping VCF padding
    back to the biological core before concordance matching.
    """
    rlen, alen = len(ref), len(alt)
    start = 0
    while start < rlen and start < alen and ref[start] == alt[start]:
        start += 1
    end = 0
    while end < (rlen - start) and end < (alen - start) and ref[-1 - end] == alt[-1 - end]:
        end += 1
    return ref[start : rlen - end], alt[start : alen - end], start


def match_snv(
    truth: Variant, chrom_index: list[tuple[int, Variant]]
) -> tuple[str, Optional[Variant]]:
    """Match a truth SNV/MNP against Lancet2 index.

    Level 0: exact POS + REF + ALT
    Level D: truth SNV decomposed into a Lancet2 MNP that spans it.
             Only MNPs qualify — CPX variants have different REF/ALT lengths,
             so positional offset mapping is invalid after the indel boundary.
    Miss: no match
    """
    # Level 0: exact match
    candidates = _search_range(chrom_index, truth.pos, 0)
    for c in candidates:
        if c.ref_seq == truth.ref_seq and c.alt_seq == truth.alt_seq:
            return ("L0", c)

    # Level D: truth SNV decomposed into a Lancet2 MNP
    # An MNP has equal-length REF and ALT, so offset i in REF corresponds
    # to offset i in ALT. We check if the truth SNV's ref/alt base appears
    # at the correct offset within the MNP.
    if truth.vtype == "SNV":
        candidates = _search_range(chrom_index, truth.pos, 50)
        for c in candidates:
            if c.vtype != "MNP" or c.ref_seq is None or c.alt_seq is None:
                continue

            mnp_start = c.pos
            mnp_end = c.pos + len(c.ref_seq) - 1
            if not (mnp_start <= truth.pos <= mnp_end):
                continue

            offset = truth.pos - c.pos
            if c.ref_seq[offset] == truth.ref_seq and c.alt_seq[offset] == truth.alt_seq:
                return ("LD", c)

    return ("MISS", None)


def match_indel(
    truth: Variant, chrom_index: list[tuple[int, Variant]]
) -> tuple[str, Optional[Variant]]:
    """Match a truth indel/CPX against Lancet2 index.

    Level 0: exact POS + REF + ALT
    Level 1: exact POS, same type, size within tolerance(0.2)
    Level 2: ±5bp, same type, size within tolerance(0.5)
    Level 3: ±50bp, same type, size within tolerance(0.5)
    Miss: no match found
    """
    abs_size = abs(truth.variant_length)

    # Level 0: exact match
    candidates = _search_range(chrom_index, truth.pos, 0)
    for c in candidates:
        if c.ref_seq == truth.ref_seq and c.alt_seq == truth.alt_seq:
            return ("L0", c)

    # Level 1: exact POS, same type, size within 20% tolerance
    tol_1 = size_tolerance(abs_size, 0.2)
    for c in candidates:
        if c.vtype == truth.vtype and abs(abs(c.variant_length) - abs_size) <= tol_1:
            return ("L1", c)

    # Level 2: ±5bp, same type, size within 50% tolerance
    tol_2 = size_tolerance(abs_size, 0.5)
    candidates = _search_range(chrom_index, truth.pos, 5)
    for c in candidates:
        if c.vtype == truth.vtype and abs(abs(c.variant_length) - abs_size) <= tol_2:
            return ("L2", c)

    # Level 3: ±50bp, same type, size within 50% tolerance
    candidates = _search_range(chrom_index, truth.pos, 50)
    for c in candidates:
        if c.vtype == truth.vtype and abs(abs(c.variant_length) - abs_size) <= tol_2:
            return ("L3", c)

    return ("MISS", None)


def match_variant(
    truth: Variant, lancet_index: dict[str, list[tuple[int, Variant]]]
) -> tuple[str, Optional[Variant]]:
    """Top-level dispatcher: routes to SNV/MNP or indel/CPX matcher.

    SNV/MNP: exact match + MNP decomposition check
    INS/DEL/CPX: tiered positional + size matching. CPX goes here because
    it has an indel component (different REF/ALT core lengths).
    """
    chrom_index = lancet_index.get(truth.chrom, [])
    if not chrom_index:
        return ("MISS", None)

    if truth.vtype in ("SNV", "MNP"):
        return match_snv(truth, chrom_index)
    return match_indel(truth, chrom_index)


def build_alt_haplotype(
    variant: Variant,
    fa: pysam.FastaFile,
    region_start: int,
    region_end: int,
) -> Optional[str]:
    """Build the ALT haplotype over a fixed reference region (0-based coords).

    Substitutes the variant's REF→ALT within the region, keeping flanks intact.
    Returns None if variant has no sequence data or REF doesn't match reference.
    """
    if variant.ref_seq is None or variant.alt_seq is None:
        return None

    ref_seq = fa.fetch(variant.chrom, region_start, region_end)
    var_offset = (variant.pos - 1) - region_start
    ref_len = len(variant.ref_seq)

    actual_ref = ref_seq[var_offset : var_offset + ref_len]
    if actual_ref.upper() != variant.ref_seq.upper():
        return None

    return ref_seq[:var_offset] + variant.alt_seq + ref_seq[var_offset + ref_len :]


def compute_haplotype_edit_distance(
    truth: Variant,
    matched: Variant,
    fa: pysam.FastaFile,
) -> Optional[int]:
    """Compute Levenshtein edit distance between ALT haplotypes of two variants.

    Builds both haplotypes over a common reference region and computes edit
    distance using edlib. Returns None if haplotype construction fails.
    Returns 0 if haplotypes are identical (same biological event).

    Flank scales with variant size: max(100, max_variant_size). For small
    variants, 100bp provides ample anchoring. For large variants (e.g., 300bp
    DEL), using 300bp flanks ensures both haplotypes have ≥50% shared context.
    """
    if truth.ref_seq is None or matched.ref_seq is None:
        return None
    if truth.alt_seq is None or matched.alt_seq is None:
        return None

    max_var_size = max(abs(truth.variant_length), abs(matched.variant_length))
    flank = max(100, max_var_size)

    t_start = truth.pos - 1
    t_end = t_start + len(truth.ref_seq)
    m_start = matched.pos - 1
    m_end = m_start + len(matched.ref_seq)

    chrom_len = fa.get_reference_length(truth.chrom)
    region_start = max(0, min(t_start, m_start) - flank)
    region_end = min(chrom_len, max(t_end, m_end) + flank)

    truth_haplo = build_alt_haplotype(truth, fa, region_start, region_end)
    matched_haplo = build_alt_haplotype(matched, fa, region_start, region_end)

    if truth_haplo is None or matched_haplo is None:
        return None
    if truth_haplo == matched_haplo:
        return 0

    result = edlib.align(truth_haplo, matched_haplo, task="distance")
    return result["editDistance"]


def run_concordance(
    truth_variants: list[Variant],
    lancet_index: dict[str, list[tuple[int, Variant]]],
    ref: Optional[Path] = None,
) -> list[ConcordanceResult]:
    """Match all truth variants against the Lancet2 index.

    For L1-L3 matches, computes haplotype edit distance when ref is provided.
    Edit distance quantifies how similar the two variants' ALT haplotypes are,
    independent of VCF representation differences.
    """
    fa = pysam.FastaFile(str(ref)) if ref is not None else None
    results: list[ConcordanceResult] = []

    for tv in tqdm(truth_variants, desc="Matching variants", unit="var"):
        level, matched = match_variant(tv, lancet_index)

        edit_dist: Optional[int] = None
        if fa is not None and matched is not None and level in ("L1", "L2", "L3"):
            edit_dist = compute_haplotype_edit_distance(tv, matched, fa)

        results.append(ConcordanceResult(
            truth=tv, level=level, matched=matched, edit_distance=edit_dist,
        ))

    if fa is not None:
        fa.close()
    return results


# ============================================================================
# Direct Read Evidence Engine
# ============================================================================


@dataclass
class ReadEvidence:
    """Per-filter breakdown of ALT read support at a variant position.

    Mirrors the read filters in Lancet2's ReadCollector (read_collector.cpp).
    Answers: how many reads carry the ALT allele, and for those filtered away,
    which specific filter removed them? This replaces opaque boolean/count-only
    functions with transparent diagnostics for both sensitivity and specificity.
    """

    raw_total_depth: int = 0       # all reads overlapping position (no filter)
    raw_alt_count: int = 0         # ALT reads (no filter)
    passing_total_depth: int = 0   # reads passing all Lancet2 filters
    passing_alt_count: int = 0     # ALT reads passing all filters
    alt_filt_duplicate: int = 0    # ALT reads filtered as duplicate
    alt_filt_mapq: int = 0         # ALT reads filtered by mapq < 20
    alt_filt_qcfail: int = 0      # ALT reads filtered as QC-fail


def _read_has_alt_allele(
    read: pysam.AlignedSegment, variant: Variant, var_pos_0: int,
) -> bool:
    """Check if read carries the ALT allele (no filter applied).

    Pure allele detection: determines whether the read's CIGAR alignment
    shows the variant's ALT sequence at the expected position. Filters
    are handled by the caller (run_batch_read_evidence) to enable per-filter
    accounting.
    """
    if read.is_unmapped or read.query_sequence is None:
        return False

    if variant.vtype == "SNV":
        for qpos, rpos in read.get_aligned_pairs():
            if rpos == var_pos_0 and qpos is not None:
                return read.query_sequence[qpos].upper() == variant.alt_seq[0].upper()
        return False

    if variant.vtype == "INS":
        expected_ins = variant.alt_seq[1:].upper()
        if not expected_ins:
            return False
        aligned_pairs = read.get_aligned_pairs()
        anchor_idx = None
        for idx, (qpos, rpos) in enumerate(aligned_pairs):
            if rpos == var_pos_0 and qpos is not None:
                anchor_idx = idx
                break
        if anchor_idx is None:
            return False
        ins_bases = []
        for qpos, rpos in aligned_pairs[anchor_idx + 1:]:
            if rpos is not None:
                break
            if qpos is not None:
                ins_bases.append(read.query_sequence[qpos].upper())
        return "".join(ins_bases) == expected_ins

    if variant.vtype == "DEL":
        del_len = len(variant.ref_seq) - 1
        if del_len <= 0:
            return False
        del_start = var_pos_0 + 1
        del_end = var_pos_0 + del_len
        deleted_count = 0
        for qpos, rpos in read.get_aligned_pairs():
            if rpos is not None and del_start <= rpos <= del_end:
                if qpos is None:
                    deleted_count += 1
        return deleted_count == del_len

    if variant.vtype in ("MNP", "CPX"):
        for qpos, rpos in read.get_aligned_pairs():
            if rpos == var_pos_0 and qpos is not None:
                return read.query_sequence[qpos].upper() == variant.alt_seq[0].upper()

    return False


def run_batch_read_evidence(
    variants: list[Variant],
    samples: list[Path],
    ref: Path,
    workers: int = 1,
) -> dict[Variant, dict[int, ReadEvidence]]:
    """Compute per-filter read evidence for a batch of variants in parallel.

    Used by BOTH the sensitivity report (§3: are missed truth variants
    supported by reads?) and the specificity report (§4: are unmatched
    Lancet calls supported by reads?). A single engine for both questions
    ensures consistent methodology and provides transparent per-filter
    breakdown everywhere.

    Opens each BAM/CRAM once per worker chunk (not once per variant) to
    amortize the ~100ms CRAM index load. Variants are sorted by (chrom, pos)
    for sequential I/O within each chunk.

    Returns dict[variant] -> dict[sample_idx, ReadEvidence].
    """
    if not variants:
        return {}

    sorted_vars = sorted(variants, key=lambda v: (v.chrom, v.pos))
    n_chunks = max(1, min(workers, len(sorted_vars)))
    chunk_size = (len(sorted_vars) + n_chunks - 1) // n_chunks
    chunks = [sorted_vars[i : i + chunk_size] for i in range(0, len(sorted_vars), chunk_size)]

    counter_lock = threading.Lock()
    counter = [0]

    def _process_chunk(chunk: list[Variant]) -> dict[Variant, dict[int, ReadEvidence]]:
        chunk_results: dict[Variant, dict[int, ReadEvidence]] = {}
        afiles = [pysam.AlignmentFile(str(s), reference_filename=str(ref)) for s in samples]
        try:
            for v in chunk:
                if v.ref_seq is None or v.alt_seq is None:
                    chunk_results[v] = {si: ReadEvidence() for si in range(len(samples))}
                    with counter_lock:
                        counter[0] += 1
                    continue

                var_pos_0 = v.pos - 1
                var_span = max(len(v.ref_seq), len(v.alt_seq))
                fetch_start = max(0, var_pos_0 - 150)
                fetch_end = var_pos_0 + var_span + 150

                per_sample: dict[int, ReadEvidence] = {}
                for si, af in enumerate(afiles):
                    ev = ReadEvidence()
                    for read in af.fetch(v.chrom, fetch_start, fetch_end):
                        ev.raw_total_depth += 1
                        supports_alt = _read_has_alt_allele(read, v, var_pos_0)
                        if supports_alt:
                            ev.raw_alt_count += 1

                        # Apply Lancet2 read filters (read_collector.cpp)
                        if read.is_unmapped or read.is_qcfail:
                            if supports_alt:
                                ev.alt_filt_qcfail += 1
                            continue
                        if read.is_duplicate:
                            if supports_alt:
                                ev.alt_filt_duplicate += 1
                            continue
                        if read.mapping_quality < 20:
                            if supports_alt:
                                ev.alt_filt_mapq += 1
                            continue

                        ev.passing_total_depth += 1
                        if supports_alt:
                            ev.passing_alt_count += 1

                    per_sample[si] = ev
                chunk_results[v] = per_sample
                with counter_lock:
                    counter[0] += 1
        finally:
            for af in afiles:
                af.close()
        return chunk_results

    results: dict[Variant, dict[int, ReadEvidence]] = {}
    with ThreadPoolExecutor(max_workers=n_chunks) as pool:
        futures = [pool.submit(_process_chunk, chunk) for chunk in chunks]

        progress = tqdm(total=len(variants), desc="read-evidence-rich", unit="var")
        while True:
            done_futures = [f for f in futures if f.done()]
            with counter_lock:
                new_count = counter[0]
            delta = new_count - progress.n
            if delta > 0:
                progress.update(delta)
            if len(done_futures) == len(futures):
                break
            time.sleep(0.1)
        with counter_lock:
            final_delta = counter[0] - progress.n
        if final_delta > 0:
            progress.update(final_delta)
        progress.close()

        for future in futures:
            results.update(future.result())

    return results


# ============================================================================
# Local Realignment Engine (bwa-mem2)
# ============================================================================


def group_nearby_variants(
    variants: list[Variant],
    max_gap: int = 1000,
) -> list[list[Variant]]:
    """Group variants within max_gap bp into batches sharing one read region.

    Sorts by (chrom, pos) and merges consecutive variants whose positions
    are within max_gap of each other. Each group is realigned as a single
    bwa-mem2 invocation against a multi-contig haplotype FASTA.
    """
    if not variants:
        return []

    sorted_vars = sorted(variants, key=lambda v: (v.chrom, v.pos))
    groups: list[list[Variant]] = [[sorted_vars[0]]]

    for v in sorted_vars[1:]:
        prev = groups[-1][-1]
        if v.chrom == prev.chrom and (v.pos - prev.pos) <= max_gap:
            groups[-1].append(v)
        else:
            groups.append([v])

    return groups


def build_haplotype_fasta(
    variant: Variant,
    ref_path: Path,
    flank: int = 500,
) -> Optional[tuple[str, str]]:
    """Build REF and ALT haplotype contig sequences for a single variant.

    REF contig: ref[pos-flank : pos+len(ref_seq)+flank]
    ALT contig: ref[pos-flank : pos] + alt_seq + ref[pos+len(ref_seq) : ...]

    Returns (ref_contig, alt_contig) or None for symbolic ALTs.
    """
    if variant.ref_seq is None or variant.alt_seq is None:
        return None

    var_pos_0 = variant.pos - 1  # 0-based

    with pysam.FastaFile(str(ref_path)) as fa:
        chrom_len = fa.get_reference_length(variant.chrom)
        start = max(0, var_pos_0 - flank)
        ref_end = var_pos_0 + len(variant.ref_seq)
        end = min(chrom_len, ref_end + flank)

        left_flank = fa.fetch(variant.chrom, start, var_pos_0)
        ref_span = fa.fetch(variant.chrom, var_pos_0, ref_end)
        right_flank = fa.fetch(variant.chrom, ref_end, end)

    ref_contig = left_flank + ref_span + right_flank
    alt_contig = left_flank + variant.alt_seq + right_flank

    return ref_contig, alt_contig


def run_tier2_realignment(
    variants: list[Variant],
    samples: list[Path],
    ref_path: Path,
    workers: int = 0,
) -> dict[tuple[str, int, int], int]:
    """Run bwa-mem2 local realignment for a list of unsupported variants.

    Groups nearby variants, builds multi-contig REF/ALT FASTA per group,
    extracts reads, realigns with bwa-mem2, and counts ALT-preferring reads.
    Groups are processed in parallel with ThreadPoolExecutor.

    Returns dict keyed by (chrom, pos, sample_idx) -> ALT support count.
    """


    if workers <= 0:
        workers = max(1, (os.cpu_count() or 4) // 2)

    scratch_dir = Path("/scratch/tier2_tmp")
    scratch_dir.mkdir(parents=True, exist_ok=True)

    groups = group_nearby_variants(variants)
    results: dict[tuple[str, int, int], int] = {}

    def process_group(gi: int, group: list[Variant]) -> dict[tuple[str, int, int], int]:
        group_dir = scratch_dir / f"group_{gi}"
        group_dir.mkdir(parents=True, exist_ok=True)
        try:
            return _run_tier2_group(group, gi, group_dir, samples, ref_path)
        finally:
            shutil.rmtree(group_dir, ignore_errors=True)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_group, gi, grp): grp for gi, grp in enumerate(groups)}
        progress = tqdm(total=len(variants), desc="realignment", unit="var")
        for future in as_completed(futures):
            grp = futures[future]
            group_results = future.result()
            results.update(group_results)
            progress.update(len(grp))
        progress.close()

    return results


def _extract_reads_to_fastq(
    sample_path: Path,
    ref_path: Path,
    chrom: str,
    start: int,
    end: int,
    fastq_path: Path,
) -> int:
    """Extract filtered reads from a region and write to FASTQ.

    Applies Lancet2-compatible read filters (qcfail, duplicate, unmapped, MAPQ<20).
    Returns the number of reads written.
    """
    count = 0
    with pysam.AlignmentFile(str(sample_path), reference_filename=str(ref_path)) as af:
        with open(fastq_path, "w") as fq:
            for read in af.fetch(chrom, start, end):
                if read.is_unmapped or read.is_qcfail or read.is_duplicate:
                    continue
                if read.mapping_quality < 20 or read.query_sequence is None:
                    continue

                quals = read.query_qualities
                qual_str = (
                    "".join(chr(q + 33) for q in quals)
                    if quals is not None
                    else "I" * len(read.query_sequence)
                )
                fq.write(f"@{read.query_name}_{count}\n{read.query_sequence}\n+\n{qual_str}\n")
                count += 1
    return count


def _align_reads(fasta_path: Path, fastq_path: Path, sam_path: Path) -> None:
    """Align reads against haplotype contigs with bwa-mem2."""


    # bwa-mem2 only supports short flags (-t, not --threads)
    cmd = ["bwa-mem2", "mem", "-t", "1", str(fasta_path), str(fastq_path)]
    with open(sam_path, "w") as sam_out:
        subprocess.run(cmd, stdout=sam_out, stderr=subprocess.DEVNULL, check=True)


def _count_alt_from_sam(sam_path: Path, group: list[Variant], sample_idx: int) -> dict[tuple[str, int, int], int]:
    """Parse SAM and count reads mapping to ALT contigs per variant.

    Contig naming convention: alt_N where N is the variant's index in group.
    Only counts reads with MAPQ > 0 that map to an alt_ contig.
    """
    counts: dict[tuple[str, int, int], int] = {}
    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        for aln in sam:
            if aln.is_unmapped or aln.mapping_quality == 0:
                continue
            ref_name = aln.reference_name
            if ref_name is None or not ref_name.startswith("alt_"):
                continue
            try:
                vi = int(ref_name.split("_", 1)[1])
            except (IndexError, ValueError):
                continue
            if vi < len(group):
                v = group[vi]
                key = (v.chrom, v.pos, sample_idx)
                counts[key] = counts.get(key, 0) + 1
    return counts


def _run_tier2_group(
    group: list[Variant],
    group_idx: int,
    group_dir: Path,
    samples: list[Path],
    ref_path: Path,
) -> dict[tuple[str, int, int], int]:
    """Process a single variant group through bwa-mem2 haplotype realignment.

    Steps: build FASTA → index → per-sample (extract → align → count).
    Returns dict keyed by (chrom, pos, sample_idx) -> ALT support count.
    """


    n_samples = len(samples)

    # Build multi-contig FASTA with REF + ALT per variant
    fasta_path = group_dir / "haplotypes.fa"
    valid_indices: list[int] = []
    with open(fasta_path, "w") as f:
        for vi, v in enumerate(group):
            haplo = build_haplotype_fasta(v, ref_path)
            if haplo is None:
                continue
            ref_contig, alt_contig = haplo
            f.write(f">ref_{vi}\n{ref_contig}\n>alt_{vi}\n{alt_contig}\n")
            valid_indices.append(vi)

    # Initialize all results to 0
    results: dict[tuple[str, int, int], int] = {
        (v.chrom, v.pos, si): 0 for v in group for si in range(n_samples)
    }

    if not valid_indices:
        return results

    # Index the haplotype FASTA
    subprocess.run(["bwa-mem2", "index", str(fasta_path)], capture_output=True, check=True)

    # Compute union fetch region: all variant spans + 500bp flanks
    fetch_chrom = group[0].chrom
    fetch_start = max(0, min(v.pos - 501 for v in group))
    fetch_end = max(
        v.pos + max(len(v.ref_seq or ""), len(v.alt_seq or "")) + 500
        for v in group
    )

    # Per-sample: extract → align → count
    for si, sample_path in enumerate(samples):
        fastq_path = group_dir / f"sample_{si}.fq"
        sam_path = group_dir / f"sample_{si}.sam"

        read_count = _extract_reads_to_fastq(
            sample_path, ref_path, fetch_chrom, fetch_start, fetch_end, fastq_path,
        )
        if read_count == 0:
            continue

        _align_reads(fasta_path, fastq_path, sam_path)
        results.update(_count_alt_from_sam(sam_path, group, si))

    return results


# ============================================================================
# Reports — Data Summary, Concordance, Sensitivity, Specificity, Forensics
# ============================================================================


def write_data_summary(
    console: Console,
    truth_small: list[Variant],
    truth_large: list[Variant],
    lancet_variants: list[Variant],
) -> None:
    """§1: Input data summary — variant counts and type distributions."""
    console.print("[bold]§1 Input Summary[/bold]\n")

    for name, variants in [
        ("Truth-Small", truth_small),
        ("Truth-Large", truth_large),
        ("Lancet2", lancet_variants),
    ]:
        if not variants:
            continue
        type_counts = Counter(v.vtype for v in variants)
        table = Table(title=f"{name}: {len(variants):,} variants")
        table.add_column("Type")
        table.add_column("Count", justify="right")
        table.add_column("Pct", justify="right")
        for vtype in ["SNV", "INS", "DEL", "MNP", "CPX"]:
            count = type_counts.get(vtype, 0)
            pct = f"{count / len(variants) * 100:.1f}%"
            table.add_row(vtype, str(count), pct)
        console.print(table)
        console.print()


def _ed_cell(tier_results: list[ConcordanceResult], level: str) -> str:
    """Format a cell: 'count (min/med/max)' edit distance for L1-L3, just count otherwise."""
    items = [r for r in tier_results if r.level == level]
    count = len(items)
    if count == 0:
        return "-"
    dists = sorted(d for r in items if (d := r.edit_distance) is not None)
    if dists:
        ed_min, ed_med, ed_max = dists[0], dists[len(dists) // 2], dists[-1]
        return f"{count} ({ed_min}/{ed_med}/{ed_max})"
    return str(count)


SIZE_TIERS = [
    ("1bp",       1,    1),
    ("2-5bp",     2,    5),
    ("6-20bp",    6,   20),
    ("21-50bp",  21,   50),
    ("51-200bp", 51,  200),
    ("201bp+",  201, 999999),
]


def _write_concordance_tsv(
    output_dir: Path,
    results: list[ConcordanceResult],
) -> None:
    """Write concordance_details.txt: ALL concordance results in BED-compatible TSV."""
    tsv_path = output_dir / "concordance_details.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\ttruth_ref\ttruth_alt\ttype\tlength\t"
            "level\tedit_dist\tmatched_chrom\tmatched_pos\tmatched_ref\t"
            "matched_alt\tmatched_type\tmatched_length\tsource\n"
        )
        for r in results:
            tv = r.truth
            # BED: 0-based start, 1-based end
            start = tv.pos - 1
            ref_len = len(tv.ref_seq) if tv.ref_seq else 1
            end = tv.pos - 1 + ref_len

            ed = str(r.edit_distance) if r.edit_distance is not None else "."

            if r.matched is not None:
                mv = r.matched
                m_fields = (
                    f"{mv.chrom}\t{mv.pos}\t"
                    f"{mv.ref_seq or '.'}\t{mv.alt_seq or '.'}\t"
                    f"{mv.vtype}\t{mv.variant_length}"
                )
            else:
                m_fields = ".\t.\t.\t.\t.\t."

            fh.write(
                f"{tv.chrom}\t{start}\t{end}\t"
                f"{tv.ref_seq or '.'}\t{tv.alt_seq or '.'}\t"
                f"{tv.vtype}\t{tv.variant_length}\t"
                f"{r.level}\t{ed}\t{m_fields}\t{tv.source}\n"
            )

    print(f"  Wrote {len(results)} rows to {tsv_path}")


def write_concordance_report(
    console: Console,
    output_dir: Path,
    results: list[ConcordanceResult],
) -> None:
    """§2: Concordance overview — match level distribution and per-type breakdowns."""
    console.print("[bold]§2 Concordance Overview[/bold]\n")

    total = len(results)
    if total == 0:
        console.print("No variants to report.")
        return

    # ── Overview: Level × Type matrix ────────────────────────────────────
    levels = ["L0", "LD", "L1", "L2", "L3", "MISS"]
    types = ["SNV", "INS", "DEL", "MNP", "CPX"]

    overview = Table(title=f"Concordance: {total:,} truth variants")
    overview.add_column("Level")
    for vt in types:
        overview.add_column(vt, justify="right")
    overview.add_column("Total", justify="right")

    level_type_counts: dict[str, Counter] = {lv: Counter() for lv in levels}
    for r in results:
        level_type_counts[r.level][r.truth.vtype] += 1

    for lv in levels:
        row = [lv]
        lv_total = 0
        for vt in types:
            c = level_type_counts[lv].get(vt, 0)
            lv_total += c
            row.append(str(c))
        row.append(str(lv_total))
        overview.add_row(*row)

    row = ["TOTAL"]
    for vt in types:
        c = sum(level_type_counts[lv].get(vt, 0) for lv in levels)
        row.append(str(c))
    row.append(str(total))
    overview.add_row(*row, style="bold")

    console.print(overview)
    console.print()

    l0_count = sum(level_type_counts["L0"].values())
    ld_count = sum(level_type_counts["LD"].values())
    matched = l0_count + ld_count + sum(
        level_type_counts[lv].get(vt, 0) for lv in ["L1", "L2", "L3"] for vt in types
    )
    miss_count = sum(level_type_counts["MISS"].values())
    console.print(f"Exact match (L0): {l0_count:,}/{total:,} = {l0_count / total * 100:.1f}%")
    console.print(f"Total matched (L0+LD+L1-L3): {matched:,}/{total:,} = {matched / total * 100:.1f}%")
    console.print(f"Missed: {miss_count:,}/{total:,} = {miss_count / total * 100:.1f}%\n")

    # ── Per-type tables ──────────────────────────────────────────────────
    # SNV table
    snv_results = [r for r in results if r.truth.vtype == "SNV"]
    if snv_results:
        snv_total = len(snv_results)
        snv_levels = Counter(r.level for r in snv_results)
        snv_tbl = Table(title=f"SNV Concordance ({snv_total:,} variants)")
        snv_tbl.add_column("L0", justify="right")
        snv_tbl.add_column("LD", justify="right")
        snv_tbl.add_column("MISS", justify="right")
        snv_tbl.add_column("Miss %", justify="right")
        miss_pct = f"{snv_levels.get('MISS', 0) / snv_total * 100:.2f}%"
        snv_tbl.add_row(
            str(snv_levels.get("L0", 0)),
            str(snv_levels.get("LD", 0)),
            str(snv_levels.get("MISS", 0)),
            miss_pct,
        )
        console.print(snv_tbl)
        console.print()

    # INS, DEL, CPX tables (with size tiers and edit distance)
    for vtype in ["INS", "DEL", "CPX"]:
        vt_results = [r for r in results if r.truth.vtype == vtype]
        if not vt_results:
            continue

        vt_total = len(vt_results)
        tbl = Table(title=f"{vtype} Concordance ({vt_total:,} variants)")
        tbl.add_column("Size", min_width=8)
        tbl.add_column("L0", justify="right")
        tbl.add_column("L1 (min/med/max ed)", justify="right")
        tbl.add_column("L2 (min/med/max ed)", justify="right")
        tbl.add_column("L3 (min/med/max ed)", justify="right")
        tbl.add_column("MISS", justify="right")
        tbl.add_column("Total", justify="right")
        tbl.add_column("Miss %", justify="right")

        for tier_name, lo, hi in SIZE_TIERS:
            tier = [r for r in vt_results if lo <= abs(r.truth.variant_length) <= hi]
            if not tier:
                continue
            tier_total = len(tier)
            tier_miss = sum(1 for r in tier if r.level == "MISS")
            miss_pct = f"{tier_miss / tier_total * 100:.1f}%" if tier_total > 0 else "-"

            tbl.add_row(
                tier_name,
                _ed_cell(tier, "L0"),
                _ed_cell(tier, "L1"),
                _ed_cell(tier, "L2"),
                _ed_cell(tier, "L3"),
                str(tier_miss) if tier_miss > 0 else "-",
                str(tier_total),
                miss_pct,
            )

        # Totals row
        vt_miss = sum(1 for r in vt_results if r.level == "MISS")
        miss_pct = f"{vt_miss / vt_total * 100:.1f}%" if vt_total > 0 else "-"
        tbl.add_row(
            "TOTAL",
            _ed_cell(vt_results, "L0"),
            _ed_cell(vt_results, "L1"),
            _ed_cell(vt_results, "L2"),
            _ed_cell(vt_results, "L3"),
            str(vt_miss) if vt_miss > 0 else "-",
            str(vt_total),
            miss_pct,
            style="bold",
        )

        console.print(tbl)
        console.print()

    # ── Write full concordance details TSV ────────────────────────────────
    _write_concordance_tsv(output_dir, results)


def write_sensitivity_report(
    console: Console,
    output_dir: Path,
    all_concordance: list[ConcordanceResult],
    samples: list[Path],
    ref: Path,
    workers: int = 32,
) -> list[Variant]:
    """S3: Sensitivity — why are truth variants missed?

    For each MISS variant, runs two tiers of read evidence:
      Tier 1: CIGAR-based ALT allele detection with per-filter breakdown
              (raw count, passing count, filtered-by-duplicate/mapq/qcfail)
      Tier 2: bwa-mem2 local realignment for ALL tier1=0 variants

    Classifies into MISSED_WITH_SUPPORT (>=2 ALT reads), MISSED_WEAK_SUPPORT
    (1 ALT read), or NO_READ_SUPPORT (0 ALT reads).

    Writes missed_variants.txt with full per-filter evidence columns.
    Returns MISSED_WITH_SUPPORT variants for probe-tracking analysis.
    """
    console.print("[bold]§3 Sensitivity Loss Analysis[/bold]\n")

    missed = [r.truth for r in all_concordance if r.level == "MISS"]
    if not missed:
        console.print("No missed variants.")
        return []

    console.print(f"Total missed variants: {len(missed):,}\n")

    # Tier 1: Rich read evidence with per-filter breakdown
    console.print("[bold]Tier 1: Direct read evidence (with per-filter breakdown)[/bold]")
    evidence_results = run_batch_read_evidence(missed, samples, ref, workers=workers)

    # classifications: variant_id -> (variant, classification, ReadEvidence, tier2_count)
    classifications: dict[int, tuple[Variant, str, ReadEvidence, int]] = {}
    tier2_candidates: list[Variant] = []

    for v in missed:
        per_sample = evidence_results.get(v, {})
        # Use the max passing_alt_count across samples as the tier1 metric
        best_ev = max(per_sample.values(), key=lambda e: e.passing_alt_count) if per_sample else ReadEvidence()
        passing_alt = best_ev.passing_alt_count

        if passing_alt >= 2:
            classifications[id(v)] = (v, "MISSED_WITH_SUPPORT", best_ev, -1)
        elif passing_alt == 1:
            classifications[id(v)] = (v, "MISSED_WEAK_SUPPORT", best_ev, -1)
        else:
            classifications[id(v)] = (v, "NO_READ_SUPPORT", best_ev, -1)
            # Tier2 on ALL tier1=0 variants (not just >5bp — rescues ~56% of SNVs)
            tier2_candidates.append(v)

    # Tier 2: Local realignment for ALL 0-support variants
    if tier2_candidates:
        console.print(f"\n[bold]Tier 2: Local realignment ({len(tier2_candidates)} variants)[/bold]")
        t2_results = run_tier2_realignment(tier2_candidates, samples, ref, workers=workers)
        for v in tier2_candidates:
            per_sample = {si: t2_results.get((v.chrom, v.pos, si), 0)
                          for si in range(len(samples))}
            max_alt = max(per_sample.values()) if per_sample else 0
            _, prev_class, ev, _ = classifications[id(v)]
            if max_alt >= 2:
                classifications[id(v)] = (v, "MISSED_WITH_SUPPORT", ev, max_alt)
            elif max_alt == 1:
                classifications[id(v)] = (v, "MISSED_WEAK_SUPPORT", ev, max_alt)
            else:
                classifications[id(v)] = (v, "NO_READ_SUPPORT", ev, max_alt)

    # ── Table 1: Classification Summary ──────────────────────────────────
    console.print()
    class_counts: dict[str, Counter] = {
        "MISSED_WITH_SUPPORT": Counter(),
        "MISSED_WEAK_SUPPORT": Counter(),
        "NO_READ_SUPPORT": Counter(),
    }
    for _, (v, cls, _, _) in classifications.items():
        class_counts[cls][v.vtype] += 1

    types = ["SNV", "INS", "DEL", "MNP", "CPX"]
    summary = Table(title=f"Table 1: Classification Summary ({len(missed):,} missed)")
    summary.add_column("Classification")
    for vt in types:
        summary.add_column(vt, justify="right")
    summary.add_column("Total", justify="right")
    summary.add_column("Pct", justify="right")

    for cls in ["MISSED_WITH_SUPPORT", "MISSED_WEAK_SUPPORT", "NO_READ_SUPPORT"]:
        row = [cls]
        cls_total = sum(class_counts[cls].values())
        for vt in types:
            row.append(str(class_counts[cls].get(vt, 0)))
        row.append(str(cls_total))
        row.append(f"{cls_total / len(missed) * 100:.1f}%")
        summary.add_row(*row)

    console.print(summary)
    console.print()

    # ── Table 2: Read Support Distribution (MISSED_WITH_SUPPORT) ─────────
    mws_all = [(v, ev, t2) for _, (v, cls, ev, t2) in classifications.items()
               if cls == "MISSED_WITH_SUPPORT"]

    support_bins = [
        ("2-3 reads",   2,   3),
        ("4-5 reads",   4,   5),
        ("6-10 reads",  6,  10),
        ("11-20 reads", 11, 20),
        ("21-50 reads", 21, 50),
        ("51+ reads",   51, 999999),
    ]

    if mws_all:
        tbl = Table(title=f"Table 2: Read Support Depth ({len(mws_all):,} MISSED_WITH_SUPPORT)")
        tbl.add_column("ALT reads")
        for vt in types:
            tbl.add_column(vt, justify="right")
        tbl.add_column("Total", justify="right")

        for bin_name, lo, hi in support_bins:
            row = [bin_name]
            bin_total = 0
            for vt in types:
                cnt = sum(1 for v, ev, t2 in mws_all
                          if v.vtype == vt and lo <= max(ev.passing_alt_count, t2 if t2 >= 0 else 0) <= hi)
                bin_total += cnt
                row.append(str(cnt))
            row.append(str(bin_total))
            if bin_total > 0:
                tbl.add_row(*row)

        console.print(tbl)
        console.print()

    # ── Table 3: Engine Attribution ───────────────────────────────────────
    engine_counts: dict[str, Counter] = {
        "Direct read evidence": Counter(),
        "Realignment rescued": Counter(),
        "No read evidence": Counter(),
    }
    for _, (v, cls, ev, t2) in classifications.items():
        if cls in ("MISSED_WITH_SUPPORT", "MISSED_WEAK_SUPPORT"):
            if ev.passing_alt_count >= 2:
                engine_counts["Direct read evidence"][v.vtype] += 1
            elif t2 >= 2:
                engine_counts["Realignment rescued"][v.vtype] += 1
        elif cls == "NO_READ_SUPPORT":
            engine_counts["No read evidence"][v.vtype] += 1

    tbl = Table(title="Table 3: Engine Attribution")
    tbl.add_column("Evidence source")
    for vt in types:
        tbl.add_column(vt, justify="right")
    tbl.add_column("Total", justify="right")
    tbl.add_column("Method")

    notes = {
        "Direct read evidence": "Reads at position directly support ALT allele",
        "Realignment rescued": "bwa-mem2 local realignment found evidence",
        "No read evidence": "Neither method found supporting reads",
    }
    for eng in ["Direct read evidence", "Realignment rescued", "No read evidence"]:
        row = [eng]
        eng_total = sum(engine_counts[eng].values())
        for vt in types:
            row.append(str(engine_counts[eng].get(vt, 0)))
        row.append(str(eng_total))
        row.append(notes[eng])
        tbl.add_row(*row)

    console.print(tbl)
    console.print()

    # ── Table 4: Size × Support Cross-tabulation ─────────────────────────
    for vtype in ["SNV", "INS", "DEL"]:
        vt_list = [(v, ev, t2) for v, ev, t2 in mws_all if v.vtype == vtype]
        if not vt_list:
            continue
        tbl = Table(title=f"Table 4: {vtype} — Size × Read Support ({len(vt_list)} variants)")
        tbl.add_column("Size")
        for bin_name, _, _ in support_bins:
            tbl.add_column(bin_name, justify="right")
        tbl.add_column("Total", justify="right")

        for tier_name, slo, shi in SIZE_TIERS:
            row = [tier_name]
            tier_total = 0
            for bin_name, blo, bhi in support_bins:
                cnt = sum(1 for v, ev, t2 in vt_list
                          if slo <= abs(v.variant_length) <= shi
                          and blo <= max(ev.passing_alt_count, t2 if t2 >= 0 else 0) <= bhi)
                tier_total += cnt
                row.append(str(cnt))
            row.append(str(tier_total))
            if tier_total > 0:
                tbl.add_row(*row)

        tbl.add_row("TOTAL", *[
            str(sum(1 for v, ev, t2 in vt_list
                    if blo <= max(ev.passing_alt_count, t2 if t2 >= 0 else 0) <= bhi))
            for _, blo, bhi in support_bins
        ], str(len(vt_list)), style="bold")

        console.print(tbl)
        console.print()

    # ── Table 5: High-Priority Misses (top 20 by read support) ───────────
    if mws_all:
        ranked = sorted(mws_all, key=lambda x: max(x[1].passing_alt_count, x[2] if x[2] >= 0 else 0), reverse=True)
        top_n = ranked[:20]

        tbl = Table(title=f"Table 5: Top {len(top_n)} High-Priority Misses (most read evidence)")
        tbl.add_column("#", justify="right")
        tbl.add_column("Chrom")
        tbl.add_column("Pos", justify="right")
        tbl.add_column("Type")
        tbl.add_column("Size", justify="right")
        tbl.add_column("Raw ALT", justify="right")
        tbl.add_column("Pass ALT", justify="right")
        tbl.add_column("Filt MAPQ", justify="right")
        tbl.add_column("Filt Dup", justify="right")
        tbl.add_column("Realign", justify="right")
        tbl.add_column("REF/ALT")

        for i, (v, ev, t2) in enumerate(top_n, 1):
            t2_str = str(t2) if t2 >= 0 else "."
            ref_alt = f"{(v.ref_seq or '.')[:20]}/{(v.alt_seq or '.')[:20]}"
            tbl.add_row(
                str(i), v.chrom, f"{v.pos:,}", v.vtype,
                str(abs(v.variant_length)),
                str(ev.raw_alt_count), str(ev.passing_alt_count),
                str(ev.alt_filt_mapq), str(ev.alt_filt_duplicate),
                t2_str, ref_alt,
            )

        console.print(tbl)
        console.print()

    # ── Write missed_variants.txt ────────────────────────────────────────
    tsv_path = output_dir / "missed_variants.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tref\talt\ttype\tlength\t"
            "classification\traw_alt_count\traw_total_depth\t"
            "passing_alt_count\tpassing_total_depth\t"
            "alt_filt_duplicate\talt_filt_mapq\talt_filt_qcfail\t"
            "tier2_alt_count\tmax_alt_count\tengine\tsource\n"
        )
        for _, (v, cls, ev, t2) in sorted(
            classifications.items(), key=lambda x: (x[1][0].chrom, x[1][0].pos)
        ):
            start = v.pos - 1
            ref_len = len(v.ref_seq) if v.ref_seq else 1
            end = v.pos - 1 + ref_len
            t2_str = str(t2) if t2 >= 0 else "."
            max_alt = max(ev.passing_alt_count, t2 if t2 >= 0 else 0)
            engine = "direct" if ev.passing_alt_count >= 2 else ("realignment" if t2 >= 2 else "none")
            fh.write(
                f"{v.chrom}\t{start}\t{end}\t"
                f"{v.ref_seq or '.'}\t{v.alt_seq or '.'}\t"
                f"{v.vtype}\t{v.variant_length}\t"
                f"{cls}\t{ev.raw_alt_count}\t{ev.raw_total_depth}\t"
                f"{ev.passing_alt_count}\t{ev.passing_total_depth}\t"
                f"{ev.alt_filt_duplicate}\t{ev.alt_filt_mapq}\t{ev.alt_filt_qcfail}\t"
                f"{t2_str}\t{max_alt}\t{engine}\t{v.source}\n"
            )

    print(f"  Wrote {len(classifications)} rows to {tsv_path}")

    mws_variants = [v for _, (v, cls, _, _) in classifications.items() if cls == "MISSED_WITH_SUPPORT"]
    return mws_variants


def write_specificity_report(
    console: Console,
    output_dir: Path,
    all_concordance: list[ConcordanceResult],
    lancet_variants: list[Variant],
    samples: list[Path],
    ref: Path,
    workers: int = 32,
) -> None:
    """S4: Specificity — are unmatched Lancet calls supported by reads?

    Identifies Lancet calls not matched to any truth variant, then runs
    the same per-filter read evidence engine used in S3 to determine
    whether each call has genuine read support or is a potential artifact.

    Writes lancet_unmatched_details.txt with full per-filter evidence columns.
    """
    console.print("[bold]§4 Specificity Audit[/bold]\n")

    # Build reverse index: Lancet2 variants consumed by concordance matches
    matched_ids = {id(r.matched) for r in all_concordance if r.matched is not None}
    unmatched = [v for v in lancet_variants if id(v) not in matched_ids]

    console.print(f"Total Lancet2 variants: {len(lancet_variants):,}")
    console.print(f"Matched by truth concordance: {len(matched_ids):,}")
    console.print(f"Unmatched (potential artifacts or novel): {len(unmatched):,}\n")

    if not unmatched:
        console.print("No unmatched variants.")
        return

    # Tier 1: Direct read evidence via pysam (per-filter breakdown)
    console.print("[bold]Tier 1: Direct read evidence (with per-filter breakdown)[/bold]")
    evidence_results = run_batch_read_evidence(unmatched, samples, ref, workers=workers)

    classifications: dict[int, tuple[Variant, str, ReadEvidence]] = {}
    tier2_candidates: list[Variant] = []

    for v in unmatched:
        per_sample = evidence_results.get(v, {})
        best_ev = max(per_sample.values(), key=lambda e: e.passing_alt_count) if per_sample else ReadEvidence()
        passing_alt = best_ev.passing_alt_count

        if passing_alt >= 2:
            classifications[id(v)] = (v, "SUPPORTED", best_ev)
        elif passing_alt == 1:
            classifications[id(v)] = (v, "WEAK_SUPPORT", best_ev)
        else:
            classifications[id(v)] = (v, "UNSUPPORTED", best_ev)
            if abs(v.variant_length) > 5:
                tier2_candidates.append(v)

    # Tier 2: Local realignment for unsupported >5bp variants
    if tier2_candidates:
        console.print(f"\n[bold]Tier 2: Local realignment ({len(tier2_candidates)} variants)[/bold]")
        t2_results = run_tier2_realignment(tier2_candidates, samples, ref, workers=workers)
        for v in tier2_candidates:
            per_sample = {si: t2_results.get((v.chrom, v.pos, si), 0)
                          for si in range(len(samples))}
            max_alt = max(per_sample.values()) if per_sample else 0
            _, _, ev = classifications[id(v)]
            if max_alt >= 2:
                classifications[id(v)] = (v, "SUPPORTED", ev)
            elif max_alt == 1:
                classifications[id(v)] = (v, "WEAK_SUPPORT", ev)

    # Summary table
    console.print()
    class_counts: dict[str, Counter] = {
        "SUPPORTED": Counter(),
        "WEAK_SUPPORT": Counter(),
        "UNSUPPORTED": Counter(),
    }
    for _, (v, cls, _) in classifications.items():
        class_counts[cls][v.vtype] += 1

    types = ["SNV", "INS", "DEL", "MNP", "CPX"]
    summary = Table(title=f"Specificity: {len(unmatched):,} unmatched Lancet2 calls")
    summary.add_column("Classification")
    for vt in types:
        summary.add_column(vt, justify="right")
    summary.add_column("Total", justify="right")
    summary.add_column("Pct", justify="right")

    for cls in ["SUPPORTED", "WEAK_SUPPORT", "UNSUPPORTED"]:
        row = [cls]
        cls_total = sum(class_counts[cls].values())
        for vt in types:
            row.append(str(class_counts[cls].get(vt, 0)))
        row.append(str(cls_total))
        row.append(f"{cls_total / len(unmatched) * 100:.1f}%")
        summary.add_row(*row)

    console.print(summary)
    console.print()

    # Size × Type cross-tab for UNSUPPORTED variants
    unsup = [(v, ev) for _, (v, cls, ev) in classifications.items() if cls == "UNSUPPORTED"]
    if unsup:
        tbl = Table(title=f"UNSUPPORTED by Size × Type ({len(unsup)} variants)")
        tbl.add_column("Size")
        for vt in types:
            tbl.add_column(vt, justify="right")
        tbl.add_column("Total", justify="right")

        for tier_name, lo, hi in SIZE_TIERS:
            tier = [(v, ev) for v, ev in unsup if lo <= abs(v.variant_length) <= hi]
            if not tier:
                continue
            row = [tier_name]
            for vt in types:
                row.append(str(sum(1 for v, _ in tier if v.vtype == vt)))
            row.append(str(len(tier)))
            tbl.add_row(*row)

        row = ["TOTAL"]
        for vt in types:
            row.append(str(sum(1 for v, _ in unsup if v.vtype == vt)))
        row.append(str(len(unsup)))
        tbl.add_row(*row, style="bold")
        console.print(tbl)
        console.print()

    # Top 20 potential artifacts (largest UNSUPPORTED variants)
    if unsup:
        ranked = sorted(unsup, key=lambda x: abs(x[0].variant_length), reverse=True)
        top_n = ranked[:20]
        tbl = Table(title=f"Top {len(top_n)} Potential Artifacts (largest UNSUPPORTED)")
        tbl.add_column("#", justify="right")
        tbl.add_column("Chrom")
        tbl.add_column("Pos", justify="right")
        tbl.add_column("Type")
        tbl.add_column("Size", justify="right")
        tbl.add_column("Raw ALT", justify="right")
        tbl.add_column("Pass ALT", justify="right")
        tbl.add_column("Filt MAPQ", justify="right")
        tbl.add_column("REF/ALT")

        for i, (v, ev) in enumerate(top_n, 1):
            ref_alt = f"{(v.ref_seq or '.')[:20]}/{(v.alt_seq or '.')[:20]}"
            tbl.add_row(
                str(i), v.chrom, f"{v.pos:,}", v.vtype,
                str(abs(v.variant_length)),
                str(ev.raw_alt_count), str(ev.passing_alt_count),
                str(ev.alt_filt_mapq), ref_alt,
            )

        console.print(tbl)
        console.print()

    # Write lancet_unmatched_details.txt (with per-filter evidence)
    tsv_path = output_dir / "lancet_unmatched_details.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tref\talt\ttype\tlength\t"
            "classification\traw_alt_count\traw_total_depth\t"
            "passing_alt_count\tpassing_total_depth\t"
            "alt_filt_duplicate\talt_filt_mapq\talt_filt_qcfail\n"
        )
        for _, (v, cls, ev) in sorted(
            classifications.items(), key=lambda x: (x[1][0].chrom, x[1][0].pos)
        ):
            start = v.pos - 1
            ref_len = len(v.ref_seq) if v.ref_seq else 1
            end = v.pos - 1 + ref_len
            fh.write(
                f"{v.chrom}\t{start}\t{end}\t"
                f"{v.ref_seq or '.'}\t{v.alt_seq or '.'}\t"
                f"{v.vtype}\t{v.variant_length}\t"
                f"{cls}\t{ev.raw_alt_count}\t{ev.raw_total_depth}\t"
                f"{ev.passing_alt_count}\t{ev.passing_total_depth}\t"
                f"{ev.alt_filt_duplicate}\t{ev.alt_filt_mapq}\t{ev.alt_filt_qcfail}\n"
            )

    print(f"  Wrote {len(classifications)} rows to {tsv_path}")

# ============================================================================
# CLI
# ============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Truth Set Concordance — Deep Analysis for Lancet2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--truth-small",
        type=Path,
        help="Small variant truth VCF (e.g., GIAB SNVs + indels)",
    )
    parser.add_argument(
        "--truth-large",
        type=Path,
        help="Large variant truth VCF (e.g., Manta PASS INS/DEL)",
    )
    parser.add_argument(
        "--lancet", type=Path, required=True, help="Lancet2 output VCF"
    )
    parser.add_argument(
        "--ref", type=Path, required=True, help="Reference FASTA (required for CRAM)"
    )
    parser.add_argument(
        "--samples",
        type=Path,
        nargs="+",
        required=True,
        help="One or more BAM/CRAM alignment files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data"),
        help="Directory for all output .txt files (default: data/)",
    )
    parser.add_argument(
        "--mode",
        choices=["small", "large", "all"],
        default="all",
        help="Which truth sets to use (default: all)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Parallel workers for read evidence and realignment analysis (default: 16)",
    )

    return parser


# ============================================================================
# Main
# ============================================================================


def _write_executive_scorecard(
    console: Console,
    all_concordance: list[ConcordanceResult],
    lancet_variants: list[Variant],
) -> None:
    """Write the executive scorecard: key metrics at a glance."""
    total_truth = len(all_concordance)
    if total_truth == 0:
        return

    types = ["SNV", "INS", "DEL"]
    matched_ids = {id(r.matched) for r in all_concordance if r.matched is not None}
    total_lancet = len(lancet_variants)
    total_unmatched = sum(1 for v in lancet_variants if id(v) not in matched_ids)

    tbl = Table(title="Executive Scorecard")
    tbl.add_column("Metric", min_width=14)
    for vt in types:
        tbl.add_column(vt, justify="right")
    tbl.add_column("Overall", justify="right")

    # Sensitivity: (total - MISS) / total
    for label, compute in [
        ("Sensitivity", lambda r, vt: (
            sum(1 for x in r if x.truth.vtype == vt and x.level != "MISS"),
            sum(1 for x in r if x.truth.vtype == vt),
        )),
        ("Exact Match", lambda r, vt: (
            sum(1 for x in r if x.truth.vtype == vt and x.level == "L0"),
            sum(1 for x in r if x.truth.vtype == vt),
        )),
    ]:
        row = [label]
        total_num = total_denom = 0
        for vt in types:
            num, denom = compute(all_concordance, vt)
            total_num += num
            total_denom += denom
            pct = f"{num / denom * 100:.1f}%" if denom else "—"
            row.append(pct)
        pct = f"{total_num / total_denom * 100:.1f}%" if total_denom else "—"
        row.append(pct)
        tbl.add_row(*row)

    # Missed count
    row = ["Missed"]
    total_missed = 0
    for vt in types:
        n = sum(1 for x in all_concordance if x.truth.vtype == vt and x.level == "MISS")
        total_missed += n
        row.append(f"{n:,}")
    row.append(f"{total_missed:,}")
    tbl.add_row(*row)

    # Truth count
    row = ["Truth Count"]
    for vt in types:
        n = sum(1 for x in all_concordance if x.truth.vtype == vt)
        row.append(f"{n:,}")
    row.append(f"{total_truth:,}")
    tbl.add_row(*row, style="dim")

    # Lancet unmatched (specificity proxy)
    row = ["Lancet Unmatched"]
    for vt in types:
        n = sum(1 for v in lancet_variants if id(v) not in matched_ids and v.vtype == vt)
        row.append(f"{n:,}")
    row.append(f"{total_unmatched:,}")
    tbl.add_row(*row)

    console.print(tbl)
    console.print()


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Unified report console ───────────────────────────────────────────
    console = make_console(args.output_dir / "truth_concordance_report.txt")
    console.print("[bold]" + "═" * 56 + "[/bold]")
    console.print("[bold] LANCET2 TRUTH CONCORDANCE REPORT[/bold]")
    console.print("[bold]" + "═" * 56 + "[/bold]\n")

    # ── 1. Load VCFs ─────────────────────────────────────────────────────
    truth_small: list[Variant] = []
    truth_large: list[Variant] = []

    if args.mode in ("small", "all") and args.truth_small:
        print(f"Loading small truth variants from {args.truth_small}...")
        truth_small = load_small_truth_variants(args.truth_small)
        print(f"  Loaded {len(truth_small):,} truth-small variants")

    if args.mode in ("large", "all") and args.truth_large:
        print(f"Loading large truth variants from {args.truth_large}...")
        truth_large = load_large_truth_variants(args.truth_large)
        print(f"  Loaded {len(truth_large):,} truth-large variants")

    print(f"Loading Lancet2 variants from {args.lancet}...")
    lancet_variants = load_lancet_variants(args.lancet)
    print(f"  Loaded {len(lancet_variants):,} Lancet2 variants")

    # ── §1 Data Summary ──────────────────────────────────────────────────
    print("\n§1: Data summary...")
    write_data_summary(console, truth_small, truth_large, lancet_variants)

    # ── §2 Concordance ───────────────────────────────────────────────────
    print("\n§2: Concordance matching...")
    lancet_index = build_lancet_index(lancet_variants)

    all_concordance: list[ConcordanceResult] = []

    if truth_small:
        all_concordance.extend(run_concordance(truth_small, lancet_index, args.ref))

    if truth_large:
        all_concordance.extend(run_concordance(truth_large, lancet_index, args.ref))

    # Executive scorecard (before concordance details)
    _write_executive_scorecard(console, all_concordance, lancet_variants)

    write_concordance_report(console, args.output_dir, all_concordance)

    # ── §3 Sensitivity ───────────────────────────────────────────────────
    print("\n§3: Sensitivity loss analysis...")
    missed_with_support = write_sensitivity_report(
        console, args.output_dir, all_concordance, args.samples, args.ref, args.workers,
    )

    # ── §4 Specificity ───────────────────────────────────────────────────
    print("\n§4: Specificity audit...")
    write_specificity_report(
        console, args.output_dir, all_concordance, lancet_variants,
        args.samples, args.ref, args.workers,
    )


    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

