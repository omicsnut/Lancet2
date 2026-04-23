#!/usr/bin/env python3
"""Truth Set Concordance — Deep Analysis for Lancet2.

Answers two questions independently:
1. Sensitivity: For every expected variant, is Lancet2 calling it? If not, why?
2. Specificity: For every Lancet2 call, do ≥2 reads in the relevant sample(s) support it?

Supports germline and somatic truth sets, single or multi-sample BAM/CRAM input.

Output files:
  truth_concordance_report.txt  Unified rich report (§1–§6, all tables)
  concordance_details.txt       Per-variant TSV: all truth variants + match level
  missed_variants.txt           Per-variant TSV: missed variants + read support
  forensics_details.txt         Per-variant TSV: pipeline stage + MSA edit %
  lancet_unmatched_details.txt  Per-variant TSV: unmatched Lancet calls

Usage (full analysis with forensics + MSA):
    pixi run -e hts-tools python3 scripts/truth_concordance.py \\
        --truth-small data/expected_small_variants_giab.chr1.vcf.gz \\
        --truth-large data/expected_large_variants_manta.chr1.vcf.gz \\
        --lancet data/post_weights.chr1.tmp.vcf.gz \\
        --ref data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz \\
        --samples data/NA12878.final.cram \\
        --log data/post_weights.full_chr1.debug_run.log \\
        --graphs data/post_weights.out_graphs \\
        --output-dir data/ \\
        --workers 64 \\
        --mode all

Arguments:
    --truth-small   Small variant truth VCF (GIAB SNVs/indels)
    --truth-large   Large variant truth VCF (e.g., Manta PASS INS/DEL)
    --lancet        Lancet2 output VCF (required)
    --ref           Reference FASTA, required for CRAM input (required)
    --samples       One or more BAM/CRAM alignment files (required)
    --output-dir    Directory for all output files (default: data/)
    --log           Lancet2 debug log file (enables §5 forensics)
    --graphs        Lancet2 output graphs directory (enables §6 MSA analysis)
    --mode          Which truth sets: small, large, or all (default: all)
    --skip-forensics  Skip §5–§6 pipeline forensics
    --workers       Parallel workers for read evidence analysis (default: 16)
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
                vtype = classify_variant(ref, alt)
                if vtype == "REF":
                    continue

                vlen = calculate_variant_length(ref, alt, vtype)
                variants.append(
                    Variant(
                        chrom=rec.chrom,
                        pos=rec.pos,
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
                            pos=rec.pos,
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
    """Load Lancet2 output VCF.

    Reads INFO/TYPE and INFO/LENGTH directly (already Sequence Core classified).
    Per-allele fields are indexed by ALT position for multi-allelic records.
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

                # Extract per-allele type
                vtype_val = _extract_info_field(vtype_raw, alt_idx)
                vtype = str(vtype_val) if vtype_val is not None else classify_variant(ref, alt)

                # Extract per-allele length
                vlen_val = _extract_info_field(vlen_raw, alt_idx)
                vlen = int(vlen_val) if vlen_val is not None else calculate_variant_length(ref, alt, vtype)

                if vtype == "REF":
                    continue

                variants.append(
                    Variant(
                        chrom=rec.chrom,
                        pos=rec.pos,
                        ref_seq=ref,
                        alt_seq=alt,
                        vtype=vtype,
                        variant_length=vlen,
                        source="lancet",
                        alt_index=alt_idx,
                    )
                )

    # Sort for binary search in Phase 2
    variants.sort(key=lambda v: (v.chrom, v.pos))
    return variants


# ============================================================================
# Positional Concordance (Phase 2)
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
    """Run Phase 2 concordance on all truth variants.

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


def read_supports_variant(read: pysam.AlignedSegment, variant: Variant) -> bool:
    """Check if a single read supports the variant ALT allele via CIGAR alignment.

    Filters: qcfail, duplicate, unmapped, MAPQ < 20 (secondary/supplementary kept).
    SNV: check query base at the aligned position.
    INS: find anchor position, check inserted bases after it.
    DEL: check for deletion spanning the expected range.
    """
    # Read-level filters
    if read.is_unmapped or read.is_qcfail or read.is_duplicate:
        return False
    if read.mapping_quality < 20:
        return False

    if variant.ref_seq is None or variant.alt_seq is None:
        return False

    # 0-based position (VCF POS is 1-based)
    var_pos_0 = variant.pos - 1

    if variant.vtype == "SNV":
        # Find query position aligned to var_pos_0
        aligned_pairs = read.get_aligned_pairs()
        for qpos, rpos in aligned_pairs:
            if rpos == var_pos_0 and qpos is not None:
                query_base = read.query_sequence[qpos].upper()
                return query_base == variant.alt_seq[0].upper()
        return False

    if variant.vtype == "INS":
        # Anchor is at var_pos_0. Look for inserted bases (rpos=None) after anchor.
        expected_ins = variant.alt_seq[1:].upper()
        if not expected_ins:
            return False

        aligned_pairs = read.get_aligned_pairs()
        # Find anchor
        anchor_qpos = None
        for idx, (qpos, rpos) in enumerate(aligned_pairs):
            if rpos == var_pos_0 and qpos is not None:
                anchor_qpos = idx
                break

        if anchor_qpos is None:
            return False

        # Collect inserted bases after anchor (rpos=None entries)
        ins_bases = []
        for qpos, rpos in aligned_pairs[anchor_qpos + 1 :]:
            if rpos is not None:
                break  # No longer in insertion
            if qpos is not None:
                ins_bases.append(read.query_sequence[qpos].upper())

        return "".join(ins_bases) == expected_ins

    if variant.vtype == "DEL":
        # Deletion: expect qpos=None for positions var_pos_0+1 .. var_pos_0+del_len
        del_len = len(variant.ref_seq) - 1  # minus anchor
        if del_len <= 0:
            return False

        aligned_pairs = read.get_aligned_pairs()
        del_start = var_pos_0 + 1
        del_end = var_pos_0 + del_len

        deleted_count = 0
        for qpos, rpos in aligned_pairs:
            if rpos is not None and del_start <= rpos <= del_end:
                if qpos is None:
                    deleted_count += 1

        return deleted_count == del_len

    # MNP/CPX: check first base (simplified)
    if variant.vtype in ("MNP", "CPX"):
        aligned_pairs = read.get_aligned_pairs()
        for qpos, rpos in aligned_pairs:
            if rpos == var_pos_0 and qpos is not None:
                query_base = read.query_sequence[qpos].upper()
                return query_base == variant.alt_seq[0].upper()

    return False


def count_alt_support_cigar(
    variant: Variant,
    samples: list[Path],
    ref: Path,
    early_exit: int = 2,
) -> dict[int, int]:
    """Count ALT-supporting reads per sample using CIGAR walker.

    Returns dict[sample_idx] -> count.
    Early exits after `early_exit` confirmed ALT reads per sample.
    """
    if variant.ref_seq is None or variant.alt_seq is None:
        return {si: 0 for si in range(len(samples))}

    var_pos_0 = variant.pos - 1
    var_span = max(len(variant.ref_seq), len(variant.alt_seq))
    # Fetch window: variant span + one read length (150bp) on each side
    fetch_start = max(0, var_pos_0 - 150)
    fetch_end = var_pos_0 + var_span + 150

    counts: dict[int, int] = {}
    for si, sample_path in enumerate(samples):
        count = 0
        with pysam.AlignmentFile(str(sample_path), reference_filename=str(ref)) as af:
            for read in af.fetch(variant.chrom, fetch_start, fetch_end):
                if read_supports_variant(read, variant):
                    count += 1
                    if count >= early_exit:
                        break
        counts[si] = count

    return counts


def run_batch_cigar(
    variants: list[Variant],
    samples: list[Path],
    ref: Path,
    early_exit: int = 1000,
    workers: int = 1,
) -> dict[Variant, dict[int, int]]:
    """Run read evidence engine on a batch of variants using chunked parallelism.

    Opens each BAM/CRAM once per worker thread (not once per variant) to
    amortize the ~100ms CRAM index load. Variants are split into contiguous
    chunks sorted by (chrom, pos) for sequential I/O.

    Progress is tracked via a thread-safe atomic counter polled by the main
    thread — workers never touch tqdm directly.
    Returns dict[variant] -> dict[sample_idx, alt_count].
    """
    if not variants:
        return {}


    # Sort by position for sequential CRAM access within each chunk
    sorted_vars = sorted(variants, key=lambda v: (v.chrom, v.pos))
    n_chunks = max(1, min(workers, len(sorted_vars)))
    chunk_size = (len(sorted_vars) + n_chunks - 1) // n_chunks
    chunks = [sorted_vars[i : i + chunk_size] for i in range(0, len(sorted_vars), chunk_size)]

    # Thread-safe counter: workers increment, main thread reads for tqdm
    counter_lock = threading.Lock()
    counter = [0]  # mutable container for closure

    def _process_chunk(chunk: list[Variant]) -> dict[Variant, dict[int, int]]:
        chunk_results: dict[Variant, dict[int, int]] = {}
        # Open each sample file once for the entire chunk
        afiles = [pysam.AlignmentFile(str(s), reference_filename=str(ref)) for s in samples]
        try:
            for v in chunk:
                if v.ref_seq is None or v.alt_seq is None:
                    chunk_results[v] = {si: 0 for si in range(len(samples))}
                    with counter_lock:
                        counter[0] += 1
                    continue

                var_pos_0 = v.pos - 1
                var_span = max(len(v.ref_seq), len(v.alt_seq))
                fetch_start = max(0, var_pos_0 - 150)
                fetch_end = var_pos_0 + var_span + 150

                counts: dict[int, int] = {}
                for si, af in enumerate(afiles):
                    count = 0
                    for read in af.fetch(v.chrom, fetch_start, fetch_end):
                        if read_supports_variant(read, v):
                            count += 1
                            if count >= early_exit:
                                break
                    counts[si] = count
                chunk_results[v] = counts
                with counter_lock:
                    counter[0] += 1
        finally:
            for af in afiles:
                af.close()
        return chunk_results

    results: dict[Variant, dict[int, int]] = {}
    with ThreadPoolExecutor(max_workers=n_chunks) as pool:
        futures = [pool.submit(_process_chunk, chunk) for chunk in chunks]

        # Main thread drives tqdm by polling the shared counter
        progress = tqdm(total=len(variants), desc="read-evidence", unit="var")
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
        # Final sync
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
    """§3: Sensitivity loss analysis — classify missed truth variants by read support.

    Returns list of MISSED_WITH_SUPPORT variants for downstream forensics.
    """
    console.print("[bold]§3 Sensitivity Loss Analysis[/bold]\n")

    missed = [r.truth for r in all_concordance if r.level == "MISS"]
    if not missed:
        console.print("No missed variants.")
        return []

    console.print(f"Total missed variants: {len(missed):,}\n")

    # Tier 1: Direct read evidence via pysam
    console.print("[bold]Tier 1: Direct read evidence[/bold]")
    cigar_results = run_batch_cigar(missed, samples, ref, workers=workers)

    classifications: dict[int, tuple[Variant, str, int, int]] = {}
    # (variant_id) -> (variant, classification, tier1_count, tier2_count)
    tier2_candidates: list[Variant] = []

    for v in missed:
        counts = cigar_results.get(v, {})
        max_alt = max(counts.values()) if counts else 0
        if max_alt >= 2:
            classifications[id(v)] = (v, "MISSED_WITH_SUPPORT", max_alt, -1)
        elif max_alt == 1:
            classifications[id(v)] = (v, "MISSED_WEAK_SUPPORT", max_alt, -1)
        else:
            classifications[id(v)] = (v, "NO_READ_SUPPORT", 0, -1)
            if abs(v.variant_length) > 5:
                tier2_candidates.append(v)

    # Tier 2: Local realignment for 0-support >5bp variants
    if tier2_candidates:
        console.print(f"\n[bold]Tier 2: Local realignment ({len(tier2_candidates)} variants)[/bold]")
        t2_results = run_tier2_realignment(tier2_candidates, samples, ref, workers=workers)
        for v in tier2_candidates:
            # t2_results is keyed by (chrom, pos, sample_idx) -> alt_count
            per_sample = {si: t2_results.get((v.chrom, v.pos, si), 0)
                          for si in range(len(samples))}
            max_alt = max(per_sample.values()) if per_sample else 0
            _, prev_class, t1_count, _ = classifications[id(v)]
            if max_alt >= 2:
                classifications[id(v)] = (v, "MISSED_WITH_SUPPORT", t1_count, max_alt)
            elif max_alt == 1:
                classifications[id(v)] = (v, "MISSED_WEAK_SUPPORT", t1_count, max_alt)
            else:
                classifications[id(v)] = (v, "NO_READ_SUPPORT", t1_count, max_alt)

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
    mws_all = [(v, t1, t2) for _, (v, cls, t1, t2) in classifications.items()
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
                cnt = sum(1 for v, t1, t2 in mws_all
                          if v.vtype == vt and lo <= max(t1, t2 if t2 >= 0 else 0) <= hi)
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
    for _, (v, cls, t1, t2) in classifications.items():
        if cls in ("MISSED_WITH_SUPPORT", "MISSED_WEAK_SUPPORT"):
            if t1 >= 2:
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
        vt_list = [(v, t1, t2) for v, t1, t2 in mws_all if v.vtype == vtype]
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
                cnt = sum(1 for v, t1, t2 in vt_list
                          if slo <= abs(v.variant_length) <= shi
                          and blo <= max(t1, t2 if t2 >= 0 else 0) <= bhi)
                tier_total += cnt
                row.append(str(cnt))
            row.append(str(tier_total))
            if tier_total > 0:
                tbl.add_row(*row)

        tbl.add_row("TOTAL", *[
            str(sum(1 for v, t1, t2 in vt_list
                    if blo <= max(t1, t2 if t2 >= 0 else 0) <= bhi))
            for _, blo, bhi in support_bins
        ], str(len(vt_list)), style="bold")

        console.print(tbl)
        console.print()

    # ── Table 5: High-Priority Misses (top 20 by read support) ───────────
    if mws_all:
        ranked = sorted(mws_all, key=lambda x: max(x[1], x[2] if x[2] >= 0 else 0), reverse=True)
        top_n = ranked[:20]

        tbl = Table(title=f"Table 5: Top {len(top_n)} High-Priority Misses (most read evidence)")
        tbl.add_column("#", justify="right")
        tbl.add_column("Chrom")
        tbl.add_column("Pos", justify="right")
        tbl.add_column("Type")
        tbl.add_column("Size", justify="right")
        tbl.add_column("Direct", justify="right")
        tbl.add_column("Realign", justify="right")
        tbl.add_column("Max", justify="right")
        tbl.add_column("REF/ALT")

        for i, (v, t1, t2) in enumerate(top_n, 1):
            t2_str = str(t2) if t2 >= 0 else "."
            max_reads = max(t1, t2 if t2 >= 0 else 0)
            ref_alt = f"{(v.ref_seq or '.')[:20]}/{(v.alt_seq or '.')[:20]}"
            tbl.add_row(
                str(i), v.chrom, f"{v.pos:,}", v.vtype,
                str(abs(v.variant_length)), str(t1), t2_str,
                str(max_reads), ref_alt,
            )

        console.print(tbl)
        console.print()

    # ── Write missed_variants.txt ────────────────────────────────────────
    tsv_path = output_dir / "missed_variants.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tref\talt\ttype\tlength\t"
            "classification\ttier1_alt_count\ttier2_alt_count\t"
            "max_alt_count\tengine\tsource\n"
        )
        for _, (v, cls, t1, t2) in sorted(
            classifications.items(), key=lambda x: (x[1][0].chrom, x[1][0].pos)
        ):
            start = v.pos - 1
            ref_len = len(v.ref_seq) if v.ref_seq else 1
            end = v.pos - 1 + ref_len
            t2_str = str(t2) if t2 >= 0 else "."
            max_alt = max(t1, t2 if t2 >= 0 else 0)
            engine = "direct" if t1 >= 2 else ("realignment" if t2 >= 2 else "none")
            fh.write(
                f"{v.chrom}\t{start}\t{end}\t"
                f"{v.ref_seq or '.'}\t{v.alt_seq or '.'}\t"
                f"{v.vtype}\t{v.variant_length}\t"
                f"{cls}\t{t1}\t{t2_str}\t{max_alt}\t{engine}\t{v.source}\n"
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
    """§4: Specificity audit — classify unmatched Lancet2 calls by read support."""
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

    # Tier 1: Direct read evidence via pysam
    console.print("[bold]Tier 1: Direct read evidence[/bold]")
    cigar_results = run_batch_cigar(unmatched, samples, ref, workers=workers)

    classifications: dict[int, tuple[Variant, str, int]] = {}
    tier2_candidates: list[Variant] = []

    for v in unmatched:
        counts = cigar_results.get(v, {})
        max_alt = max(counts.values()) if counts else 0

        if max_alt >= 2:
            classifications[id(v)] = (v, "SUPPORTED", max_alt)
        elif max_alt == 1:
            classifications[id(v)] = (v, "WEAK_SUPPORT", max_alt)
        else:
            classifications[id(v)] = (v, "UNSUPPORTED", 0)
            if abs(v.variant_length) > 5:
                tier2_candidates.append(v)

    # Tier 2: Local realignment for unsupported >5bp variants
    if tier2_candidates:
        console.print(f"\n[bold]Tier 2: Local realignment ({len(tier2_candidates)} variants)[/bold]")
        t2_results = run_tier2_realignment(tier2_candidates, samples, ref, workers=workers)
        for v in tier2_candidates:
            # t2_results is keyed by (chrom, pos, sample_idx) -> alt_count
            per_sample = {si: t2_results.get((v.chrom, v.pos, si), 0)
                          for si in range(len(samples))}
            max_alt = max(per_sample.values()) if per_sample else 0
            if max_alt >= 2:
                classifications[id(v)] = (v, "SUPPORTED", max_alt)
            elif max_alt == 1:
                classifications[id(v)] = (v, "WEAK_SUPPORT", max_alt)

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
    unsup = [(v, alt) for _, (v, cls, alt) in classifications.items() if cls == "UNSUPPORTED"]
    if unsup:
        tbl = Table(title=f"UNSUPPORTED by Size × Type ({len(unsup)} variants)")
        tbl.add_column("Size")
        for vt in types:
            tbl.add_column(vt, justify="right")
        tbl.add_column("Total", justify="right")

        for tier_name, lo, hi in SIZE_TIERS:
            tier = [(v, a) for v, a in unsup if lo <= abs(v.variant_length) <= hi]
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
        tbl.add_column("REF/ALT")

        for i, (v, _) in enumerate(top_n, 1):
            ref_alt = f"{(v.ref_seq or '.')[:20]}/{(v.alt_seq or '.')[:20]}"
            tbl.add_row(
                str(i), v.chrom, f"{v.pos:,}", v.vtype,
                str(abs(v.variant_length)), ref_alt,
            )

        console.print(tbl)
        console.print()

    # Write lancet_unmatched_details.txt
    tsv_path = output_dir / "lancet_unmatched_details.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tref\talt\ttype\tlength\t"
            "classification\tmax_alt_count\n"
        )
        for _, (v, cls, max_alt) in sorted(
            classifications.items(), key=lambda x: (x[1][0].chrom, x[1][0].pos)
        ):
            start = v.pos - 1
            ref_len = len(v.ref_seq) if v.ref_seq else 1
            end = v.pos - 1 + ref_len
            fh.write(
                f"{v.chrom}\t{start}\t{end}\t"
                f"{v.ref_seq or '.'}\t{v.alt_seq or '.'}\t"
                f"{v.vtype}\t{v.variant_length}\t"
                f"{cls}\t{max_alt}\n"
            )

    print(f"  Wrote {len(classifications)} rows to {tsv_path}")


# ── Forensics helpers ────────────────────────────────────────────────────

# Lancet2 pipeline stages, ordered by pipeline depth.
# Detection uses hierarchical approach: log messages → MSA haplotype analysis.
# Naming: S{stage_number}_{ABBREV}. S_UNK is unnumbered (not a real stage).
PIPELINE_STAGES = {
    "S0_NONLY":    "Window skipped: all-N reference",
    "S0_REPEAT":  "Window skipped: repeated max-k-mers in ref",
    "S0_INACTIVE":"Window skipped: no mutation evidence / low coverage",
    "S1_LOWCOV":  "Read collection: coverage below MinAnchorCov",
    "S2_NOSRC":   "No valid source/sink anchor found",
    "S3_CYCLE":   "Graph has cycle, all k-values exhausted",
    "S3_COMPLEX": "Graph too complex, all k-values exhausted",
    "S3_NOASM":   "All k-values exhausted, no haplotypes",
    "S4_NOPATH":  "Variant not in assembled MSA haplotypes",
    "S5_NOMSA":   "SPOA MSA produced no variant bubbles",
    "S6_GENO":    "Variant in haplotype but not called after genotyping",
    "S_UNK":      "Window completed, cause undetermined (no --graphs)",
}

# Ordering for "which window got furthest" tiebreaking.
# Higher index = deeper in pipeline = more informative failure.
_STAGE_ORDER = list(PIPELINE_STAGES.keys())


def deduce_window_geometry(log_path: Path) -> tuple[int, int]:
    """Deduce window_length and step_size from Lancet2 debug log.

    Scans for 'Processing window chr:START-END' lines. Window length
    is end - start from the first line. Step size is the minimum
    difference between consecutive same-chrom window starts (sorted).

    Handles interleaved multi-threaded logs where windows from different
    chromosomes and non-consecutive positions appear out of order.
    """
    # Region format: chr1:S-E or {HLA:name}:S-E or HLA*12:A:S-E
    # Greedy \S+ backtracks to find the rightmost :\d+-\d+ (see reference.cpp:L204)
    pattern = re.compile(r"Processing window (\S+):(\d+)-(\d+)")
    per_chrom: dict[str, set[int]] = defaultdict(set)
    window_length = 0

    # Collect enough windows to reliably compute step_size.
    # 200 unique starts per chrom is more than enough.
    max_per_chrom = 200
    with open(log_path) as fh:
        for line in fh:
            m = pattern.search(line)
            if m is None:
                continue
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            if window_length == 0:
                window_length = end - start
            if len(per_chrom[chrom]) < max_per_chrom:
                per_chrom[chrom].add(start)
            # Stop once any chrom has enough
            if any(len(s) >= max_per_chrom for s in per_chrom.values()):
                break

    if window_length == 0:
        raise ValueError(f"No 'Processing window' lines found in {log_path}")

    # Compute step_size: minimum gap between consecutive sorted starts
    step_size = window_length  # default: non-overlapping
    for starts_set in per_chrom.values():
        sorted_starts = sorted(starts_set)
        if len(sorted_starts) < 2:
            continue
        for i in range(1, len(sorted_starts)):
            diff = sorted_starts[i] - sorted_starts[i - 1]
            if 0 < diff < step_size:
                step_size = diff
        break  # one chrom is sufficient

    return window_length, step_size


def _build_window_index(
    log_regions: dict[str, list[str]],
) -> dict[str, list[tuple[int, int, str]]]:
    """Build per-chrom sorted interval list from indexed log regions.

    Returns dict[chrom] -> sorted list of (start, end, region_key).
    Uses only regions that match the 'chr:start-end' format.
    """
    region_re = re.compile(r"^(\S+):(\d+)-(\d+)$")
    per_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    seen: set[str] = set()

    for region_key in log_regions:
        if region_key in seen:
            continue
        seen.add(region_key)
        m = region_re.match(region_key)
        if m is None:
            continue
        chrom = m.group(1)
        start, end = int(m.group(2)), int(m.group(3))
        per_chrom[chrom].append((start, end, region_key))

    for chrom in per_chrom:
        per_chrom[chrom].sort()

    return per_chrom


def _windows_covering(
    chrom: str,
    pos: int,
    window_index: dict[str, list[tuple[int, int, str]]],
) -> list[str]:
    """Find all actual log windows that cover a given 1-based position.

    Uses bisect on the per-chrom sorted interval list for O(log N) lookup.
    A window (start, end) covers pos when start <= pos <= end.
    """
    intervals = window_index.get(chrom, [])
    if not intervals:
        return []

    # Binary search: find rightmost window whose start <= pos
    hi = bisect.bisect_right(intervals, (pos, float("inf"), "")) - 1
    results = []
    for i in range(max(0, hi - 20), min(len(intervals), hi + 2)):
        start, end, region_key = intervals[i]
        if start <= pos <= end:
            results.append(region_key)

    return results


def _index_log_by_region(log_path: Path) -> dict[str, list[str]]:
    """Index debug log lines by window region string for fast lookup."""
    # Greedy \S+ backtracks to match rightmost :\d+-\d+ — handles chr1:S-E,
    # {HLA:name}:S-E, and HLA*12:A:S-E (see reference.cpp:L204)
    region_pattern = re.compile(r"(\S+:\d+-\d+)")
    index: dict[str, list[str]] = defaultdict(list)
    with open(log_path) as fh:
        for line in fh:
            m = region_pattern.search(line)
            if m:
                index[m.group(1)].append(line.rstrip())
    return index


def _classify_window_lines(lines: list[str]) -> str:
    """Classify failure stage from a single window's log lines.

    Sorts lines by ISO 8601 timestamp (lexicographically sortable), then
    walks in chronological order. The last matched lifecycle marker is the
    deepest pipeline stage reached by this window.

    Uses hierarchical detection matching the pipeline order in
    variant_builder.cpp (L191-L307) and graph.cpp (L74-L207).

    Log message sources (exact strings from C++ source):
      S0: "only N bases in reference" (variant_builder.cpp:L199)
      S0: "reference has repeat" (variant_builder.cpp:L206)
      S0: "no evidence of mutation" (variant_builder.cpp:L214)
      S1: "total sample coverage" + "Skipping" (variant_builder.cpp:L226)
      S5: "source/sink was not found" (graph.cpp:L132)
      S7: "Cycle found" (graph.cpp:L157)
      S7: "Graph too complex" (graph.cpp:L168)
      S8: "Assembled path sequence(s)" (graph.cpp:L785)
      S10: "Building MSA" (variant_builder.cpp:L257)
      S11: "No variants found in graph component" (variant_builder.cpp:L273)
      OK: "Genotyped N variant(s)" (variant_builder.cpp:L305)
    """
    stage = "S_UNK"

    for line in sorted(lines):
        if "only N bases in reference" in line:
            stage = "S0_NONLY"
        elif "reference has repeat" in line:
            stage = "S0_REPEAT"
        elif "no evidence of mutation" in line:
            stage = "S0_INACTIVE"
        elif "total sample coverage" in line and "Skipping" in line:
            stage = "S1_LOWCOV"
        elif "source/sink was not found" in line:
            stage = "S2_NOSRC"
        elif "Cycle found" in line:
            stage = "S3_CYCLE"
        elif "Graph too complex" in line:
            stage = "S3_COMPLEX"
        elif "No variants found" in line:
            stage = "S5_NOMSA"
        elif "SKIPPED_NOASM_HAPLOTYPE" in line:
            stage = "S3_NOASM"
        elif "MISSING_NO_MSA_VARIANTS" in line:
            stage = "S5_NOMSA"
        elif "FOUND_GENOTYPED_VARIANT" in line:
            stage = "S_UNK"
        elif "SKIPPED_INACTIVE_REGION" in line:
            stage = "S0_INACTIVE"

    return stage


def classify_failure_stage(
    log_lines: dict[str, list[str]],
    variant: Variant,
    window_index: dict[str, list[tuple[int, int, str]]],
) -> tuple[str, str, str]:
    """Classify the pipeline failure stage for a missed variant.

    A variant may be covered by multiple overlapping windows. The variant
    is missed only if ALL covering windows failed. We report the stage of
    the window that got furthest in the pipeline (most informative failure).
    """
    covering = _windows_covering(variant.chrom, variant.pos, window_index)

    best_stage = None
    best_order = -1
    best_region = "unknown"

    for region in covering:
        lines = log_lines.get(region, [])
        if not lines:
            continue

        stage = _classify_window_lines(lines)
        order = _STAGE_ORDER.index(stage)
        if order > best_order:
            best_stage = stage
            best_order = order
            best_region = region

    if best_stage is None:
        best_stage = "S_UNK"

    return best_stage, PIPELINE_STAGES[best_stage], best_region


# ── MSA haplotype variant check ──────────────────────────────────────────


def _build_msa_index(poa_dir: Path) -> dict[tuple[str, int, int], list[Path]]:
    """Pre-build spatial index of MSA FASTA files.

    Returns dict mapping (chrom, start, end) to sorted list of FASTA paths.
    Filename format: msa__{chrom}_{start}_{end}__c{component}.fasta
    """
    msa_re = re.compile(r"msa__(.+?)_(\d+)_(\d+)__c(\d+)\.fasta")
    index: dict[tuple[str, int, int], list[Path]] = defaultdict(list)
    for f in poa_dir.iterdir():
        m = msa_re.match(f.name)
        if m:
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            index[(chrom, start, end)].append(f)
    for key in index:
        index[key].sort()
    return dict(index)


def _parse_msa_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Parse MSA FASTA, return list of (name, gap-stripped sequence)."""
    haplotypes: list[tuple[str, str]] = []
    try:
        text = fasta_path.read_text().strip()
    except OSError:
        return haplotypes
    if not text:
        return haplotypes
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            name = lines[i][1:].strip()
            seq_lines: list[str] = []
            i += 1
            while i < len(lines) and not lines[i].startswith(">"):
                seq_lines.append(lines[i].strip())
                i += 1
            haplotypes.append((name, "".join(seq_lines).replace("-", "")))
        else:
            i += 1
    return haplotypes


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def _check_variant_in_msa(
    variant: Variant,
    window_index: dict[str, list[tuple[int, int, str]]],
    msa_index: dict[tuple[str, int, int], list[Path]],
    ref_fa: pysam.FastaFile,
) -> float:
    """Check if a variant's ALT allele appears in any MSA haplotype.

    Builds an ALT pseudo-haplotype (ref flanking + ALT allele, 50bp each
    side) and searches for it within each non-ref MSA haplotype using
    edlib semi-global (HW) alignment.

    Returns the best (minimum) edit distance percentage across all
    haplotypes across all covering windows. Returns 100.0 if no MSA
    files exist for any covering window.
    """
    chrom, pos = variant.chrom, variant.pos
    var_ref = variant.ref_seq or ""
    var_alt = variant.alt_seq or ""

    # Build ALT pseudo-haplotype: 50bp flanking + ALT allele
    flank = 50
    center = pos - 1  # 0-based
    left = ref_fa.fetch(chrom, max(0, center - flank), center)
    right = ref_fa.fetch(chrom, center + len(var_ref), center + len(var_ref) + flank)
    alt_pseudo = left + var_alt + right
    alt_pseudo_rc = _revcomp(alt_pseudo)
    alt_len = len(alt_pseudo)

    best_edit_pct = 100.0

    # Find all covering windows
    covering = _windows_covering(chrom, pos, window_index)
    for region in covering:
        m = re.match(r"^(\S+):(\d+)-(\d+)$", region)
        if not m:
            continue
        wchrom, wstart, wend = m.group(1), int(m.group(2)), int(m.group(3))
        msa_paths = msa_index.get((wchrom, wstart, wend), [])
        for fpath in msa_paths:
            for name, seq in _parse_msa_fasta(fpath):
                if name.startswith("ref"):
                    continue
                aln_fwd = edlib.align(alt_pseudo, seq, mode="HW", task="distance")
                aln_rc = edlib.align(alt_pseudo_rc, seq, mode="HW", task="distance")
                edit_pct = min(aln_fwd["editDistance"], aln_rc["editDistance"]) / alt_len * 100
                if edit_pct < best_edit_pct:
                    best_edit_pct = edit_pct

    return best_edit_pct


def _percentile(data: list[float], p: float) -> float:
    """Compute p-th percentile of sorted data using linear interpolation."""
    s = sorted(data)
    k = (len(s) - 1) * p / 100
    lo, hi = int(k), min(int(k) + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def write_forensics_report(
    console: Console,
    output_dir: Path,
    missed_with_support: list[Variant],
    all_concordance: list[ConcordanceResult],
    log_path: Path,
    graphs_dir: Optional[Path] = None,
    ref_path: Optional[Path] = None,
) -> None:
    """§5–6: Pipeline forensics and MSA haplotype analysis.

    For each variant, maps position to covering window(s), parses all log
    messages for those windows, and classifies the failure stage using
    hierarchical detection of the window lifecycle messages.

    When --graphs is provided:
      1. Resolves S_UNK variants via MSA haplotype check:
         - Variant found in haplotype (edit < 10%) → S6_GENO
         - Variant not found → S4_NOPATH
      2. Runs MSA check on called variants as comparison baseline,
         showing the expected edit distance profile for variants that
         ARE successfully called.
    """
    console.print("[bold]§5 Forensics: Pipeline Stage Attribution[/bold]\n")

    if not missed_with_support:
        console.print("No MISSED_WITH_SUPPORT variants to analyze.")
        return

    window_length, step_size = deduce_window_geometry(log_path)
    console.print(f"Window geometry: length={window_length}, step={step_size}")
    console.print(f"Overlap: {window_length - step_size}bp\n")

    print("  Indexing debug log by region...")
    log_lines = _index_log_by_region(log_path)
    print(f"  Indexed {len(log_lines)} unique regions")

    print("  Building window interval index...")
    window_index = _build_window_index(log_lines)
    total_windows = sum(len(v) for v in window_index.values())
    print(f"  {total_windows} windows across {len(window_index)} chromosomes")

    # Build MSA index if graphs directory is provided
    msa_index: dict[tuple[str, int, int], list[Path]] = {}
    ref_fa: Optional[pysam.FastaFile] = None
    poa_dir: Optional[Path] = None
    if graphs_dir and ref_path:
        poa_dir = graphs_dir / "poa_graph"
        if poa_dir.is_dir():
            print("  Building MSA file index...")
            msa_index = _build_msa_index(poa_dir)
            n_files = sum(len(v) for v in msa_index.values())
            print(f"  Indexed {n_files} MSA files")
            ref_fa = pysam.FastaFile(str(ref_path))
        else:
            print(f"  Warning: MSA directory {poa_dir} not found, skipping MSA checks")

    # ── MSA: check called variants as baseline ───────────────────────────
    # Run MSA check on called variants (L0/LD/L1) as baseline. This
    # establishes the expected edit distance profile for variants that
    # ARE in the MSA, validating the 10% threshold empirically.
    called_msa: list[tuple[Variant, float]] = []  # (variant, edit_pct)
    if msa_index and ref_fa and all_concordance:
        called_variants = [r.truth for r in all_concordance if r.level in ("L0", "LD", "L1")]
        if called_variants:
            print(f"  MSA baseline: checking {len(called_variants)} called variants...")
            for v in tqdm(called_variants, desc="msa-baseline", unit="var"):
                epct = _check_variant_in_msa(v, window_index, msa_index, ref_fa)
                called_msa.append((v, epct))

    # ── Classify missed-with-support variants ────────────────────────────
    # classifications: (variant, stage, desc, region, msa_edit_pct)
    classifications: list[tuple[Variant, str, str, str, float]] = []
    for v in tqdm(missed_with_support, desc="forensics", unit="var"):
        stage, desc, region = classify_failure_stage(log_lines, v, window_index)
        msa_edit_pct = -1.0  # sentinel: not checked

        # Resolve S_UNK via MSA haplotype check when graphs are available
        if stage == "S_UNK" and msa_index and ref_fa:
            msa_edit_pct = _check_variant_in_msa(v, window_index, msa_index, ref_fa)
            if msa_edit_pct < 10.0:
                stage = "S6_GENO"
                desc = PIPELINE_STAGES["S6_GENO"]
            else:
                stage = "S4_NOPATH"
                desc = PIPELINE_STAGES["S4_NOPATH"]

        classifications.append((v, stage, desc, region, msa_edit_pct))

    if ref_fa:
        ref_fa.close()

    # Summary: stage × type distribution
    stage_counts: dict[str, Counter] = defaultdict(Counter)
    for v, stage, desc, _, _ in classifications:
        stage_counts[stage][v.vtype] += 1

    types = ["SNV", "INS", "DEL", "MNP", "CPX"]
    summary = Table(title=f"Forensics: {len(classifications):,} MISSED_WITH_SUPPORT variants")
    summary.add_column("Stage")
    summary.add_column("Description", min_width=30)
    for vt in types:
        summary.add_column(vt, justify="right")
    summary.add_column("Total", justify="right")
    summary.add_column("Pct", justify="right")

    for stage in sorted(stage_counts.keys(), key=lambda s: _STAGE_ORDER.index(s)):
        desc = PIPELINE_STAGES[stage]
        row = [stage, desc]
        stage_total = sum(stage_counts[stage].values())
        for vt in types:
            row.append(str(stage_counts[stage].get(vt, 0)))
        row.append(str(stage_total))
        row.append(f"{stage_total / len(classifications) * 100:.1f}%")
        summary.add_row(*row)

    console.print(summary)
    console.print()

    # ── MSA: Called vs Missed ─────────────────────────────────────────────
    # Compare edit distance profiles between called and missed variants
    # to validate the 10% threshold and show baseline comparison.
    if called_msa:
        import statistics

        console.print("[bold]§6 MSA Haplotype Analysis[/bold]\n")
        console.print("[bold]Called vs Missed Edit Distance[/bold]")
        console.print(
            f"  Called: {len(called_msa)} variants checked  |  "
            f"Missed (MSA-resolved): "
            f"{sum(1 for _, _, _, _, e in classifications if e >= 0)}\n"
        )



        # Table: Type × (Called stats | Missed stats)
        cal_tbl = Table(title="MSA Edit Distance: Called vs Missed (by type)")
        cal_tbl.add_column("Type")
        cal_tbl.add_column("Set")
        cal_tbl.add_column("N", justify="right")
        cal_tbl.add_column("Median", justify="right")
        cal_tbl.add_column("P95", justify="right")
        cal_tbl.add_column("Max", justify="right")
        cal_tbl.add_column("<10%", justify="right")

        missed_msa = [(v, e) for v, _, _, _, e in classifications if e >= 0]

        for vtype in ["SNV", "INS", "DEL"]:
            for label, data in [("Called", called_msa), ("Missed", missed_msa)]:
                edits = [e for v, e in data if v.vtype == vtype]
                if not edits:
                    continue
                n_under10 = sum(1 for e in edits if e < 10.0)
                cal_tbl.add_row(
                    vtype if label == "Called" else "",
                    label,
                    str(len(edits)),
                    f"{statistics.median(edits):.1f}%",
                    f"{_percentile(edits, 95):.1f}%",
                    f"{max(edits):.1f}%",
                    f"{n_under10}/{len(edits)} ({n_under10/len(edits)*100:.1f}%)",
                )

        # Overall
        for label, data in [("Called", called_msa), ("Missed", missed_msa)]:
            edits = [e for _, e in data]
            if edits:
                n_under10 = sum(1 for e in edits if e < 10.0)
                cal_tbl.add_row(
                    "ALL" if label == "Called" else "",
                    label,
                    str(len(edits)),
                    f"{statistics.median(edits):.1f}%",
                    f"{_percentile(edits, 95):.1f}%",
                    f"{max(edits):.1f}%",
                    f"{n_under10}/{len(edits)} ({n_under10/len(edits)*100:.1f}%)",
                    style="bold",
                )

        console.print(cal_tbl)
        console.print()

        # Size-stratified comparison for INS and DEL
        for vtype in ["INS", "DEL"]:
            c_items = [(v, e) for v, e in called_msa if v.vtype == vtype]
            m_items = [(v, e) for v, e in missed_msa if v.vtype == vtype]
            if len(c_items) < 5 and len(m_items) < 5:
                continue

            sz_tbl = Table(title=f"Called vs Missed by Size: {vtype}")
            sz_tbl.add_column("Size Tier")
            sz_tbl.add_column("Set")
            sz_tbl.add_column("N", justify="right")
            sz_tbl.add_column("Median", justify="right")
            sz_tbl.add_column("P95", justify="right")
            sz_tbl.add_column("<10%", justify="right")

            for tier_name, lo, hi in SIZE_TIERS:
                for label, data in [("Called", c_items), ("Missed", m_items)]:
                    edits = [e for v, e in data if lo <= abs(v.variant_length) <= hi]
                    if not edits:
                        continue
                    n_under10 = sum(1 for e in edits if e < 10.0)
                    sz_tbl.add_row(
                        tier_name if label == "Called" else "",
                        label,
                        str(len(edits)),
                        f"{statistics.median(edits):.1f}%",
                        f"{_percentile(edits, 95):.1f}%",
                        f"{n_under10}/{len(edits)} ({n_under10/len(edits)*100:.1f}%)",
                    )

            console.print(sz_tbl)
            console.print()

    # ── MSA Haplotype Analysis ────────────────────────────────────────────
    # For S_UNK variants resolved via MSA check, provide multi-faceted
    # breakdown to identify where assembly vs genotyping failures concentrate
    # and validate that the 10% threshold is appropriate across all sizes.

    msa_checked = [(v, stage, epct) for v, stage, _, _, epct in classifications if epct >= 0]
    if msa_checked:
        import statistics

        console.print("[bold]MSA Haplotype Analysis[/bold]")
        console.print(
            f"  {len(msa_checked)} variants resolved via MSA check "
            f"(S_UNK → S6_GENO or S4_NOPATH)\n"
        )

        # ── Table A: Resolution by Type ──────────────────────────────────
        # Top-level split: how many variants per type are assembly failures
        # (S4_NOPATH = variant not in haplotype) vs genotyping failures
        # (S6_GENO = variant IS in haplotype but not called)?
        res_tbl = Table(title="MSA Resolution by Variant Type")
        res_tbl.add_column("Type")
        res_tbl.add_column("S6_GENO\n(in hap)", justify="right")
        res_tbl.add_column("S4_NOPATH\n(not in hap)", justify="right")
        res_tbl.add_column("Total", justify="right")
        res_tbl.add_column("% In Hap", justify="right")

        for vtype in ["SNV", "INS", "DEL"]:
            vt_items = [(v, s, e) for v, s, e in msa_checked if v.vtype == vtype]
            if not vt_items:
                continue
            n_geno = sum(1 for _, s, _ in vt_items if s == "S6_GENO")
            n_nopath = sum(1 for _, s, _ in vt_items if s == "S4_NOPATH")
            total = len(vt_items)
            pct = n_geno / total * 100 if total else 0
            res_tbl.add_row(vtype, str(n_geno), str(n_nopath), str(total), f"{pct:.1f}%")

        totals = len(msa_checked)
        total_geno = sum(1 for _, s, _ in msa_checked if s == "S6_GENO")
        total_nopath = sum(1 for _, s, _ in msa_checked if s == "S4_NOPATH")
        res_tbl.add_row(
            "ALL", str(total_geno), str(total_nopath), str(totals),
            f"{total_geno / totals * 100:.1f}%", style="bold",
        )
        console.print(res_tbl)
        console.print()

        # ── Table B: Size-stratified resolution per type ─────────────────
        # For each variant type, break down S6_GENO vs S4_NOPATH by size
        # tier. This reveals size-specific bottlenecks (e.g., "genotyper
        # rejects 80% of 1bp INS but only 40% of 21-50bp INS").
        for vtype in ["SNV", "INS", "DEL"]:
            vt_items = [(v, s, e) for v, s, e in msa_checked if v.vtype == vtype]
            if len(vt_items) < 3:
                continue

            size_tbl = Table(title=f"MSA Resolution by Size: {vtype}")
            size_tbl.add_column("Size Tier")
            size_tbl.add_column("S6_GENO", justify="right")
            size_tbl.add_column("S4_NOPATH", justify="right")
            size_tbl.add_column("Total", justify="right")
            size_tbl.add_column("% In Hap", justify="right")
            size_tbl.add_column("Med Edit%\n(S6_GENO)", justify="right")
            size_tbl.add_column("Med Edit%\n(S4_NOPATH)", justify="right")

            for tier_name, lo, hi in SIZE_TIERS:
                tier = [(v, s, e) for v, s, e in vt_items if lo <= abs(v.variant_length) <= hi]
                if not tier:
                    continue
                n_g = sum(1 for _, s, _ in tier if s == "S6_GENO")
                n_p = sum(1 for _, s, _ in tier if s == "S4_NOPATH")
                pct = n_g / len(tier) * 100 if tier else 0
                geno_edits = [e for _, s, e in tier if s == "S6_GENO"]
                nopath_edits = [e for _, s, e in tier if s == "S4_NOPATH"]
                med_g = f"{statistics.median(geno_edits):.1f}" if geno_edits else "—"
                med_p = f"{statistics.median(nopath_edits):.1f}" if nopath_edits else "—"
                size_tbl.add_row(tier_name, str(n_g), str(n_p), str(len(tier)),
                                 f"{pct:.1f}%", med_g, med_p)

            # Totals row
            all_g = [e for _, s, e in vt_items if s == "S6_GENO"]
            all_p = [e for _, s, e in vt_items if s == "S4_NOPATH"]
            n_g_total = len(all_g)
            n_p_total = len(all_p)
            pct_total = n_g_total / len(vt_items) * 100 if vt_items else 0
            med_g_t = f"{statistics.median(all_g):.1f}" if all_g else "—"
            med_p_t = f"{statistics.median(all_p):.1f}" if all_p else "—"
            size_tbl.add_row(
                "ALL", str(n_g_total), str(n_p_total), str(len(vt_items)),
                f"{pct_total:.1f}%", med_g_t, med_p_t, style="bold",
            )
            console.print(size_tbl)
            console.print()

        # ── Table C: Edit Distance Profile ───────────────────────────────
        # Validates the 10% threshold: shows edit distance distribution for
        # S6_GENO (should cluster near 0%) and S4_NOPATH (should be well
        # above 10%) across all types. If any type shows overlap near 10%,
        # the threshold needs adjustment for that type.
        ed_tbl = Table(title="MSA Edit Distance Profile (threshold validation)")
        ed_tbl.add_column("Type × Resolution")
        ed_tbl.add_column("N", justify="right")
        ed_tbl.add_column("Min", justify="right")
        ed_tbl.add_column("Median", justify="right")
        ed_tbl.add_column("P95", justify="right")
        ed_tbl.add_column("Max", justify="right")



        for vtype in ["SNV", "INS", "DEL"]:
            for stage_label, stage_key in [("S6_GENO", "S6_GENO"), ("S4_NOPATH", "S4_NOPATH")]:
                edits = [e for v, s, e in msa_checked if v.vtype == vtype and s == stage_key]
                if not edits:
                    continue
                ed_tbl.add_row(
                    f"{vtype} {stage_label}",
                    str(len(edits)),
                    f"{min(edits):.1f}%",
                    f"{statistics.median(edits):.1f}%",
                    f"{_percentile(edits, 95):.1f}%",
                    f"{max(edits):.1f}%",
                )

        # Overall rows
        for stage_label, stage_key in [("S6_GENO", "S6_GENO"), ("S4_NOPATH", "S4_NOPATH")]:
            edits = [e for _, s, e in msa_checked if s == stage_key]
            if edits:
                ed_tbl.add_row(
                    f"ALL {stage_label}",
                    str(len(edits)),
                    f"{min(edits):.1f}%",
                    f"{statistics.median(edits):.1f}%",
                    f"{_percentile(edits, 95):.1f}%",
                    f"{max(edits):.1f}%",
                    style="bold",
                )

        console.print(ed_tbl)
        console.print()

    # Size × Stage cross-tab per variant type
    # One table per type showing where in the pipeline each size tier fails.
    active_stages = sorted(stage_counts.keys(), key=lambda s: _STAGE_ORDER.index(s))
    for vtype in ["SNV", "INS", "DEL"]:
        vt_class = [(v, s) for v, s, _, _, _ in classifications if v.vtype == vtype]
        if len(vt_class) < 5:
            continue

        # Only show stages that have at least one variant of this type
        vt_stages = [s for s in active_stages if any(st == s for _, st in vt_class)]
        if not vt_stages:
            continue

        tbl = Table(title=f"Failure Stage × Size: {vtype} ({len(vt_class)} variants)")
        tbl.add_column("Size")
        for s in vt_stages:
            tbl.add_column(s, justify="right")
        tbl.add_column("Total", justify="right")

        for tier_name, lo, hi in SIZE_TIERS:
            tier = [(v, s) for v, s in vt_class if lo <= abs(v.variant_length) <= hi]
            if not tier:
                continue
            row = [tier_name]
            for s in vt_stages:
                n = sum(1 for _, st in tier if st == s)
                row.append(str(n) if n else "·")
            row.append(str(len(tier)))
            tbl.add_row(*row)

        row = ["TOTAL"]
        for s in vt_stages:
            row.append(str(sum(1 for _, st in vt_class if st == s)))
        row.append(str(len(vt_class)))
        tbl.add_row(*row, style="bold")
        console.print(tbl)
        console.print()

    # Write forensics_details.txt
    tsv_path = output_dir / "forensics_details.txt"
    with open(tsv_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tref\talt\ttype\tlength\t"
            "failure_stage\tstage_description\twindow_region\t"
            "msa_edit_pct\tsource\n"
        )
        for v, stage, desc, region, msa_epct in sorted(
            classifications, key=lambda x: (x[0].chrom, x[0].pos)
        ):
            start = v.pos - 1
            ref_len = len(v.ref_seq) if v.ref_seq else 1
            end = v.pos - 1 + ref_len
            epct_str = f"{msa_epct:.1f}" if msa_epct >= 0 else "."
            fh.write(
                f"{v.chrom}\t{start}\t{end}\t"
                f"{v.ref_seq or '.'}\t{v.alt_seq or '.'}\t"
                f"{v.vtype}\t{v.variant_length}\t"
                f"{stage}\t{desc}\t{region}\t"
                f"{epct_str}\t{v.source}\n"
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
        "--log", type=Path, help="Lancet2 debug log file (enables §5 forensics)"
    )
    parser.add_argument(
        "--graphs", type=Path, help="Lancet2 output graphs dir (enables §6 MSA analysis)"
    )
    parser.add_argument(
        "--mode",
        choices=["small", "large", "all"],
        default="all",
        help="Which truth sets to use (default: all)",
    )
    parser.add_argument(
        "--skip-forensics",
        action="store_true",
        help="Skip §5–§6 pipeline forensics + MSA analysis",
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

    # ── §5-6 Forensics + MSA ─────────────────────────────────────────────
    if not args.skip_forensics and args.log:
        print("\n§5-6: Forensics + MSA analysis...")
        write_forensics_report(
            console, args.output_dir, missed_with_support, all_concordance, args.log,
            graphs_dir=args.graphs, ref_path=args.ref,
        )
    elif not args.log:
        print("\nSkipping forensics (no --log provided)")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

