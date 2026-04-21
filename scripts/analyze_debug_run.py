#!/usr/bin/env python3
"""Analyze Lancet2 debug run output (--verbose + --out-graphs).

Usage:
    pixi run -e hts-tools python3 scripts/analyze_debug_run.py \\
        --log run.log --graphs out_graphs --vcf output.vcf.gz
    pixi run -e hts-tools python3 scripts/analyze_debug_run.py \\
        --log run.log --graphs out_graphs --vcf output.vcf.gz --view coverage

Prerequisites:
    - pixi hts-tools environment (bcftools, rich)
    - Debug run with --verbose and --out-graphs flags
"""

from __future__ import annotations

import argparse
import os
import re
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from tqdm import tqdm

REPO_ROOT = Path(__file__).resolve().parent.parent
console = Console(record=True)

VIEWS = ("overview", "coverage", "msa", "variants", "suspects", "all")
VIEW_HELP = """\
Views:
  overview     Window counts, k-value distribution, assembly stats
  coverage     Node coverage by type, SAMPLE/SHARED ratios per component
  msa          Haplotype counts, REF anchor lengths, boundary gap analysis
  variants     VCF variant type counts, length distributions (requires --vcf)
  suspects     Deep dive into low SAMPLE/SHARED ratio components (requires --graphs)
  all          All of the above (default)
"""


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Data Model
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


@dataclass(frozen=True, slots=True)
class WindowInfo:
    """A successfully genotyped window with its final k-value and components."""

    region: str
    start: int
    end: int
    k: int
    comps: frozenset[int]


@dataclass(slots=True)
class DotNode:
    """A single node parsed from a DOT graph file."""

    color: str
    coverage: int
    length: int

    @property
    def node_type(self) -> str:
        if self.color == "lightblue":
            return "REF"
        if self.color == "orchid":
            return "SHARED"
        return "SAMPLE"


@dataclass(slots=True)
class ComponentStats:
    """Aggregated coverage stats for one graph component."""

    window: str
    comp: int
    k: int
    ref_covs: list[int] = field(default_factory=list)
    shared_covs: list[int] = field(default_factory=list)
    sample_covs: list[int] = field(default_factory=list)

    @property
    def med_ref(self) -> float:
        return statistics.median(self.ref_covs) if self.ref_covs else 0

    @property
    def med_shared(self) -> float:
        return statistics.median(self.shared_covs) if self.shared_covs else 0

    @property
    def med_sample(self) -> float:
        return statistics.median(self.sample_covs) if self.sample_covs else 0

    @property
    def sample_shared_ratio(self) -> float:
        denom = self.med_shared
        return self.med_sample / denom if denom > 0 else 0.0


@dataclass(frozen=True, slots=True)
class MsaInfo:
    """Summary of one MSA FASTA file."""

    path: str
    n_haplotypes: int
    ref_length: int
    left_gap: int
    right_gap: int
    alt_diffs: tuple[int, ...]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Parsing
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_RE_GENOTYPED = re.compile(
    r"Genotyped (\d+) variant\(s\) for window (chr\w+):(\d+)-(\d+)"
)
_RE_ASSEMBLED = re.compile(
    r"Assembled (\d+)bp path sequence for (chr\w+):(\d+)-(\d+) comp=(\d+) with k=(\d+)"
)
_RE_DOT_NODE = re.compile(
    r'\[shape=\w+\s+fillcolor=(\w+)\s+label="([^"]+)"\]'
)


def parse_log(log_path: Path) -> list[WindowInfo]:
    """Extract successful windows and their final k/component info from log."""
    text = log_path.read_text()

    genotyped = set()
    for m in _RE_GENOTYPED.finditer(text):
        region = f"{m.group(2)}:{m.group(3)}-{m.group(4)}"
        genotyped.add(region)

    window_map: dict[str, dict] = {}
    for m in _RE_ASSEMBLED.finditer(text):
        region = f"{m.group(2)}:{m.group(3)}-{m.group(4)}"
        if region not in genotyped:
            continue
        start, end = int(m.group(3)), int(m.group(4))
        comp, k = int(m.group(5)), int(m.group(6))
        if region not in window_map:
            window_map[region] = {"start": start, "end": end, "k": k, "comps": set()}
        window_map[region]["k"] = max(window_map[region]["k"], k)
        window_map[region]["comps"].add(comp)

    return [
        WindowInfo(region=r, start=d["start"], end=d["end"],
                   k=d["k"], comps=frozenset(d["comps"]))
        for r, d in window_map.items()
    ]


def _parse_dot_job(args: tuple) -> tuple | None:
    """Worker function for parallel DOT parsing (must be top-level for pickling)."""
    window, comp, k, dot_path_str = args
    try:
        with open(dot_path_str) as f:
            text = f.read()
    except FileNotFoundError:
        return None
    ref_covs, shared_covs, sample_covs = [], [], []
    for m in _RE_DOT_NODE.finditer(text):
        color = m.group(1)
        cov_m = re.search(r"coverage=(\d+)", m.group(2))
        cov = int(cov_m.group(1)) if cov_m else 0
        if color == "lightblue":
            ref_covs.append(cov)
        elif color == "orchid":
            shared_covs.append(cov)
        else:
            sample_covs.append(cov)
    return (window, comp, k, ref_covs, shared_covs, sample_covs)


def load_component_stats(
    windows: list[WindowInfo], graphs_dir: Path, n_workers: int = 16
) -> list[ComponentStats]:
    """Load coverage stats for each component from DOT files (parallel)."""
    dbg_dir = graphs_dir / "dbg_graph"
    jobs = []
    for win in windows:
        chrom_start_end = win.region.replace(":", "_").replace("-", "_")
        for comp in win.comps:
            dot_path = dbg_dir / f"dbg__{chrom_start_end}__fully_pruned__k{win.k}__comp{comp}.dot"
            jobs.append((win.region, comp, win.k, str(dot_path)))

    results = []
    with Pool(n_workers) as pool:
        for result in tqdm(pool.imap(_parse_dot_job, jobs, chunksize=256),
                           total=len(jobs), desc="DOT files", unit="file"):
            if result is None:
                continue
            window, comp, k, ref_covs, shared_covs, sample_covs = result
            stats = ComponentStats(window=window, comp=comp, k=k)
            stats.ref_covs = ref_covs
            stats.shared_covs = shared_covs
            stats.sample_covs = sample_covs
            results.append(stats)
    return results


def _parse_msa_job(fpath_str: str) -> MsaInfo | None:
    """Worker function for parallel MSA parsing (must be top-level for pickling)."""
    try:
        with open(fpath_str) as f:
            lines = f.read().strip().split("\n")
    except FileNotFoundError:
        return None
    if len(lines) < 2:
        return None
    ref_msa = lines[1]
    ref_raw = len(ref_msa.replace("-", ""))
    left_gap = len(ref_msa) - len(ref_msa.lstrip("-"))
    right_gap = len(ref_msa) - len(ref_msa.rstrip("-"))
    diffs = []
    for i in range(3, len(lines), 2):
        alt_raw = len(lines[i].replace("-", ""))
        diffs.append(alt_raw - ref_raw)
    fname = fpath_str.rsplit("/", 1)[-1]
    return MsaInfo(path=fname, n_haplotypes=len(lines) // 2,
                   ref_length=ref_raw, left_gap=left_gap, right_gap=right_gap,
                   alt_diffs=tuple(diffs))


def load_msa_info(graphs_dir: Path, n_workers: int = 16) -> list[MsaInfo]:
    """Parse all MSA FASTA files for boundary gaps and haplotype counts (parallel)."""
    poa_dir = graphs_dir / "poa_graph"
    fasta_files = sorted(str(p) for p in poa_dir.glob("msa__*.fasta"))
    results = []
    with Pool(n_workers) as pool:
        for info in tqdm(pool.imap(_parse_msa_job, fasta_files, chunksize=256),
                         total=len(fasta_files), desc="MSA files", unit="file"):
            if info is not None:
                results.append(info)
    return results


def parse_vcf(vcf_path: Path) -> list[dict]:
    """Extract variant records with allele support using bcftools query.

    For multi-allelic sites (AD = ref,alt1,alt2,...), each ALT allele gets
    its own record with ref_ad = AD[0] and alt_ad = AD[i+1].
    """
    result = subprocess.run(
        ["pixi", "run", "-e", "hts-tools", "bcftools", "query", "-f",
         "%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/LENGTH"
         "[\t%AD\t%DP]\n", str(vcf_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        console.print(f"[red]ERROR:[/] bcftools failed: {result.stderr.strip()}")
        return []
    variants = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        # parts: CHROM, POS, REF, ALT(s), TYPE(s), LENGTH(s), AD, DP
        alts = parts[3].split(",")
        types = parts[4].split(",")
        lengths = parts[5].split(",")
        ad_str = parts[6]  # e.g. "9,2" or "9,2,3" for multi-allelic
        dp = int(parts[7]) if parts[7] != "." else 0
        ad_vals = [int(x) for x in ad_str.split(",") if x != "."]
        ref_ad = ad_vals[0] if ad_vals else 0

        for i, (alt, vtype, vlen) in enumerate(zip(alts, types, lengths)):
            alt_ad = ad_vals[i + 1] if (i + 1) < len(ad_vals) else 0
            vaf = alt_ad / dp if dp > 0 else 0.0
            variants.append({
                "chrom": parts[0], "pos": int(parts[1]),
                "ref": parts[2], "alt": alt,
                "type": vtype, "len": int(vlen),
                "ref_ad": ref_ad, "alt_ad": alt_ad, "dp": dp, "vaf": vaf,
            })
    return variants


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Report Views
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def make_bar(fraction: float, width: int = 20, style: str = "white") -> Text:
    """Render a horizontal bar using Unicode block characters."""
    fraction = max(0.0, min(1.0, fraction))
    filled = int(fraction * width)
    bar_str = f"{'█' * filled}".ljust(width, "░")
    return Text(bar_str, style=style)


def render_overview(windows: list[WindowInfo]) -> None:
    """Window counts, k-value distribution, assembly stats."""
    k_dist = Counter(w.k for w in windows)
    comp_counts = [len(w.comps) for w in windows]

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Key", style="bold")
    table.add_column("Value")
    table.add_row("Genotyped windows", f"[bold]{len(windows):,}[/]")
    table.add_row("Total components", f"{sum(comp_counts):,}")
    table.add_row("Comps/window", f"median={statistics.median(comp_counts):.0f}, "
                  f"max={max(comp_counts)}")
    console.print(Panel(table, title="[bold]Assembly Overview[/]", border_style="cyan"))

    kt = Table(title="K-value Distribution", title_style="bold",
               border_style="dim")
    kt.add_column("k", justify="right")
    kt.add_column("Windows", justify="right")
    kt.add_column("Pct", justify="right")
    kt.add_column("Bar", no_wrap=True)
    for k in sorted(k_dist):
        pct = 100 * k_dist[k] / len(windows)
        kt.add_row(str(k), f"{k_dist[k]:,}", f"{pct:.1f}%",
                    make_bar(pct / 100, style="cyan"))
    console.print(kt)


def render_coverage(comp_stats: list[ComponentStats]) -> None:
    """Node coverage by type and SAMPLE/SHARED ratio analysis."""
    all_ref = [c for cs in comp_stats for c in cs.ref_covs]
    all_shared = [c for cs in comp_stats for c in cs.shared_covs]
    all_sample = [c for cs in comp_stats for c in cs.sample_covs]

    ct = Table(title="Node Coverage by Type", title_style="bold",
               border_style="dim")
    ct.add_column("Type")
    ct.add_column("Count", justify="right")
    ct.add_column("Min", justify="right")
    ct.add_column("Median", justify="right")
    ct.add_column("Mean", justify="right")
    ct.add_column("Max", justify="right")
    for name, vals, style in [
        ("REF", all_ref, "blue"), ("SHARED", all_shared, "magenta"),
        ("SAMPLE", all_sample, "green"),
    ]:
        if not vals:
            continue
        ct.add_row(Text(name, style=style), f"{len(vals):,}", str(min(vals)),
                   f"{statistics.median(vals):.0f}", f"{statistics.mean(vals):.0f}",
                   str(max(vals)))
    console.print(ct)

    low_cov = sum(1 for c in all_sample if c < 5)
    console.print(f"  SAMPLE nodes with cov < 5: [yellow]{low_cov:,}[/]/"
                  f"{len(all_sample):,} ({100 * low_cov / len(all_sample):.1f}%)\n")

    # SAMPLE/SHARED ratio distribution
    ratios = [cs.sample_shared_ratio for cs in comp_stats if cs.med_shared > 0]
    if not ratios:
        return

    rt = Table(title="SAMPLE/SHARED Coverage Ratio per Component",
               title_style="bold", border_style="dim",
               caption=f"n={len(ratios)} components | "
                       f"median={statistics.median(ratios):.3f} | "
                       f"mean={statistics.mean(ratios):.3f}")
    rt.add_column("Range")
    rt.add_column("Count", justify="right")
    rt.add_column("Pct", justify="right")
    rt.add_column("Bar", no_wrap=True)
    bins = [(0, 0.1, "<0.10"), (0.1, 0.25, "0.10–0.25"),
            (0.25, 0.5, "0.25–0.50"), (0.5, 0.75, "0.50–0.75"),
            (0.75, 1.0, "0.75–1.00"), (1.0, 2.0, "1.00–2.00"),
            (2.0, 999, ">2.00")]
    for lo, hi, label in bins:
        count = sum(1 for r in ratios if lo <= r < hi)
        if count == 0:
            continue
        pct = 100 * count / len(ratios)
        style = "red" if lo < 0.1 else ("yellow" if lo < 0.25 else "white")
        rt.add_row(label, f"{count:,}", f"{pct:.1f}%",
                   make_bar(pct / 100, style=style))
    console.print(rt)

    # Denominator stability comparison
    both = [cs for cs in comp_stats if cs.ref_covs and cs.shared_covs]
    if len(both) > 10:
        ref_meds = [cs.med_ref for cs in both]
        shared_meds = [cs.med_shared for cs in both]
        all_meds = [statistics.median(cs.ref_covs + cs.shared_covs + cs.sample_covs)
                    for cs in both]
        dt = Table(title="Denominator Stability (CV = stdev/mean, lower = more stable)",
                   title_style="bold", border_style="dim")
        dt.add_column("Denominator")
        dt.add_column("Median", justify="right")
        dt.add_column("Mean", justify="right")
        dt.add_column("CV", justify="right")
        for name, vals in [("REF median", ref_meds),
                           ("SHARED median", shared_meds),
                           ("ALL median", all_meds)]:
            cv = statistics.stdev(vals) / statistics.mean(vals) if statistics.mean(vals) > 0 else 0
            dt.add_row(name, f"{statistics.median(vals):.0f}",
                       f"{statistics.mean(vals):.1f}", f"{cv:.3f}")
        console.print(dt)


def render_msa(msa_info: list[MsaInfo]) -> None:
    """Haplotype counts, REF lengths, and boundary gap analysis."""
    if not msa_info:
        console.print("[dim]No MSA files found.[/]")
        return

    hap_dist = Counter(m.n_haplotypes for m in msa_info)
    ref_lens = [m.ref_length for m in msa_info]
    left_gaps = sum(1 for m in msa_info if m.left_gap > 0)
    right_gaps = sum(1 for m in msa_info if m.right_gap > 0)
    any_gap = sum(1 for m in msa_info if m.left_gap > 0 or m.right_gap > 0)

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Key", style="bold")
    table.add_column("Value")
    table.add_row("MSA components", f"{len(msa_info):,}")
    table.add_row("REF anchor length",
                  f"min={min(ref_lens)}, median={statistics.median(ref_lens):.0f}, "
                  f"max={max(ref_lens)}")
    table.add_row("REF with 5' boundary gaps", f"{left_gaps:,}")
    table.add_row("REF with 3' boundary gaps", f"{right_gaps:,}")
    affected_pct = 100 * any_gap / len(msa_info)
    style = "red" if affected_pct > 5 else ("yellow" if affected_pct > 1 else "green")
    table.add_row("REF with ANY boundary gap",
                  Text(f"{any_gap:,} ({affected_pct:.1f}%)", style=style))
    console.print(Panel(table, title="[bold]MSA Summary[/]", border_style="cyan"))

    ht = Table(title="Haplotype Count Distribution", title_style="bold",
               border_style="dim")
    ht.add_column("Haplotypes", justify="right")
    ht.add_column("Components", justify="right")
    for k in sorted(hap_dist):
        ht.add_row(str(k), f"{hap_dist[k]:,}")
    console.print(ht)


def render_variants(variants: list[dict]) -> None:
    """VCF variant type counts and length distributions."""
    if not variants:
        console.print("[dim]No VCF data available.[/]")
        return

    type_counts = Counter(v["type"] for v in variants)
    total = len(variants)

    vt = Table(title="Variant Counts by Type", title_style="bold",
               border_style="dim")
    vt.add_column("Type")
    vt.add_column("Count", justify="right")
    vt.add_column("Pct", justify="right")
    for vtype in sorted(type_counts):
        pct = 100 * type_counts[vtype] / total
        vt.add_row(vtype, f"{type_counts[vtype]:,}", f"{pct:.1f}%")
    vt.add_row(Text("TOTAL", style="bold"), Text(f"{total:,}", style="bold"), "")
    console.print(vt)

    ins_c = type_counts.get("INS", 0)
    del_c = type_counts.get("DEL", 0)
    if del_c > 0:
        ratio = ins_c / del_c
        style = "green" if 0.8 <= ratio <= 1.5 else "red"
        console.print(f"  INS/DEL ratio: [{style}]{ratio:.2f}[/]  "
                      f"(expected ~1.0–1.2)\n")

    # Length distributions for INS and DEL
    bin_ranges = [(1, 5), (6, 10), (11, 20), (21, 50), (51, 100),
                  (101, 150), (151, 200), (201, 500), (501, 99999)]
    bin_labels = ["1–5", "6–10", "11–20", "21–50", "51–100",
                  "101–150", "151–200", "201–500", ">500"]

    for vtype in ("INS", "DEL"):
        lens = [v["len"] for v in variants if v["type"] == vtype]
        if not lens:
            continue
        lt = Table(title=f"{vtype} Length Distribution", title_style="bold",
                   border_style="dim")
        lt.add_column("Bin")
        lt.add_column("Count", justify="right")
        lt.add_column("Pct", justify="right")
        lt.add_column("Bar", no_wrap=True)
        for (lo, hi), label in zip(bin_ranges, bin_labels):
            count = sum(1 for l in lens if lo <= l <= hi)
            if count == 0:
                continue
            pct = 100 * count / len(lens)
            style = "red" if label == "101–150" and pct > 30 else "white"
            lt.add_row(label, f"{count:,}", f"{pct:.1f}%",
                       make_bar(pct / 100, style=style))
        console.print(lt)


def render_suspects(
    comp_stats: list[ComponentStats],
    ratio_threshold: float,
    variants: list[dict] | None = None,
) -> None:
    """Deep dive into components with low SAMPLE/SHARED coverage ratio.

    These components have sample-only nodes whose coverage is a small fraction
    of the total locus depth (SHARED nodes).  They are candidates for:
      - Sequencing error bubbles that survived pruning
      - Low-frequency mosaic/somatic variants
      - Graph artifacts from misaligned reads in repetitive regions
    """
    suspects = [
        cs for cs in comp_stats
        if cs.med_shared > 0 and cs.sample_shared_ratio < ratio_threshold
    ]
    total = sum(1 for cs in comp_stats if cs.med_shared > 0)

    if not suspects:
        console.print(f"[green]No components below ratio threshold {ratio_threshold}[/]")
        return

    console.print(Panel(
        f"[bold]{len(suspects):,}[/] of {total:,} components "
        f"({100 * len(suspects) / total:.1f}%) have "
        f"SAMPLE/SHARED ratio < {ratio_threshold}",
        title="[bold]Suspect Components[/]", border_style="yellow",
    ))

    # ── K-value distribution for suspects vs all ──────────────────────────
    suspect_k = Counter(cs.k for cs in suspects)
    all_k = Counter(cs.k for cs in comp_stats)
    kt = Table(title="K-value: Suspects vs All Components", title_style="bold",
               border_style="dim")
    kt.add_column("k", justify="right")
    kt.add_column("Suspect", justify="right")
    kt.add_column("All", justify="right")
    kt.add_column("Enrichment", justify="right")
    for k in sorted(set(suspect_k) | set(all_k)):
        s_count = suspect_k.get(k, 0)
        a_count = all_k.get(k, 0)
        s_pct = 100 * s_count / len(suspects) if suspects else 0
        a_pct = 100 * a_count / len(comp_stats) if comp_stats else 0
        enrichment = s_pct / a_pct if a_pct > 0 else 0
        style = "red" if enrichment > 2 else ("yellow" if enrichment > 1.5 else "white")
        kt.add_row(str(k), f"{s_count:,} ({s_pct:.1f}%)",
                   f"{a_count:,} ({a_pct:.1f}%)",
                   Text(f"{enrichment:.2f}x", style=style))
    console.print(kt)

    # ── Coverage breakdown for suspects ───────────────────────────────────
    s_shared = [cs.med_shared for cs in suspects]
    s_sample = [cs.med_sample for cs in suspects]
    s_ref = [cs.med_ref for cs in suspects if cs.ref_covs]
    s_ratios = [cs.sample_shared_ratio for cs in suspects]
    s_n_sample = [len(cs.sample_covs) for cs in suspects]
    s_n_ref = [len(cs.ref_covs) for cs in suspects]

    ct = Table(title="Suspect Coverage Profile", title_style="bold",
               border_style="dim")
    ct.add_column("Metric")
    ct.add_column("Min", justify="right")
    ct.add_column("Median", justify="right")
    ct.add_column("Mean", justify="right")
    ct.add_column("Max", justify="right")
    for name, vals in [
        ("SHARED median cov", s_shared),
        ("SAMPLE median cov", s_sample),
        ("REF median cov", s_ref),
        ("SAMPLE/SHARED ratio", s_ratios),
        ("# SAMPLE nodes", s_n_sample),
        ("# REF nodes", s_n_ref),
    ]:
        if not vals:
            continue
        ct.add_row(name, f"{min(vals):.1f}", f"{statistics.median(vals):.1f}",
                   f"{statistics.mean(vals):.1f}", f"{max(vals):.1f}")
    console.print(ct)

    # ── SHARED coverage bins: are suspects in low-coverage regions? ───────
    shared_bins = [(0, 10), (10, 20), (20, 30), (30, 50), (50, 9999)]
    shared_labels = ["<10", "10–20", "20–30", "30–50", ">50"]
    st = Table(title="Suspect SHARED Coverage Distribution (total locus depth)",
               title_style="bold", border_style="dim")
    st.add_column("SHARED cov")
    st.add_column("Count", justify="right")
    st.add_column("Pct", justify="right")
    st.add_column("Bar", no_wrap=True)
    for (lo, hi), label in zip(shared_bins, shared_labels):
        count = sum(1 for v in s_shared if lo <= v < hi)
        if count == 0:
            continue
        pct = 100 * count / len(s_shared)
        style = "yellow" if lo < 10 else "white"
        st.add_row(label, f"{count:,}", f"{pct:.1f}%",
                   make_bar(pct / 100, style=style))
    console.print(st)

    # ── Cross-reference with VCF: allele support in suspect vs non-suspect ──
    # Use a sorted interval index for O(V log S) instead of O(V × S)
    if variants:
        from bisect import bisect_right

        # Build sorted interval list per chromosome: [(start, end), ...]
        suspect_intervals: dict[str, list[tuple[int, int]]] = {}
        for cs in suspects:
            parts = cs.window.split(":")
            chrom = parts[0]
            start, end = (int(x) for x in parts[1].split("-"))
            suspect_intervals.setdefault(chrom, []).append((start, end))
        for chrom in suspect_intervals:
            suspect_intervals[chrom].sort()

        # Pre-compute sorted start arrays for bisect
        suspect_starts = {
            chrom: [iv[0] for iv in ivs]
            for chrom, ivs in suspect_intervals.items()
        }

        suspect_hits = []
        non_suspect_hits = []
        for v in variants:
            chrom = v["chrom"]
            in_suspect = False
            if chrom in suspect_intervals:
                pos = v["pos"]
                starts = suspect_starts[chrom]
                ivs = suspect_intervals[chrom]
                idx = bisect_right(starts, pos) - 1
                if idx >= 0 and ivs[idx][0] <= pos <= ivs[idx][1]:
                    in_suspect = True
            if in_suspect:
                suspect_hits.append(v)
            else:
                non_suspect_hits.append(v)

        # ── Variant type counts: suspect vs non-suspect ──
        s_types = Counter(v["type"] for v in suspect_hits)
        n_types = Counter(v["type"] for v in non_suspect_hits)
        tt = Table(title="Variant Counts: Suspect vs Non-Suspect Windows",
                   title_style="bold", border_style="dim")
        tt.add_column("Type")
        tt.add_column("Suspect", justify="right")
        tt.add_column("Non-Suspect", justify="right")
        for vtype in sorted(set(s_types) | set(n_types)):
            tt.add_row(vtype, f"{s_types.get(vtype, 0):,}",
                       f"{n_types.get(vtype, 0):,}")
        tt.add_row(Text("TOTAL", style="bold"),
                   Text(f"{len(suspect_hits):,}", style="bold"),
                   Text(f"{len(non_suspect_hits):,}", style="bold"))
        console.print(tt)

        # ── Allele support comparison by variant type ──
        for vtype in ("SNV", "INS", "DEL"):
            s_vars = [v for v in suspect_hits if v["type"] == vtype]
            n_vars = [v for v in non_suspect_hits if v["type"] == vtype]
            if not s_vars and not n_vars:
                continue

            at = Table(
                title=f"{vtype} Allele Support: Suspect vs Non-Suspect",
                title_style="bold", border_style="dim",
            )
            at.add_column("Metric")
            at.add_column("Suspect", justify="right")
            at.add_column("Non-Suspect", justify="right")

            for label, s_vals, n_vals in [
                ("Count", [1] * len(s_vars), [1] * len(n_vars)),
                ("REF AD", [v["ref_ad"] for v in s_vars],
                           [v["ref_ad"] for v in n_vars]),
                ("ALT AD", [v["alt_ad"] for v in s_vars],
                           [v["alt_ad"] for v in n_vars]),
                ("DP", [v["dp"] for v in s_vars],
                       [v["dp"] for v in n_vars]),
                ("VAF", [v["vaf"] for v in s_vars],
                        [v["vaf"] for v in n_vars]),
            ]:
                if label == "Count":
                    at.add_row(label, f"{len(s_vars):,}", f"{len(n_vars):,}")
                    continue
                s_str = (f"med={statistics.median(s_vals):.1f} "
                         f"mean={statistics.mean(s_vals):.1f}") if s_vals else "—"
                n_str = (f"med={statistics.median(n_vals):.1f} "
                         f"mean={statistics.mean(n_vals):.1f}") if n_vals else "—"
                at.add_row(label, s_str, n_str)
            console.print(at)

    # ── Top 20 lowest-ratio components ────────────────────────────────────
    worst = sorted(suspects, key=lambda cs: cs.sample_shared_ratio)[:20]
    wt = Table(title="Top 20 Lowest SAMPLE/SHARED Ratio Components",
               title_style="bold", border_style="dim")
    wt.add_column("Window")
    wt.add_column("Comp", justify="right")
    wt.add_column("k", justify="right")
    wt.add_column("Ratio", justify="right")
    wt.add_column("SHARED", justify="right")
    wt.add_column("SAMPLE", justify="right")
    wt.add_column("REF", justify="right")
    wt.add_column("#S nodes", justify="right")
    for cs in worst:
        wt.add_row(
            cs.window, str(cs.comp), str(cs.k),
            f"{cs.sample_shared_ratio:.3f}",
            f"{cs.med_shared:.0f}", f"{cs.med_sample:.0f}",
            f"{cs.med_ref:.0f}", str(len(cs.sample_covs)),
        )
    console.print(wt)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CLI
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze Lancet2 debug run output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=VIEW_HELP,
    )
    parser.add_argument("--log", type=Path, required=True,
                        help="Path to verbose log file")
    parser.add_argument("--graphs", type=Path, required=True,
                        help="Path to out_graphs directory")
    parser.add_argument("--vcf", type=Path, default=None,
                        help="Path to output VCF (optional)")
    parser.add_argument("--workers", type=int, default=16,
                        help="Number of parallel workers (default: 16)")
    parser.add_argument("--ratio-threshold", type=float, default=0.10,
                        help="SAMPLE/SHARED ratio cutoff for suspect components "
                             "(default: 0.10)")
    parser.add_argument("--view", choices=VIEWS, default="all",
                        help="Report view to render (default: all)")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if not args.log.exists():
        console.print(f"[red]ERROR:[/] Log file not found: {args.log}")
        sys.exit(1)
    if not args.graphs.exists():
        console.print(f"[red]ERROR:[/] Graphs directory not found: {args.graphs}")
        sys.exit(1)

    view = args.view

    # Always parse the log
    console.print("[dim]Parsing log file...[/]")
    windows = parse_log(args.log)

    if view in ("overview", "all"):
        render_overview(windows)

    # comp_stats is needed by both coverage and suspects; load once
    comp_stats = None
    if view in ("coverage", "suspects", "all"):
        console.print("[dim]Loading DOT graph coverage data...[/]")
        comp_stats = load_component_stats(windows, args.graphs, args.workers)

    if view in ("coverage", "all"):
        render_coverage(comp_stats)

    if view in ("msa", "all"):
        console.print("[dim]Loading MSA FASTA files...[/]")
        msa_info = load_msa_info(args.graphs, args.workers)
        render_msa(msa_info)

    variants = None
    if view in ("variants", "suspects", "all"):
        if args.vcf and args.vcf.exists():
            console.print("[dim]Parsing VCF...[/]")
            variants = parse_vcf(args.vcf)
        elif view == "variants":
            console.print("[red]ERROR:[/] --vcf required for variants view")
            sys.exit(1)

    if view in ("variants", "all") and variants:
        render_variants(variants)

    if view in ("suspects", "all") and comp_stats:
        render_suspects(comp_stats, args.ratio_threshold, variants)


if __name__ == "__main__":
    main()
