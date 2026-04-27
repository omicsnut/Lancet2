#!/usr/bin/env python3
"""Analyze Lancet2 probe variant results — forensic pipeline attribution.

Reads probe_results.tsv (from Lancet2 --probe-variants) cross-referenced
against missed_variants.txt and concordance_details.txt (from truth_concordance.py)
to answer: for every missed variant, exactly where and why did the pipeline lose it?

The C++ side emits one raw fact row per (probe_id, window, comp_id, k) attempt
with no attribution. This script is the sole source of lost_at_stage derivation,
using a bottom-up cascade that checks the deepest pipeline evidence first.

Output files:
  probe_analysis_report.txt   Unified rich report (§1–§7, all tables)
  probe_stage_attribution.txt Per-probe TSV: final lost_at, type, size, tier1 reads
  probe_survival_matrix.txt   Per-(probe, window, comp, k) TSV: all 6 survival counts

Usage:
    pixi run -e hts-tools python3 scripts/analyze_probe_results.py \\
        --probe-results data/probe_missing_variants.results.tsv \\
        --missed-variants data/missed_variants.txt \\
        --concordance-details data/concordance_details.txt \\
        --log data/probe_missing_vars.chr1.debug_run.log \\
        --output-dir data/ \\
        --view all
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import polars as pl
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


REPO_ROOT = Path(__file__).resolve().parent.parent

VIEWS = ("scorecard", "funnel", "survival", "breakdown", "genotyper", "targets", "deepdive", "windows", "all")

VIEW_HELP = """\
Views:
  scorecard    Executive summary: coverage validation, vital signs
  funnel       Stage attribution: where variants are lost (21-level cascade)
  survival     Pruning survival: k-mer attrition through 6 graph stages
  breakdown    Type × Size × Stage cross-tabulation
  genotyper    Genotyper forensics: stolen reads, MSA extraction
  targets      High-priority inspection targets (closest to success)
  deepdive     Deep dive into the top 2 most common loss stages
  windows      Window coverage: cross-window attribution, boundary effects
  all          All of the above (default)
"""

# All possible lost_at_stage values — computed in Python from raw TSV facts.
# Ordered by pipeline position for display. Grouped into categories.
STAGE_ORDER = [
    "not_processed",
    "variant_in_anchor", "no_anchor", "short_anchor",
    "graph_has_cycle", "graph_too_complex",
    "pruned_at_build", "pruned_at_lowcov1", "pruned_at_compress1",
    "pruned_at_lowcov2", "pruned_at_compress2", "pruned_at_tips",
    "bfs_exhausted", "no_path",
    "msa_not_extracted", "msa_shifted", "msa_subsumed",
    "geno_no_overlap", "geno_zero_alt_reads", "geno_stolen",
    "survived",
]

STAGE_CATEGORY = {
    "not_processed": "Not processed",
    "variant_in_anchor": "Structural", "no_anchor": "Structural",
    "short_anchor": "Structural", "graph_has_cycle": "Structural",
    "graph_too_complex": "Structural",
    "pruned_at_build": "Pruning", "pruned_at_lowcov1": "Pruning",
    "pruned_at_compress1": "Pruning", "pruned_at_lowcov2": "Pruning",
    "pruned_at_compress2": "Pruning", "pruned_at_tips": "Pruning",
    "bfs_exhausted": "Path enum", "no_path": "Path enum",
    "msa_not_extracted": "MSA", "msa_shifted": "MSA",
    "msa_subsumed": "MSA",
    "geno_no_overlap": "Genotyper", "geno_zero_alt_reads": "Genotyper",
    "geno_stolen": "Genotyper",
    "survived": "Survived",
}

CATEGORY_ORDER = [
    "Not processed", "Structural", "Pruning", "Path enum",
    "MSA", "Genotyper", "Survived",
]

CATEGORY_STYLES = {
    "Not processed": "dim",
    "Structural": "red",
    "Pruning": "yellow",
    "Path enum": "magenta",
    "MSA": "cyan",
    "Genotyper": "blue",
    "Survived": "green",
}

SIZE_TIERS = [
    ("1bp",       1,    1),
    ("2-5bp",     2,    5),
    ("6-20bp",    6,   20),
    ("21-50bp",  21,   50),
    ("51-200bp", 51,  200),
    ("201bp+",  201, 999999),
]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Output Helpers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


class DualConsole:
    """Console wrapper that writes to both a file (clean) and terminal (rich)."""

    def __init__(self, output_path: Path) -> None:
        self._file = open(output_path, "w")
        self._file_console = Console(
            file=self._file, width=200, no_color=True, highlight=False,
        )
        self._term_console = Console(stderr=True)

    def print(self, *args, **kwargs) -> None:
        self._file_console.print(*args, **kwargs)
        self._term_console.print(*args, **kwargs)


def make_bar(fraction: float, width: int = 25, style: str = "white") -> Text:
    """Render a horizontal bar using Unicode block characters."""
    fraction = max(0.0, min(1.0, fraction))
    filled = int(fraction * width)
    half = "▌" if (fraction * width - filled) > 0.5 else ""
    bar_str = f"{'█' * filled}{half}".ljust(width, "░")
    return Text(bar_str, style=style)


def severity_style(pct: float) -> str:
    """Pick a rich style based on percentage severity."""
    if pct >= 20.0:
        return "red"
    if pct >= 10.0:
        return "yellow"
    if pct >= 2.0:
        return "white"
    return "dim"


def pct_str(count: int, total: int) -> str:
    """Format a percentage string, guarding against division by zero."""
    return f"{count / total * 100:.1f}%" if total > 0 else "—"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Data Loading
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def load_probe_results(path: Path) -> pl.DataFrame:
    """Load probe_results.tsv — multi-row per probe (one per k-value)."""
    return pl.read_csv(path, separator="\t", infer_schema_length=5000)


def _load_commented_tsv(path: Path) -> pl.DataFrame:
    """Load a TSV file whose header line starts with '#'.

    Reads the first line to extract column names (stripping the leading '#'),
    then reads the remaining data rows with those names. This avoids polars'
    comment_prefix eating the header.
    """
    with open(path) as fh:
        header_line = fh.readline().strip()
    columns = [c.lstrip("#") for c in header_line.split("\t")]
    return pl.read_csv(
        path, separator="\t",
        has_header=False, skip_rows=1, new_columns=columns,
        infer_schema_length=5000,
    )


def load_missed_variants(path: Path) -> pl.DataFrame:
    """Load missed_variants.txt — one row per missed truth variant."""
    return _load_commented_tsv(path)


def load_concordance_details(path: Path) -> pl.DataFrame:
    """Load concordance_details.txt — one row per truth variant (all levels)."""
    return _load_commented_tsv(path)


_RE_WINDOW_STATUS = re.compile(
    r"(chr\w+):(\d+)-(\d+) done with (\w+)"
)


def load_debug_log_window_statuses(path: Path) -> pl.DataFrame:
    """Parse window completion statuses from the Lancet2 debug log.

    Extracts lines matching 'chrX:START-END done with STATUS_CODE'
    and returns a DataFrame with columns: chrom, start, end, status.
    """
    statuses: list[dict] = []
    with open(path) as fh:
        for line in fh:
            if "done with" not in line:
                continue
            match = _RE_WINDOW_STATUS.search(line)
            if match:
                statuses.append({
                    "chrom": match.group(1),
                    "win_start": int(match.group(2)),
                    "win_end": int(match.group(3)),
                    "win_status": match.group(4),
                })
    return pl.DataFrame(statuses) if statuses else pl.DataFrame(
        schema={"chrom": pl.Utf8, "win_start": pl.Int64,
                "win_end": pl.Int64, "win_status": pl.Utf8}
    )


def derive_lost_at(df: pl.DataFrame) -> pl.DataFrame:
    """Derive lost_at_stage for every row from observed diagnostic facts.

    Bottom-up cascade: check deepest pipeline evidence first. The first
    matching condition wins. This replaces the old C++ DeriveLostAt which
    used a top-down priority cascade and operated on a single "best" record.

    Genotyper outcomes checked first (deepest), then MSA, paths, pruning
    stages (in reverse pipeline order), and finally structural flags.
    """
    return df.with_columns(
        # ── Genotyper (deepest) ──────────────────────────────────────
        pl.when(pl.col("is_geno_has_alt_support") == 1)
          .then(pl.lit("survived"))
        .when((pl.col("is_msa_exact_match") == 1) & (pl.col("is_geno_no_overlap") == 1))
          .then(pl.lit("geno_no_overlap"))
        .when((pl.col("is_msa_exact_match") == 1) & (pl.col("is_geno_has_alt_support") == 0)
              & (pl.col("n_geno_stolen_to_ref") + pl.col("n_geno_stolen_to_wrong_alt") > 0))
          .then(pl.lit("geno_stolen"))
        .when((pl.col("is_msa_exact_match") == 1) & (pl.col("is_geno_has_alt_support") == 0))
          .then(pl.lit("geno_zero_alt_reads"))
        # ── MSA extraction ───────────────────────────────────────────
        .when((pl.col("hap_indices") != ".") & (pl.col("is_msa_shifted") == 1))
          .then(pl.lit("msa_shifted"))
        .when((pl.col("hap_indices") != ".") & (pl.col("is_msa_subsumed") == 1))
          .then(pl.lit("msa_subsumed"))
        .when(pl.col("hap_indices") != ".")
          .then(pl.lit("msa_not_extracted"))
        # ── Path enumeration ─────────────────────────────────────────
        .when((pl.col("n_surviving_tips") > 0) & (pl.col("is_traversal_limited") == 1))
          .then(pl.lit("bfs_exhausted"))
        .when(pl.col("n_surviving_tips") > 0)
          .then(pl.lit("no_path"))
        # ── Pruning (reverse pipeline order: first stage that zeroed) ─
        .when((pl.col("n_surviving_compress2") > 0) & (pl.col("n_surviving_tips") == 0))
          .then(pl.lit("pruned_at_tips"))
        .when((pl.col("n_surviving_lowcov2") > 0) & (pl.col("n_surviving_compress2") == 0))
          .then(pl.lit("pruned_at_compress2"))
        .when((pl.col("n_surviving_compress1") > 0) & (pl.col("n_surviving_lowcov2") == 0))
          .then(pl.lit("pruned_at_lowcov2"))
        .when((pl.col("n_surviving_lowcov1") > 0) & (pl.col("n_surviving_compress1") == 0))
          .then(pl.lit("pruned_at_compress1"))
        .when((pl.col("n_surviving_build") > 0) & (pl.col("n_surviving_lowcov1") == 0))
          .then(pl.lit("pruned_at_lowcov1"))
        .when((pl.col("n_expected_alt_kmers") > 0) & (pl.col("n_surviving_build") == 0))
          .then(pl.lit("pruned_at_build"))
        # ── Structural (shallowest — only when nothing above fired) ──
        .when(pl.col("is_graph_cycle") == 1).then(pl.lit("graph_has_cycle"))
        .when(pl.col("is_graph_complex") == 1).then(pl.lit("graph_too_complex"))
        .when(pl.col("is_no_anchor") == 1).then(pl.lit("no_anchor"))
        .when(pl.col("is_short_anchor") == 1).then(pl.lit("short_anchor"))
        .when(pl.col("is_variant_in_anchor") == 1).then(pl.lit("variant_in_anchor"))
        .otherwise(pl.lit("not_processed"))
        .alias("lost_at_stage")
    )


# Depth scores: higher = deeper in the pipeline = more interesting.
# Used to select the "best" record per probe_id.
_DEPTH_MAP = {s: i for i, s in enumerate(STAGE_ORDER)}


def compute_depth(df: pl.DataFrame) -> pl.DataFrame:
    """Add _depth column: integer score for how deep each attempt reached."""
    return df.with_columns(
        pl.col("lost_at_stage")
        .replace_strict(_DEPTH_MAP, default=0)
        .alias("_depth")
    )


def select_best_per_probe(df: pl.DataFrame) -> pl.DataFrame:
    """One row per probe_id — the attempt that reached the deepest stage.

    When multiple records tie on depth, the highest kmer_size wins.
    This replaces the old build_final_attribution which blindly took the
    highest-k record regardless of how far it actually progressed.
    """
    return (
        df
        .sort(["_depth", "kmer_size"], descending=[True, True])
        .group_by("probe_id")
        .first()
        .sort("probe_id")
    )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §1 Executive Scorecard
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_scorecard(
    console: DualConsole,
    probes: pl.DataFrame,
    attribution: pl.DataFrame,
    missed: pl.DataFrame,
    concordance: Optional[pl.DataFrame],
    window_statuses: Optional[pl.DataFrame],
) -> None:
    """§1: Executive scorecard — coverage validation and vital signs."""
    console.print("[bold]§1 Executive Scorecard[/bold]\n")

    n_input = missed.height
    n_rows = probes.height
    n_unique_probes = probes["probe_id"].n_unique()
    n_attribution = attribution.height

    # Coverage gap: probe IDs in missed but not in probe_results
    coverage_gap = n_input - n_unique_probes
    gap_style = "green" if coverage_gap == 0 else "red bold"
    gap_text = f"[{gap_style}]{coverage_gap}[/]"
    if coverage_gap == 0:
        gap_text += f" [{gap_style}]✓[/]"

    k_values = sorted(probes["kmer_size"].unique().to_list())
    k_range_str = f"{len(k_values)} ({min(k_values)}..{max(k_values)})" if k_values else "—"

    # Stage distribution for the scorecard
    stage_counts = attribution["lost_at_stage"].value_counts().sort("count", descending=True)

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Key", style="bold", min_width=35)
    table.add_column("Value")

    # Concordance context (if available)
    if concordance is not None:
        n_truth = concordance.height
        level_counts = concordance.get_column(concordance.columns[7]).value_counts()
        n_matched = n_truth - n_input
        table.add_row("Truth variants (concordance_details)", f"{n_truth:,}")
        table.add_row("Matched by Lancet2 (L0+LD+L1-L3)",
                       f"[green]{n_matched:,}[/] ({pct_str(n_matched, n_truth)})")
        table.add_row("Missed (input to probe system)",
                       f"[yellow]{n_input:,}[/] ({pct_str(n_input, n_truth)})")
        table.add_row("", "")

    table.add_row("Input variants (missed_variants.txt)", f"{n_input:,}")
    table.add_row("Probe results rows (multi-k)", f"{n_rows:,}")
    table.add_row("Unique probe IDs in results", f"{n_unique_probes:,}")
    table.add_row("Coverage gap (missing probes)", gap_text)
    table.add_row("Distinct k-values observed", k_range_str)

    # Top 5 stages
    table.add_row("", "")
    for row in stage_counts.head(5).iter_rows(named=True):
        stage = row[stage_counts.columns[0]]
        count = row["count"]
        cat = STAGE_CATEGORY.get(stage, "Other")
        style = CATEGORY_STYLES.get(cat, "white")
        table.add_row(f"  {stage}", f"[{style}]{count:,}[/] ({pct_str(count, n_attribution)})")

    # Window status summary from debug log
    if window_statuses is not None and window_statuses.height > 0:
        table.add_row("", "")
        ws_counts = window_statuses["win_status"].value_counts().sort("count", descending=True)
        n_windows = window_statuses.height
        table.add_row("Total windows in debug log", f"{n_windows:,}")
        for row in ws_counts.head(5).iter_rows(named=True):
            status = row[ws_counts.columns[0]]
            count = row["count"]
            table.add_row(f"  {status}", f"{count:,} ({pct_str(count, n_windows)})")

    console.print(Panel(table, title="[bold]Probe Forensics Overview[/]", border_style="cyan"))
    console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §2 Stage Attribution Funnel
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_stage_funnel(
    console: DualConsole,
    attribution: pl.DataFrame,
    output_dir: Path,
) -> None:
    """§2: Stage attribution funnel — where variants are lost (21-level cascade)."""
    console.print("[bold]§2 Stage Attribution Funnel[/bold]\n")

    total = attribution.height
    stage_counts = dict(
        attribution["lost_at_stage"]
        .value_counts()
        .iter_rows()
    )

    # Per-stage table in pipeline order
    table = Table(
        title=f"Pipeline Attribution ({total:,} probes)",
        title_style="bold", border_style="dim",
    )
    table.add_column("Stage")
    table.add_column("Category")
    table.add_column("Count", justify="right")
    table.add_column("Pct", justify="right")
    table.add_column("Cum%", justify="right")
    table.add_column("Bar", no_wrap=True)

    cumulative = 0
    for stage in STAGE_ORDER:
        count = stage_counts.get(stage, 0)
        if count == 0:
            continue
        cumulative += count
        pct = count / total * 100
        cum_pct = cumulative / total * 100
        cat = STAGE_CATEGORY.get(stage, "Other")
        style = CATEGORY_STYLES.get(cat, "white")
        bar = make_bar(pct / 100, width=25, style=style)
        table.add_row(
            Text(stage, style=style), Text(cat, style=style),
            Text(f"{count:,}", style=style), Text(f"{pct:.1f}%", style=style),
            f"{cum_pct:.1f}%", bar,
        )

    console.print(table)
    console.print()

    # Category-level summary
    cat_counts: dict[str, int] = {}
    for stage in STAGE_ORDER:
        count = stage_counts.get(stage, 0)
        cat = STAGE_CATEGORY.get(stage, "Other")
        cat_counts[cat] = cat_counts.get(cat, 0) + count

    cat_table = Table(
        title="Category Summary", title_style="bold", border_style="dim",
    )
    cat_table.add_column("Category")
    cat_table.add_column("Count", justify="right")
    cat_table.add_column("Pct", justify="right")
    cat_table.add_column("Bar", no_wrap=True)

    for cat in CATEGORY_ORDER:
        count = cat_counts.get(cat, 0)
        if count == 0:
            continue
        pct = count / total * 100
        style = CATEGORY_STYLES.get(cat, "white")
        bar = make_bar(pct / 100, width=25, style=style)
        cat_table.add_row(
            Text(cat, style=style), Text(f"{count:,}", style=style),
            Text(f"{pct:.1f}%", style=style), bar,
        )

    console.print(cat_table)
    console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §3 Pruning Survival Funnel
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SURVIVAL_COLS = [
    ("n_surviving_build", "build"),
    ("n_surviving_lowcov1", "lowcov1"),
    ("n_surviving_compress1", "compress1"),
    ("n_surviving_lowcov2", "lowcov2"),
    ("n_surviving_compress2", "compress2"),
    ("n_surviving_tips", "tips"),
]


def render_survival_analysis(
    console: DualConsole,
    probes: pl.DataFrame,
) -> None:
    """§3: Pruning survival — k-mer attrition through 6 graph stages."""
    console.print("[bold]§3 Pruning Survival Funnel[/bold]\n")

    # Only analyze probes that were actually processed (kmer_size > 0)
    processed = probes.filter(pl.col("kmer_size") > 0)
    if processed.height == 0:
        console.print("[dim]No processed probes to analyze.[/]")
        return

    n_probes = processed.height

    table = Table(
        title=f"K-mer Survival Through Pruning ({n_probes:,} probe×k records)",
        title_style="bold", border_style="dim",
    )
    table.add_column("Stage")
    table.add_column("Median", justify="right")
    table.add_column("Mean", justify="right")
    table.add_column("Min", justify="right")
    table.add_column("Max", justify="right")
    table.add_column("All-zero", justify="right")
    table.add_column("Drop→0 %", justify="right")

    prev_zero = 0
    for col, label in SURVIVAL_COLS:
        vals = processed[col]
        median = vals.median()
        mean = vals.mean()
        minimum = vals.min()
        maximum = vals.max()
        zero_count = (vals == 0).sum()
        new_zeros = zero_count - prev_zero
        drop_pct = new_zeros / n_probes * 100

        style = "red" if drop_pct > 10 else ("yellow" if drop_pct > 2 else "white")
        table.add_row(
            label,
            f"{median:.0f}", f"{mean:.1f}",
            str(minimum), str(maximum),
            f"{zero_count:,}",
            Text(f"{drop_pct:.1f}%", style=style),
        )
        prev_zero = zero_count

    console.print(table)
    console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §4 Type × Size × Stage Cross-Tabulation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_type_size_breakdown(
    console: DualConsole,
    attribution: pl.DataFrame,
    missed: pl.DataFrame,
) -> None:
    """§4: Type × Size × Stage cross-tabulation."""
    console.print("[bold]§4 Type × Size × Stage Breakdown[/bold]\n")

    # Join attribution with missed variants to get type and length.
    # Use (pos, ref, alt) as composite key to avoid duplicate rows
    # when multiple missed variants share the same genomic position.
    joined = attribution.join(
        missed.select("start", "ref", "alt", "type", "length").rename({"ref": "m_ref", "alt": "m_alt"}),
        left_on=["pos", "ref", "alt"], right_on=["start", "m_ref", "m_alt"], how="left",
    )

    # Per-type stage distribution
    types = ["SNV", "INS", "DEL", "MNP", "CPX"]
    all_stages = ["not_processed", "no_anchor", "variant_in_anchor",
                  "graph_has_cycle", "pruned_at_build", "pruned_at_tips",
                  "no_path", "msa_subsumed", "geno_zero_alt_reads", "survived"]

    for vtype in types:
        vt_data = joined.filter(pl.col("type") == vtype)
        if vt_data.height == 0:
            continue

        # Only include stages with non-zero counts to keep table width manageable
        full_dist = dict(vt_data["lost_at_stage"].value_counts().iter_rows())
        active_stages = [s for s in all_stages if full_dist.get(s, 0) > 0]

        table = Table(
            title=f"{vtype} — Stage × Size ({vt_data.height:,} probes)",
            title_style="bold", border_style="dim",
        )
        table.add_column("Size")
        for stage in active_stages:
            table.add_column(stage, justify="right", no_wrap=True)
        table.add_column("Other", justify="right", no_wrap=True)
        table.add_column("Total", justify="right", no_wrap=True)

        for tier_name, lo, hi in SIZE_TIERS:
            tier = vt_data.filter(
                (pl.col("length").abs() >= lo) & (pl.col("length").abs() <= hi)
            )
            if tier.height == 0:
                continue
            stage_dist = dict(tier["lost_at_stage"].value_counts().iter_rows())
            row = [tier_name]
            other = tier.height
            for stage in active_stages:
                count = stage_dist.get(stage, 0)
                other -= count
                row.append(str(count) if count > 0 else "·")
            row.append(str(other) if other > 0 else "·")
            row.append(str(tier.height))
            table.add_row(*row)

        # Totals row
        row = ["TOTAL"]
        other = vt_data.height
        for stage in active_stages:
            count = full_dist.get(stage, 0)
            other -= count
            row.append(str(count) if count > 0 else "·")
        row.append(str(other) if other > 0 else "·")
        row.append(str(vt_data.height))
        table.add_row(*row, style="bold")

        console.print(table)
        console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §5 Genotyper Forensics
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_genotyper_forensics(
    console: DualConsole,
    attribution: pl.DataFrame,
) -> None:
    """§5: Genotyper forensics — MSA extraction and read assignment analysis."""
    console.print("[bold]§5 Genotyper Forensics[/bold]\n")

    # MSA extraction outcomes
    msa_stages = ["msa_not_extracted", "msa_shifted", "msa_subsumed"]
    geno_stages = ["geno_no_overlap", "geno_zero_alt_reads", "geno_stolen", "survived"]
    relevant = attribution.filter(
        pl.col("lost_at_stage").is_in(msa_stages + geno_stages)
    )

    if relevant.height == 0:
        console.print("[dim]No probes reached MSA/Genotyper stages.[/]")
        return

    stage_dist = dict(relevant["lost_at_stage"].value_counts().iter_rows())

    table = Table(title="MSA + Genotyper Outcomes", title_style="bold", border_style="dim")
    table.add_column("Stage", min_width=20)
    table.add_column("Count", justify="right", width=8)
    table.add_column("Pct", justify="right", width=7)
    table.add_column("Bar", width=26, no_wrap=True)

    total = relevant.height
    for stage in msa_stages + geno_stages:
        count = stage_dist.get(stage, 0)
        if count == 0:
            continue
        pct = count / total * 100
        cat = STAGE_CATEGORY.get(stage, "Other")
        style = CATEGORY_STYLES.get(cat, "white")
        bar = make_bar(pct / 100, width=25, style=style)
        table.add_row(
            Text(stage, style=style), Text(f"{count:,}", style=style),
            Text(f"{pct:.1f}%", style=style), bar,
        )

    console.print(table)
    console.print()

    # Read assignment stats for genotyper-stage probes
    geno_probes = attribution.filter(
        pl.col("lost_at_stage").is_in(geno_stages)
    )
    if geno_probes.height == 0:
        return

    geno_cols = [
        ("n_geno_true_alt_reads", "True ALT reads"),
        ("n_geno_total_ref_reads", "Total REF reads"),
        ("n_geno_stolen_to_ref", "Stolen to REF"),
        ("n_geno_stolen_to_wrong_alt", "Stolen to wrong ALT"),
        ("n_geno_non_overlapping", "Non-overlapping"),
    ]

    stats_table = Table(
        title=f"Genotyper Read Assignment ({geno_probes.height:,} probes)",
        title_style="bold", border_style="dim",
    )
    stats_table.add_column("Metric", min_width=22)
    stats_table.add_column("Median", justify="right", width=8)
    stats_table.add_column("Mean", justify="right", width=8)
    stats_table.add_column("P95", justify="right", width=6)
    stats_table.add_column("Max", justify="right", width=6)

    for col, label in geno_cols:
        vals = geno_probes[col]
        stats_table.add_row(
            label,
            f"{vals.median():.0f}", f"{vals.mean():.1f}",
            f"{vals.quantile(0.95):.0f}", f"{vals.max()}",
        )

    console.print(stats_table)
    console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §6 High-Priority Inspection Targets
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Priority of stages — higher = closer to success, more interesting to inspect.
STAGE_PRIORITY = {s: i for i, s in enumerate(STAGE_ORDER)}


def render_inspection_targets(
    console: DualConsole,
    attribution: pl.DataFrame,
    missed: pl.DataFrame,
) -> None:
    """§6: High-priority inspection targets — variants closest to success."""
    console.print("[bold]§6 High-Priority Inspection Targets[/bold]\n")

    # Join to get type, tier1 reads (composite key to avoid duplicates)
    joined = attribution.join(
        missed.select("start", "ref", "alt", "type", "length", "tier1_alt_count").rename({"ref": "m_ref", "alt": "m_alt"}),
        left_on=["pos", "ref", "alt"], right_on=["start", "m_ref", "m_alt"], how="left",
    )

    # Exclude survived and not_processed — focus on real losses
    targets = joined.filter(
        ~pl.col("lost_at_stage").is_in(["survived", "not_processed"])
    )

    if targets.height == 0:
        console.print("[dim]No non-trivial losses to inspect.[/]")
        return

    # Add priority score for sorting: higher stage priority + more tier1 reads
    targets = targets.with_columns(
        pl.col("lost_at_stage").replace_strict(STAGE_PRIORITY, default=0).alias("_priority"),
    ).sort(["_priority", "n_tier1_reads"], descending=[True, True])

    top_n = min(30, targets.height)
    top = targets.head(top_n)

    table = Table(
        title=f"Top {top_n} Inspection Targets (closest to success, most read evidence)",
        title_style="bold", border_style="dim",
    )
    table.add_column("#", justify="right")
    table.add_column("Chrom")
    table.add_column("Pos", justify="right")
    table.add_column("Type")
    table.add_column("Size", justify="right")
    table.add_column("T1", justify="right")
    table.add_column("Lost at")
    table.add_column("Tips", justify="right")
    table.add_column("MSA?")
    table.add_column("REF/ALT")

    for i, row in enumerate(top.iter_rows(named=True), 1):
        lost = row["lost_at_stage"]
        cat = STAGE_CATEGORY.get(lost, "Other")
        style = CATEGORY_STYLES.get(cat, "white")
        vtype = row.get("type", "?") or "?"
        vlen = row.get("length", 0) or 0
        ref = (row.get("ref", ".") or ".")[:15]
        alt = (row.get("alt", ".") or ".")[:15]
        msa = "✓" if row.get("is_msa_exact_match", 0) else "·"

        table.add_row(
            str(i), row["chrom"], f"{row['pos']:,}",
            vtype, str(abs(vlen)),
            str(row.get("n_tier1_reads", 0)),
            Text(lost, style=style),
            str(row.get("n_surviving_tips", 0)),
            msa, f"{ref}/{alt}",
        )

    console.print(table)
    console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §7 Top Loss Stage Deep Dive
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_top_loss_deep_dive(
    console: DualConsole,
    attribution: pl.DataFrame,
    missed: pl.DataFrame,
) -> None:
    """§7: Deep dive into the top 2 most common loss stages."""
    console.print("[bold]§7 Top Loss Stage Deep Dive[/bold]\n")

    # Exclude survived and not_processed — focus on actual pipeline losses
    losses = attribution.filter(
        ~pl.col("lost_at_stage").is_in(["survived", "not_processed"])
    )
    if losses.height == 0:
        console.print("[dim]No pipeline losses to analyze.[/]")
        return

    # Find the top 2 stages by count
    stage_counts = (
        losses["lost_at_stage"]
        .value_counts()
        .sort("count", descending=True)
        .head(2)
    )
    top_stages = stage_counts[stage_counts.columns[0]].to_list()

    # Join with missed variants for type/length/tier1 data
    joined = losses.join(
        missed.select("start", "ref", "alt", "type", "length",
                       "tier1_alt_count", "classification")
        .rename({"ref": "m_ref", "alt": "m_alt"}),
        left_on=["pos", "ref", "alt"],
        right_on=["start", "m_ref", "m_alt"],
        how="left",
    )

    for stage in top_stages:
        stage_data = joined.filter(pl.col("lost_at_stage") == stage)
        n = stage_data.height
        cat = STAGE_CATEGORY.get(stage, "Other")
        style = CATEGORY_STYLES.get(cat, "white")

        console.print(
            Panel(
                f"[bold]{n:,}[/] probes ({pct_str(n, attribution.height)} of all probes)",
                title=f"[bold {style}]{stage}[/]",
                border_style=style,
            )
        )

        # ── Type breakdown ────────────────────────────────────────────────
        type_dist = stage_data["type"].value_counts().sort("count", descending=True)
        type_table = Table(
            title="Variant Type Breakdown", title_style="bold", border_style="dim",
        )
        type_table.add_column("Type")
        type_table.add_column("Count", justify="right")
        type_table.add_column("Pct", justify="right")
        type_table.add_column("Bar", no_wrap=True)

        for row in type_dist.iter_rows(named=True):
            vtype = row[type_dist.columns[0]]
            count = row["count"]
            pct = count / n * 100
            bar = make_bar(pct / 100, width=20, style=style)
            type_table.add_row(vtype, f"{count:,}", f"{pct:.1f}%", bar)
        console.print(type_table)

        # ── Size tier breakdown ───────────────────────────────────────────
        size_table = Table(
            title="Size Tier Breakdown", title_style="bold", border_style="dim",
        )
        size_table.add_column("Size")
        size_table.add_column("Count", justify="right")
        size_table.add_column("Pct", justify="right")
        size_table.add_column("Bar", no_wrap=True)

        for tier_name, lo, hi in SIZE_TIERS:
            tier = stage_data.filter(
                (pl.col("length").abs() >= lo) & (pl.col("length").abs() <= hi)
            )
            if tier.height == 0:
                continue
            pct = tier.height / n * 100
            bar = make_bar(pct / 100, width=20, style=style)
            size_table.add_row(tier_name, f"{tier.height:,}", f"{pct:.1f}%", bar)
        console.print(size_table)

        # ── Tier-1 read support distribution ──────────────────────────────
        t1 = stage_data["tier1_alt_count"]
        t1_numeric = t1.cast(pl.Int64, strict=False).fill_null(0)

        read_table = Table(
            title="Tier-1 ALT Read Support", title_style="bold", border_style="dim",
        )
        read_table.add_column("Metric")
        read_table.add_column("Value", justify="right")

        read_table.add_row("Median", f"{t1_numeric.median():.0f}")
        read_table.add_row("Mean", f"{t1_numeric.mean():.1f}")
        read_table.add_row("P25", f"{t1_numeric.quantile(0.25):.0f}")
        read_table.add_row("P75", f"{t1_numeric.quantile(0.75):.0f}")
        read_table.add_row("Max", f"{t1_numeric.max()}")

        # Read support bins
        bins = [(0, 0, "0"), (1, 3, "1-3"), (4, 10, "4-10"),
                (11, 25, "11-25"), (26, 999999, "26+")]
        for label, lo, hi in [(b[2], b[0], b[1]) for b in bins]:
            count = ((t1_numeric >= lo) & (t1_numeric <= hi)).sum()
            read_table.add_row(f"  reads {label}", f"{count:,} ({pct_str(count, n)})")

        console.print(read_table)

        # ── Classification breakdown ──────────────────────────────────────
        cls_dist = stage_data["classification"].value_counts().sort("count", descending=True)
        cls_table = Table(
            title="Truth Classification", title_style="bold", border_style="dim",
        )
        cls_table.add_column("Classification")
        cls_table.add_column("Count", justify="right")
        cls_table.add_column("Pct", justify="right")

        for row in cls_dist.iter_rows(named=True):
            cls = row[cls_dist.columns[0]]
            count = row["count"]
            cls_table.add_row(cls, f"{count:,}", pct_str(count, n))
        console.print(cls_table)

        # ── Top 15 example variants ───────────────────────────────────────
        examples = stage_data.sort("n_tier1_reads", descending=True).head(15)

        ex_table = Table(
            title=f"Top 15 {stage} Variants (by tier-1 reads)",
            title_style="bold", border_style="dim",
        )
        ex_table.add_column("#", justify="right")
        ex_table.add_column("Chrom")
        ex_table.add_column("Pos", justify="right")
        ex_table.add_column("Type")
        ex_table.add_column("Size", justify="right")
        ex_table.add_column("T1 reads", justify="right")
        ex_table.add_column("Classification")
        ex_table.add_column("k", justify="right")
        ex_table.add_column("Tips", justify="right")
        ex_table.add_column("REF/ALT")

        for i, row in enumerate(examples.iter_rows(named=True), 1):
            ref = (row.get("ref", ".") or ".")[:12]
            alt = (row.get("alt", ".") or ".")[:12]
            vtype = row.get("type", "?") or "?"
            vlen = row.get("length", 0) or 0
            ex_table.add_row(
                str(i), row["chrom"], f"{row['pos']:,}",
                vtype, str(abs(vlen)),
                str(row.get("n_tier1_reads", 0)),
                str(row.get("classification", "?")),
                str(row.get("kmer_size", 0)),
                str(row.get("n_surviving_tips", 0)),
                f"{ref}/{alt}",
            )

        console.print(ex_table)
        console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# §8 Window Coverage Analysis
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_window_analysis(
    console: DualConsole,
    probes: pl.DataFrame,
    attribution: pl.DataFrame,
) -> None:
    """§8: Window coverage — cross-window attribution and boundary effects."""
    console.print("[bold]§8 Window Coverage Analysis[/bold]\n")

    if "window" not in probes.columns or probes.height == 0:
        console.print("[dim]No window data in probe results.[/]")
        return

    # ── Attempt distribution per probe ────────────────────────────────────
    attempts = probes.group_by("probe_id").agg(
        pl.len().alias("n_attempts"),
        pl.col("window").n_unique().alias("n_windows"),
        pl.col("comp_id").n_unique().alias("n_components"),
        pl.col("kmer_size").n_unique().alias("n_k_values"),
        pl.col("_depth").max().alias("max_depth"),
        pl.col("_depth").min().alias("min_depth"),
    )

    n_probes = attempts.height
    table = Table(
        title=f"Attempt Distribution ({n_probes:,} probes)",
        title_style="bold", border_style="dim",
    )
    table.add_column("Metric", min_width=35)
    table.add_column("Value", justify="right")

    for col, label in [
        ("n_attempts", "Attempts per probe"),
        ("n_windows", "Distinct windows per probe"),
        ("n_components", "Distinct components per probe"),
        ("n_k_values", "Distinct k-values per probe"),
    ]:
        vals = attempts[col]
        table.add_row(
            f"  {label} (median / mean / max)",
            f"{vals.median():.0f} / {vals.mean():.1f} / {vals.max()}",
        )

    console.print(table)
    console.print()

    # ── Boundary sensitivity: probes with mixed outcomes ──────────────────
    # survived_depth is the depth score for "survived"
    survived_depth = _DEPTH_MAP.get("survived", len(STAGE_ORDER) - 1)
    mixed = attempts.filter(
        (pl.col("max_depth") >= survived_depth) & (pl.col("min_depth") < survived_depth)
    )
    n_mixed = mixed.height

    boundary_table = Table(
        title="Boundary Sensitivity", title_style="bold", border_style="dim",
    )
    boundary_table.add_column("Metric", min_width=45)
    boundary_table.add_column("Value", justify="right")
    boundary_table.add_row(
        "Probes that survived in ≥1 attempt but failed in another",
        f"{n_mixed:,} ({pct_str(n_mixed, n_probes)})",
    )
    if n_mixed > 0:
        boundary_table.add_row(
            "  → These probes are boundary-sensitive",
            "[yellow]window placement affects outcome[/]",
        )

    console.print(boundary_table)
    console.print()

    # ── Per-window pathology hotspots ─────────────────────────────────────
    window_stats = (
        probes
        .filter(pl.col("window") != ".")
        .group_by("window")
        .agg(
            pl.len().alias("n_records"),
            pl.col("probe_id").n_unique().alias("n_probes"),
            (pl.col("lost_at_stage") == "survived").sum().alias("n_survived"),
            (pl.col("lost_at_stage").is_in(
                ["graph_has_cycle", "graph_too_complex", "no_anchor", "short_anchor"]
            )).sum().alias("n_structural_failures"),
        )
        .with_columns(
            (pl.col("n_survived") / pl.col("n_records") * 100).alias("survival_pct"),
        )
        .sort("n_structural_failures", descending=True)
    )

    if window_stats.height > 0:
        top_n = min(15, window_stats.height)
        top_windows = window_stats.head(top_n)

        win_table = Table(
            title=f"Top {top_n} Windows by Structural Failures",
            title_style="bold", border_style="dim",
        )
        win_table.add_column("Window")
        win_table.add_column("Records", justify="right")
        win_table.add_column("Probes", justify="right")
        win_table.add_column("Survived", justify="right")
        win_table.add_column("Structural", justify="right")
        win_table.add_column("Surv%", justify="right")

        for row in top_windows.iter_rows(named=True):
            style = "red" if row["n_structural_failures"] > 5 else "white"
            win_table.add_row(
                row["window"],
                str(row["n_records"]),
                str(row["n_probes"]),
                str(row["n_survived"]),
                Text(str(row["n_structural_failures"]), style=style),
                f"{row['survival_pct']:.0f}%",
            )

        console.print(win_table)
        console.print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Output Writers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def write_stage_attribution_tsv(
    output_dir: Path,
    attribution: pl.DataFrame,
    missed: pl.DataFrame,
) -> None:
    """Write per-probe final attribution TSV joined with variant metadata."""
    joined = attribution.select(
        "probe_id", "chrom", "pos", "ref", "alt",
        "n_tier1_reads", "kmer_size", "lost_at_stage",
        "n_surviving_tips", "is_msa_exact_match",
        "n_geno_true_alt_reads", "n_geno_stolen_to_ref",
    ).join(
        missed.select("start", "ref", "alt", "type", "length", "classification").rename({"ref": "m_ref", "alt": "m_alt"}),
        left_on=["pos", "ref", "alt"], right_on=["start", "m_ref", "m_alt"], how="left",
    )
    path = output_dir / "probe_stage_attribution.txt"
    joined.write_csv(path, separator="\t")
    print(f"  Wrote {joined.height} rows to {path}", file=sys.stderr)


def write_survival_matrix_tsv(
    output_dir: Path,
    probes: pl.DataFrame,
) -> None:
    """Write per-(probe_id, kmer_size) survival matrix TSV."""
    cols = ["probe_id", "chrom", "pos", "ref", "alt", "window", "comp_id",
            "kmer_size", "n_expected_alt_kmers", "n_alt_kmers_in_reads",
            "n_surviving_build", "n_surviving_lowcov1",
            "n_surviving_compress1", "n_surviving_lowcov2",
            "n_surviving_compress2", "n_surviving_tips",
            "lost_at_stage"]
    matrix = probes.select([c for c in cols if c in probes.columns])
    path = output_dir / "probe_survival_matrix.txt"
    matrix.write_csv(path, separator="\t")
    print(f"  Wrote {matrix.height} rows to {path}", file=sys.stderr)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CLI
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze Lancet2 probe variant results — forensic pipeline attribution.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=VIEW_HELP,
    )
    parser.add_argument("--probe-results", type=Path, required=True,
                        help="Path to probe_results.tsv from Lancet2 probe run")
    parser.add_argument("--missed-variants", type=Path, required=True,
                        help="Path to missed_variants.txt from truth_concordance.py")
    parser.add_argument("--concordance-details", type=Path, default=None,
                        help="Path to concordance_details.txt (optional, enriches §1)")
    parser.add_argument("--log", type=Path, default=None,
                        help="Path to Lancet2 debug log (enables window status analysis)")
    parser.add_argument("--output-dir", type=Path, default=Path("data"),
                        help="Directory for output files (default: data/)")
    parser.add_argument("--view", choices=VIEWS, default="all",
                        help="Report section to render (default: all)")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    for path, label in [
        (args.probe_results, "--probe-results"),
        (args.missed_variants, "--missed-variants"),
    ]:
        if not path.exists():
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    # Load data
    print("Loading probe results...", file=sys.stderr)
    probes = load_probe_results(args.probe_results)

    print("Loading missed variants...", file=sys.stderr)
    missed = load_missed_variants(args.missed_variants)

    concordance = None
    if args.concordance_details and args.concordance_details.exists():
        print("Loading concordance details...", file=sys.stderr)
        concordance = load_concordance_details(args.concordance_details)

    window_statuses = None
    if args.log and args.log.exists():
        print("Parsing debug log window statuses...", file=sys.stderr)
        window_statuses = load_debug_log_window_statuses(args.log)

    # Derive lost_at_stage for every row from observed diagnostic facts.
    # This is the sole source of attribution — C++ emits raw data only.
    print("Computing per-row attribution...", file=sys.stderr)
    probes = derive_lost_at(probes)
    probes = compute_depth(probes)

    # Select best record per probe (deepest stage, highest k on tie)
    attribution = select_best_per_probe(probes)

    # Setup output
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "probe_analysis_report.txt"
    console = DualConsole(report_path)

    view = args.view

    if view in ("scorecard", "all"):
        render_scorecard(console, probes, attribution, missed, concordance, window_statuses)

    if view in ("funnel", "all"):
        render_stage_funnel(console, attribution, args.output_dir)

    if view in ("survival", "all"):
        render_survival_analysis(console, probes)

    if view in ("breakdown", "all"):
        render_type_size_breakdown(console, attribution, missed)

    if view in ("genotyper", "all"):
        render_genotyper_forensics(console, attribution)

    if view in ("targets", "all"):
        render_inspection_targets(console, attribution, missed)

    if view in ("deepdive", "all"):
        render_top_loss_deep_dive(console, attribution, missed)

    if view in ("windows", "all"):
        render_window_analysis(console, probes, attribution)

    # Write TSV output files
    write_stage_attribution_tsv(args.output_dir, attribution, missed)
    write_survival_matrix_tsv(args.output_dir, probes)

    print(f"\n  Wrote report to {report_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
