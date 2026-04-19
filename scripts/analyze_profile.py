#!/usr/bin/env python3
"""Analyze gperftools CPU profile data via pprof for Lancet2.

Usage:
    python3 scripts/analyze_profile.py PROFILE_FILE
    python3 scripts/analyze_profile.py PROFILE_FILE --binary path/to/Lancet2
    python3 scripts/analyze_profile.py PROFILE_FILE --view components
    python3 scripts/analyze_profile.py PROFILE_FILE --view top --top 50
    python3 scripts/analyze_profile.py PROFILE_FILE --diff-base OLD_PROFILE
    python3 scripts/analyze_profile.py PROFILE_FILE --html report.html
    python3 scripts/analyze_profile.py PROFILE_FILE --list "BuildGraph"
    python3 scripts/analyze_profile.py PROFILE_FILE --save-summary v2.10-pre
    python3 scripts/analyze_profile.py --history

Prerequisites:
    - pprof binary on PATH (pixi run -e profiling ensure-pprof)
    - Binary for symbolization (auto-detected from profile directory)
    - pixi profiling environment for rich/jinja2 (pixi run -e profiling ...)
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

REPO_ROOT = Path(__file__).resolve().parent.parent

# Single console instance; record=True enables save_html() for free.
console = Console(record=True)

VIEWS = ("overview", "top", "modules", "components", "tree", "hotpaths", "lines", "all")

VIEW_HELP = """\
Views:
  overview     Total time, sampling rate, top-level stats
  top          Flat and cumulative hottest functions
  modules      Group by Lancet2 layer and external deps
  components   High-level component attribution
  tree         Caller/callee tree for hottest functions
  hotpaths     Focused drill-down on top Lancet2 functions
  lines        Source-line attribution (NOT included in 'all')
  all          All of the above except 'lines' (default)
"""

HISTORY_FILE = REPO_ROOT / "profiling" / "history.jsonl"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Data Model
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


@dataclass(frozen=True, slots=True)
class FunctionEntry:
    """One row from pprof -text output."""

    flat_sec: float
    flat_pct: float
    sum_pct: float
    cum_sec: float
    cum_pct: float
    name: str


@dataclass(slots=True)
class ProfileMeta:
    """Metadata parsed from the pprof header block."""

    filename: str = ""
    profile_type: str = ""
    total_sec: float = 0.0
    dropped_nodes: int = 0
    drop_threshold_sec: float = 0.0
    shown_nodes: int = 0
    total_nodes: int = 0


@dataclass(frozen=True, slots=True)
class ProfileData:
    """Complete parsed profile: metadata + flat-sorted + cumulative-sorted entries."""

    meta: ProfileMeta
    by_flat: list[FunctionEntry]
    by_cum: list[FunctionEntry]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Function Classification
#
# Every profiled function is mapped to a *module* (fine-grained), then to a
# *component* (coarse-grained).  Rules are evaluated in order — first match wins.
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# (prefix, module) — order matters: more specific prefixes first.
_MODULE_RULES: list[tuple[str, str]] = [
    # Lancet2 layers
    ("lancet::base::", "lancet/base"),
    ("lancet::hts::", "lancet/hts"),
    ("lancet::cbdg::", "lancet/cbdg"),
    ("lancet::caller::", "lancet/caller"),
    ("lancet::core::", "lancet/core"),
    ("lancet::cli::", "lancet/cli"),
    # Third-party libraries
    ("spoa::", "spoa"),
    ("ksw_", "minimap2/ksw2"),
    ("mm_", "minimap2"),
    ("mg_", "minimap2"),
    ("radix_sort_128x", "minimap2"),
    ("idx_read", "htslib"),
    ("bgzf_", "htslib"),
    ("hts_", "htslib"),
    ("sam_", "htslib"),
    ("fai_", "htslib"),
    ("hread", "htslib"),
    ("fd_read", "htslib"),
    ("deflate_", "htslib/zlib"),
    ("inflate_", "htslib/zlib"),
    ("absl::", "abseil"),
    # System and runtime
    ("__madvise", "sys/madvise"),
    ("__read", "sys/io"),
    ("__write", "sys/io"),
    ("open64", "sys/io"),
    ("close", "sys/io"),
    ("__memcpy", "sys/memcpy"),
    ("__memmove", "sys/memmove"),
    ("__memset", "sys/memset"),
    ("clock_gettime", "sys/clock"),
    ("__libc_", "sys/libc"),
    ("mi_", "mimalloc"),
    ("_mi_", "mimalloc"),
    ("operator new", "allocator"),
    ("operator delete", "allocator"),
    ("std::", "libstdc++"),
    ("clone", "sys/thread"),
    ("start_thread", "sys/thread"),
    ("execute_native_thread", "sys/thread"),
]

# module → component (coarse grouping for the component view).
_COMPONENT_MAP: dict[str, str] = {
    "lancet/base": "Lancet2 (graph assembly)",
    "lancet/cbdg": "Lancet2 (graph assembly)",
    "lancet/caller": "Lancet2 (genotyping)",
    "lancet/core": "Lancet2 (pipeline)",
    "lancet/cli": "Lancet2 (pipeline)",
    "lancet/hts": "Lancet2 (HTS I/O)",
    "minimap2": "minimap2",
    "minimap2/ksw2": "minimap2",
    "spoa": "SPOA (MSA)",
    "htslib": "HTSlib",
    "htslib/zlib": "HTSlib",
    "abseil": "Abseil",
    "mimalloc": "Memory management",
    "allocator": "Memory management",
    "libstdc++": "C++ stdlib",
}

# All sys/* modules map to "System".
_SYSTEM_PREFIX = "sys/"

# Stable display order and per-component color (rich markup style names).
_COMPONENT_STYLES: dict[str, str] = {
    "Lancet2 (graph assembly)": "magenta",
    "Lancet2 (genotyping)": "blue",
    "Lancet2 (pipeline)": "cyan",
    "Lancet2 (HTS I/O)": "cyan",
    "minimap2": "red",
    "SPOA (MSA)": "yellow",
    "HTSlib": "green",
    "Memory management": "dim",
    "Abseil": "dim",
    "System": "dim",
    "C++ stdlib": "dim",
    "Other": "white",
}


def classify_module(name: str) -> str:
    """Map a function name to a fine-grained module via prefix matching."""
    for prefix, module in _MODULE_RULES:
        if name.startswith(prefix):
            return module
    return "other"


def classify_component(module: str) -> str:
    """Map a module to a coarse component group."""
    if module.startswith(_SYSTEM_PREFIX):
        return "System"
    return _COMPONENT_MAP.get(module, "Other")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# pprof Driver
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def find_pprof() -> str:
    """Locate the pprof binary, preferring the pixi-managed installation."""
    path = shutil.which("pprof")
    if path:
        return path

    for env in ("profiling", "default"):
        pixi_path = REPO_ROOT / ".pixi" / "envs" / env / "bin" / "pprof"
        if pixi_path.exists():
            return str(pixi_path)

    console.print("[red]ERROR:[/] pprof not found on PATH.\nRun: pixi run -e profiling ensure-pprof")
    sys.exit(1)


def invoke_pprof(
    pprof_bin: str,
    profile: str,
    binary: str | None,
    flags: list[str],
) -> str:
    """Run pprof with the given flags and return its stdout."""
    cmd = [pprof_bin, *flags]
    if binary:
        cmd.append(binary)
    cmd.append(profile)

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        console.print(f"[red]ERROR:[/] pprof failed (exit {result.returncode}): "
                       f"{result.stderr.strip()}")
        sys.exit(1)
    return result.stdout


def auto_detect_binary(profile: str) -> str | None:
    """Search common locations for the Lancet2 binary."""
    candidates = [
        Path(profile).resolve().parent / "Lancet2",
        REPO_ROOT / "cmake-build-relwithdebinfo" / "Lancet2",
        REPO_ROOT / "cmake-build-release" / "Lancet2",
        REPO_ROOT / "cmake-build-debug" / "Lancet2",
    ]
    for path in candidates:
        if path.exists():
            return str(path)
    return None


_BUILD_DIRS = [
    "cmake-build-relwithdebinfo",
    "cmake-build-release",
    "cmake-build-debug",
]


def resolve_profile(raw: str) -> Path:
    """Resolve a profile path, searching common build directories if needed.

    pixi always sets cwd to the repo root, so a bare filename like
    ``Lancet.cpu_profile.*.bin`` typed from inside a build directory won't
    resolve as-is.  This function tries the literal path first, then searches
    each ``cmake-build-*`` directory under REPO_ROOT.
    """
    path = Path(raw)
    if path.exists():
        return path.resolve()

    for build_dir in _BUILD_DIRS:
        candidate = REPO_ROOT / build_dir / raw
        if candidate.exists():
            return candidate.resolve()

    return path  # Return as-is; caller will report the error


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Parsing pprof Text Output
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_RE_FUNCTION = re.compile(
    r"^\s*"
    r"(?P<flat>[\d.]+)[a-z]*\s+"
    r"(?P<flat_pct>[\d.]+)%\s+"
    r"(?P<sum_pct>[\d.]+)%\s+"
    r"(?P<cum>[\d.]+)[a-z]*\s+"
    r"(?P<cum_pct>[\d.]+)%\s+"
    r"(?P<name>.+)$"
)
_RE_FILE = re.compile(r"^File:\s+(.+)$")
_RE_TYPE = re.compile(r"^Type:\s+(.+)$")
_RE_TOTAL = re.compile(r"of\s+([\d.]+)s\s+total")
_RE_DROPPED = re.compile(r"Dropped\s+(\d+)\s+nodes\s+\(cum\s+<=\s+([\d.]+)s\)")
_RE_SHOWING = re.compile(r"Showing\s+top\s+(\d+)\s+nodes\s+out\s+of\s+(\d+)")

# Diff mode: pprof output can have negative values like "-12.34s".
_RE_FUNCTION_DIFF = re.compile(
    r"^\s*"
    r"(?P<flat>-?[\d.]+)[a-z]*\s+"
    r"(?P<flat_pct>-?[\d.]+)%\s+"
    r"(?P<sum_pct>-?[\d.]+)%\s+"
    r"(?P<cum>-?[\d.]+)[a-z]*\s+"
    r"(?P<cum_pct>-?[\d.]+)%\s+"
    r"(?P<name>.+)$"
)


def _parse_one_text_output(raw: str, *, allow_negative: bool = False) -> tuple[ProfileMeta, list[FunctionEntry]]:
    """Parse a single pprof -text invocation into metadata and entries."""
    meta = ProfileMeta()
    entries: list[FunctionEntry] = []
    regex = _RE_FUNCTION_DIFF if allow_negative else _RE_FUNCTION

    for line in raw.splitlines():
        stripped = line.strip()

        if m := _RE_FILE.match(stripped):
            meta.filename = m.group(1)
        elif m := _RE_TYPE.match(stripped):
            meta.profile_type = m.group(1)
        elif m := _RE_TOTAL.search(stripped):
            meta.total_sec = float(m.group(1))
        elif m := _RE_DROPPED.search(stripped):
            meta.dropped_nodes = int(m.group(1))
            meta.drop_threshold_sec = float(m.group(2))
        elif m := _RE_SHOWING.search(stripped):
            meta.shown_nodes = int(m.group(1))
            meta.total_nodes = int(m.group(2))
        elif m := regex.match(stripped):
            entries.append(FunctionEntry(
                flat_sec=float(m["flat"]),
                flat_pct=float(m["flat_pct"]),
                sum_pct=float(m["sum_pct"]),
                cum_sec=float(m["cum"]),
                cum_pct=float(m["cum_pct"]),
                name=m["name"].strip(),
            ))

    # When nodecount >= actual nodes, pprof omits the "Showing top" line.
    if meta.shown_nodes == 0 and entries:
        meta.shown_nodes = len(entries)
        meta.total_nodes = len(entries) + meta.dropped_nodes

    return meta, entries


def load_profile(
    pprof_bin: str, profile: str, binary: str | None, nodecount: int
) -> ProfileData:
    """Fetch flat-sorted and cum-sorted views from pprof and parse them."""
    common_flags = ["-text", f"-nodecount={nodecount}"]

    raw_flat = invoke_pprof(pprof_bin, profile, binary, [*common_flags, "-flat"])
    raw_cum = invoke_pprof(pprof_bin, profile, binary, [*common_flags, "-cum"])

    meta, by_flat = _parse_one_text_output(raw_flat)
    _, by_cum = _parse_one_text_output(raw_cum)

    if meta.total_sec == 0:
        console.print("[red]ERROR:[/] Could not parse profile. Is the file a valid gperftools profile?")
        sys.exit(1)

    return ProfileData(meta=meta, by_flat=by_flat, by_cum=by_cum)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Aggregation Helpers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


@dataclass(slots=True)
class GroupStats:
    """Accumulated self-time and constituent functions for one group."""

    flat_sec: float = 0.0
    functions: list[FunctionEntry] | None = None

    def add(self, entry: FunctionEntry) -> None:
        self.flat_sec += entry.flat_sec
        if self.functions is None:
            self.functions = []
        self.functions.append(entry)

    def top_function(self) -> FunctionEntry:
        assert self.functions
        return max(self.functions, key=lambda e: e.flat_sec)


def aggregate_by(entries: list[FunctionEntry], key_fn) -> dict[str, GroupStats]:
    """Group entries by key_fn(entry.name) and accumulate self-time."""
    groups: dict[str, GroupStats] = {}
    for entry in entries:
        key = key_fn(entry.name)
        if key not in groups:
            groups[key] = GroupStats()
        groups[key].add(entry)
    return groups


def sorted_groups(groups: dict[str, GroupStats]) -> list[tuple[str, GroupStats]]:
    """Return groups sorted by descending self-time."""
    return sorted(groups.items(), key=lambda kv: -kv[1].flat_sec)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Rich Bar Helper
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def make_bar(fraction: float, width: int = 30, style: str = "white") -> Text:
    """Render a horizontal bar using Unicode block characters as rich Text."""
    fraction = max(0.0, min(1.0, fraction))
    filled = int(fraction * width)
    half = "▌" if (fraction * width - filled) > 0.5 else ""
    bar_str = f"{'█' * filled}{half}".ljust(width, "░")
    return Text(bar_str, style=style)


def severity_style(pct: float, *, thresholds: tuple[float, float] = (5.0, 2.0)) -> str:
    """Pick a rich style based on percentage severity thresholds."""
    high, med = thresholds
    if pct >= high:
        return "red"
    if pct >= med:
        return "yellow"
    if pct >= 0.5:
        return "white"
    return "dim"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Report Views
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def render_overview(data: ProfileData) -> None:
    """Profile summary: total time, function counts, quick stats."""
    meta = data.meta
    total = meta.total_sec
    accounted = sum(e.flat_sec for e in data.by_flat)
    top5_pct = sum(e.flat_sec for e in data.by_flat[:5]) / total * 100

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Key", style="bold")
    table.add_column("Value")
    table.add_row("Binary", meta.filename)
    table.add_row("Profile type", meta.profile_type)
    table.add_row("Total CPU time", f"[bold]{total:.2f}s[/] ({total / 60:.1f} min)")
    table.add_row("Functions shown", f"{meta.shown_nodes} of {meta.total_nodes}")
    table.add_row("Dropped nodes", f"{meta.dropped_nodes} (cum ≤ {meta.drop_threshold_sec:.2f}s)")
    table.add_row("Accounted flat", f"{accounted:.2f}s ({accounted / total * 100:.1f}% of total)")
    table.add_row("Top 5 self-time", f"[yellow]{top5_pct:.1f}%[/] of total CPU")

    console.print(Panel(table, title="[bold]Profile Overview[/]", border_style="cyan"))


def _render_function_table(
    entries: list[FunctionEntry],
    total: float,
    limit: int,
    *,
    primary: str,
    thresholds: tuple[float, float],
    title: str,
) -> None:
    """Render a ranked function table (shared logic for flat and cum views)."""
    is_flat = primary == "flat"

    table = Table(title=title, title_style="bold", border_style="dim", expand=True)
    table.add_column("Rank", justify="right", style="dim", width=5)
    if is_flat:
        table.add_column("Self (s)", justify="right", width=10)
        table.add_column("Self%", justify="right", width=7)
        table.add_column("Cum (s)", justify="right", style="dim", width=10)
        table.add_column("Cum%", justify="right", style="dim", width=7)
    else:
        table.add_column("Cum (s)", justify="right", width=10)
        table.add_column("Cum%", justify="right", width=7)
        table.add_column("Self (s)", justify="right", style="dim", width=10)
        table.add_column("Self%", justify="right", style="dim", width=7)
    table.add_column("Bar", width=16, no_wrap=True)
    table.add_column("Function", no_wrap=True)

    for rank, entry in enumerate(entries[:limit], 1):
        pct = (entry.flat_sec if is_flat else entry.cum_sec) / total * 100
        style = severity_style(pct, thresholds=thresholds)
        bar = make_bar(pct / 100, width=15, style=style)

        if is_flat:
            table.add_row(
                str(rank),
                Text(f"{entry.flat_sec:8.2f}s", style=style),
                Text(f"{pct:5.1f}%", style=style),
                f"{entry.cum_sec:8.2f}s",
                f"{entry.cum_pct:5.1f}%",
                bar,
                entry.name,
            )
        else:
            table.add_row(
                str(rank),
                Text(f"{entry.cum_sec:8.2f}s", style=style),
                Text(f"{pct:5.1f}%", style=style),
                f"{entry.flat_sec:8.2f}s",
                f"{entry.flat_pct:5.1f}%",
                bar,
                entry.name,
            )

    console.print(table)


def render_top(data: ProfileData, limit: int) -> None:
    """Top functions by self-time and cumulative time."""
    total = data.meta.total_sec
    n_flat = min(limit, len(data.by_flat))
    n_cum = min(limit, len(data.by_cum))

    _render_function_table(
        data.by_flat, total, limit, primary="flat", thresholds=(5.0, 2.0),
        title=f"Top {n_flat} Functions by Self-Time (flat)",
    )
    console.print()
    _render_function_table(
        data.by_cum, total, limit, primary="cum", thresholds=(20.0, 5.0),
        title=f"Top {n_cum} Functions by Cumulative Time (cum)",
    )


def render_modules(data: ProfileData) -> None:
    """Functions grouped by fine-grained module with per-Lancet2 drill-down."""
    total = data.meta.total_sec
    groups = aggregate_by(data.by_flat, classify_module)

    table = Table(title="Module Breakdown (by self-time)", title_style="bold",
                  border_style="dim", expand=True)
    table.add_column("Module", width=30)
    table.add_column("Self (s)", justify="right", width=10)
    table.add_column("Self%", justify="right", width=7)
    table.add_column("Bar", width=21, no_wrap=True)
    table.add_column("Top Function", no_wrap=True, style="dim")

    for mod, stats in sorted_groups(groups):
        if stats.flat_sec < 0.01:
            continue
        pct = stats.flat_sec / total * 100
        style = severity_style(pct)
        bar = make_bar(pct / 100, width=20, style=style)
        top_fn = stats.top_function().name

        table.add_row(
            Text(mod, style=style), Text(f"{stats.flat_sec:8.2f}s", style=style),
            Text(f"{pct:5.1f}%", style=style), bar, top_fn,
        )

    console.print(table)

    # Drill down into each significant Lancet2 module
    for mod, stats in sorted_groups(groups):
        if not mod.startswith("lancet/") or stats.flat_sec < 0.5:
            continue

        detail = Table(title=f"{mod} — function detail", title_style="bold",
                       border_style="dim")
        detail.add_column("Self (s)", justify="right", width=10)
        detail.add_column("Self%", justify="right", width=7)
        detail.add_column("Cum (s)", justify="right", width=10)
        detail.add_column("Function")

        funcs = sorted(stats.functions or [], key=lambda e: -e.flat_sec)
        for fn in funcs[:10]:
            pct = fn.flat_sec / total * 100
            detail.add_row(f"{fn.flat_sec:8.2f}s", f"{pct:5.1f}%",
                           f"{fn.cum_sec:8.2f}s", fn.name)

        console.print(detail)


def render_components(data: ProfileData) -> None:
    """High-level component attribution with internal-vs-external summary."""
    total = data.meta.total_sec

    # Two-level aggregation: entries → modules → components
    mod_groups = aggregate_by(data.by_flat, classify_module)
    comp_flat: dict[str, float] = {}
    comp_subs: dict[str, dict[str, float]] = {}

    for mod, stats in mod_groups.items():
        comp = classify_component(mod)
        comp_flat[comp] = comp_flat.get(comp, 0.0) + stats.flat_sec
        comp_subs.setdefault(comp, {})[mod] = stats.flat_sec

    accounted = sum(comp_flat.values())
    unaccounted = total - accounted

    table = Table(title="Component Attribution (self-time)", title_style="bold",
                  border_style="dim", expand=True,
                  caption=f"Total: {total:.2f}s │ Accounted: {accounted:.2f}s │ "
                          f"Unaccounted: {unaccounted:.2f}s ({unaccounted / total * 100:.1f}%)")
    table.add_column("Component", width=35)
    table.add_column("Self (s)", justify="right", width=10)
    table.add_column("Self%", justify="right", width=7)
    table.add_column("Bar", width=31, no_wrap=True)

    for comp, sec in sorted(comp_flat.items(), key=lambda kv: -kv[1]):
        pct = sec / total * 100
        if pct < 0.01:
            continue
        style = _COMPONENT_STYLES.get(comp, "white")
        bar = make_bar(pct / 100, width=30, style=style)

        table.add_row(Text(comp, style=style), Text(f"{sec:8.2f}s", style=style),
                      Text(f"{pct:5.1f}%", style=style), bar)

        for mod, mod_sec in sorted(comp_subs[comp].items(), key=lambda kv: -kv[1]):
            if mod_sec < 1.0:
                continue
            table.add_row(f"  └─ {mod}", f"{mod_sec:8.2f}s",
                          f"{mod_sec / total * 100:5.1f}%", "")

    console.print(table)

    # Internal vs External summary
    lancet = sum(v for k, v in comp_flat.items() if k.startswith("Lancet2"))
    mm2 = comp_flat.get("minimap2", 0.0)
    htslib = comp_flat.get("HTSlib", 0.0)
    spoa = comp_flat.get("SPOA (MSA)", 0.0)
    other = total - lancet - mm2 - htslib - spoa

    summary = Table(title="Internal vs External", title_style="bold", border_style="dim")
    summary.add_column("Category", width=35)
    summary.add_column("Time (s)", justify="right", width=10)
    summary.add_column("Pct", justify="right", width=7)
    summary.add_column("Bar", width=26, no_wrap=True)

    for label, sec, style in [
        ("Lancet2 own code", lancet, "cyan"),
        ("minimap2", mm2, "red"),
        ("HTSlib", htslib, "green"),
        ("SPOA", spoa, "yellow"),
        ("Other (alloc, system, etc)", other, "dim"),
    ]:
        pct = sec / total * 100
        bar = make_bar(pct / 100, width=25, style=style)
        summary.add_row(Text(label, style=style), Text(f"{sec:8.2f}s", style=style),
                        Text(f"{pct:5.1f}%", style=style), bar)

    console.print(summary)


def render_tree(pprof_bin: str, profile: str, binary: str | None, focus: str | None) -> None:
    """Render pprof's caller/callee tree, optionally focused on a regex."""
    flags = ["-tree", "-nodecount=25"]
    if focus:
        flags.append(f"-focus={focus}")

    raw = invoke_pprof(pprof_bin, profile, binary, flags)

    tree = Tree("[bold]Call Tree (caller → callee)[/]")
    for line in raw.splitlines():
        if line.startswith("File:") or line.startswith("Type:"):
            continue
        stripped = line.strip()
        if any(kw in stripped for kw in ("accounting for", "Dropped", "Showing top")):
            tree.add(f"[dim]{stripped}[/]")
        elif line.startswith("---"):
            tree.add("[cyan]" + "─" * 70 + "[/]")
        else:
            tree.add(line)

    console.print(Panel(tree, border_style="cyan"))


def render_hotpaths(
    pprof_bin: str, profile: str, binary: str | None, data: ProfileData
) -> None:
    """Drill into the top Lancet2 functions by cumulative time."""
    targets = [
        e for e in data.by_cum
        if e.name.startswith("lancet::") and e.cum_pct >= 5.0
    ][:5]

    if not targets:
        console.print(Panel("[dim]No Lancet2 functions with ≥5% cumulative time found.[/]",
                            title="[bold]Hot Path Drill-Down[/]", border_style="cyan"))
        return

    for entry in targets:
        flags = ["-tree", "-nodecount=10", f"-focus={re.escape(entry.name)}"]
        raw = invoke_pprof(pprof_bin, profile, binary, flags)

        tree = Tree(f"[bold]{entry.name}[/] — {entry.cum_sec:.1f}s ({entry.cum_pct:.1f}% cum)")
        for line in raw.splitlines():
            if line.startswith("File:") or line.startswith("Type:"):
                continue
            stripped = line.strip()
            if "accounting for" in stripped or "Dropped" in stripped:
                continue
            if "Showing top" in stripped:
                tree.add(f"[dim]{stripped}[/]")
            elif line.startswith("---"):
                tree.add("[dim]" + "─" * 60 + "[/]")
            else:
                tree.add(line)

        console.print(Panel(tree, title="[bold]Hot Path Drill-Down[/]", border_style="cyan"))


def render_lines(pprof_bin: str, profile: str, binary: str | None, nodecount: int) -> None:
    """Source-line level attribution (not included in 'all' view)."""
    raw = invoke_pprof(pprof_bin, profile, binary,
                       ["-text", "-lines", f"-nodecount={nodecount}", "-flat"])

    _, entries = _parse_one_text_output(raw)
    if not entries:
        console.print("[dim]No line-level data available.[/]")
        return

    tree = Tree("[bold]Source-Line Attribution[/]")
    # Group by function (everything before the last space-separated filename:line)
    current_fn: str | None = None
    fn_branch: Tree | None = None

    for entry in entries[:100]:
        # pprof -lines output format: "func file:line"
        parts = entry.name.rsplit(" ", 1)
        fn_name = parts[0] if len(parts) == 2 else entry.name
        location = parts[1] if len(parts) == 2 else ""

        if fn_name != current_fn:
            current_fn = fn_name
            fn_branch = tree.add(f"[bold]{fn_name}[/]")

        style = severity_style(entry.flat_pct)
        assert fn_branch is not None
        fn_branch.add(f"[{style}]{entry.flat_sec:8.2f}s ({entry.flat_pct:5.1f}%)[/] {location}")

    console.print(Panel(tree, border_style="cyan"))


def render_list(pprof_bin: str, profile: str, binary: str | None, pattern: str) -> None:
    """Annotated source view via pprof -list."""
    raw = invoke_pprof(pprof_bin, profile, binary, [f"-list={pattern}"])

    if not raw.strip():
        console.print(f"[yellow]No source found matching pattern:[/] {pattern}")
        return

    try:
        from rich.syntax import Syntax
        # pprof -list output is annotated C++ source, render with syntax highlighting
        console.print(Panel(
            Syntax(raw, "cpp", theme="monokai", line_numbers=False),
            title=f"[bold]Annotated Source: {pattern}[/]",
            border_style="cyan",
        ))
    except ImportError:
        console.print(Panel(raw, title=f"[bold]Annotated Source: {pattern}[/]",
                            border_style="cyan"))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Diff Mode
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def load_diff(
    pprof_bin: str, base_profile: str, new_profile: str,
    binary: str | None, nodecount: int,
) -> tuple[ProfileData, ProfileData, ProfileData]:
    """Load base, new, and diff profiles for comparison."""
    base = load_profile(pprof_bin, base_profile, binary, nodecount)
    new = load_profile(pprof_bin, new_profile, binary, nodecount)

    # Diff profile: values can be negative (improvement)
    diff_flags = ["-text", f"-nodecount={nodecount}", "-flat",
                  f"-diff_base={base_profile}"]
    raw_diff = invoke_pprof(pprof_bin, new_profile, binary, diff_flags)
    diff_meta, diff_entries = _parse_one_text_output(raw_diff, allow_negative=True)
    diff_data = ProfileData(meta=diff_meta, by_flat=diff_entries, by_cum=[])

    return base, new, diff_data


def render_diff(base: ProfileData, new: ProfileData, diff: ProfileData) -> None:
    """Render the diff report between two profiles."""
    base_total = base.meta.total_sec
    new_total = new.meta.total_sec
    delta = new_total - base_total
    delta_pct = delta / base_total * 100 if base_total > 0 else 0

    verdict_style = "red bold" if delta > 0 else "green bold"
    verdict = "REGRESSION ▲" if delta > 0 else "IMPROVEMENT ▼"

    # Overview
    overview = Table(show_header=False, box=None, padding=(0, 2))
    overview.add_column("Key", style="bold")
    overview.add_column("Value")
    overview.add_row("Base total", f"{base_total:.2f}s ({base_total / 60:.1f} min)")
    overview.add_row("New total", f"{new_total:.2f}s ({new_total / 60:.1f} min)")
    overview.add_row("Delta", Text(f"{delta:+.2f}s ({delta_pct:+.1f}%)", style=verdict_style))
    overview.add_row("Verdict", Text(verdict, style=verdict_style))

    console.print(Panel(overview, title="[bold]Profile Diff Overview[/]", border_style="cyan"))

    # Build lookup for base functions
    base_lookup: dict[str, FunctionEntry] = {e.name: e for e in base.by_flat}
    new_lookup: dict[str, FunctionEntry] = {e.name: e for e in new.by_flat}
    all_names = sorted(set(base_lookup) | set(new_lookup),
                       key=lambda n: abs(new_lookup.get(n, FunctionEntry(0, 0, 0, 0, 0, n)).flat_sec
                                         - base_lookup.get(n, FunctionEntry(0, 0, 0, 0, 0, n)).flat_sec),
                       reverse=True)

    # Top changes table
    table = Table(title="Top Changes by Self-Time Delta", title_style="bold",
                  border_style="dim", expand=True)
    table.add_column("Function", no_wrap=True)
    table.add_column("Base (s)", justify="right", width=10)
    table.add_column("New (s)", justify="right", width=10)
    table.add_column("Δ (s)", justify="right", width=10)
    table.add_column("Δ%", justify="right", width=8)
    table.add_column("Status", width=14)

    regressions: list[str] = []
    for name in all_names[:30]:
        b = base_lookup.get(name)
        n = new_lookup.get(name)
        b_sec = b.flat_sec if b else 0.0
        n_sec = n.flat_sec if n else 0.0
        d_sec = n_sec - b_sec
        d_pct = d_sec / b_sec * 100 if b_sec > 0 else float("inf") if d_sec > 0 else 0

        if abs(d_sec) < 0.01:
            continue

        if d_sec > 0:
            style = "red"
            status = "🔴 regressed"
            if d_sec > 1.0 and (b_sec == 0 or d_pct > 10):
                regressions.append(f"{name}: {d_sec:+.2f}s ({d_pct:+.1f}%)")
        else:
            style = "green"
            status = "🟢 improved"

        d_pct_str = f"{d_pct:+.1f}%" if abs(d_pct) < 1e6 else "new"
        table.add_row(
            name,
            f"{b_sec:8.2f}s",
            f"{n_sec:8.2f}s",
            Text(f"{d_sec:+8.2f}s", style=style),
            Text(d_pct_str, style=style),
            Text(status, style=style),
        )

    console.print(table)

    # Regression warnings
    if regressions:
        warning_text = "\n".join(f"  • {r}" for r in regressions)
        console.print(Panel(
            f"[red bold]Functions with significant regressions (>10% AND >1s):[/]\n{warning_text}",
            title="[red bold]⚠ Regressions Detected[/]",
            border_style="red",
        ))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Git Version Metadata
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def _git(*args: str, check: bool = True) -> str:
    """Run a git command and return its stripped stdout."""
    try:
        return subprocess.check_output(
            ["git", *args], cwd=REPO_ROOT, stderr=subprocess.DEVNULL
        ).decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def get_git_version() -> dict:
    """Compute git version metadata matching build_push_image.sh's generate_tag.

    Format: VERSION_TAG-BRANCH-COMMIT[-dirty]
    Example: v2.9.0-main-a3b4c5d6e7-dirty
    """
    raw_tag = _git("describe", "--abbrev=0", "--tags")
    # Keep only major.minor.patch if tag has extra segments
    version_tag = ".".join(raw_tag.split(".")[:3])
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    commit = _git("rev-parse", "--short=10", "--verify", "HEAD")

    try:
        subprocess.check_call(
            ["git", "diff", "--quiet"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL
        )
        dirty = False
    except (subprocess.CalledProcessError, FileNotFoundError):
        dirty = True

    dirty_suffix = "-dirty" if dirty else ""
    full = f"{version_tag}-{branch}-{commit}{dirty_suffix}"

    return {
        "full": full,
        "version_tag": version_tag,
        "branch": branch,
        "commit": commit,
        "dirty": dirty,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Profile History
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def save_history_entry(tag: str, data: ProfileData, profile_path: str,
                       binary_path: str | None) -> None:
    """Append a JSON summary to the history file.

    Captures the complete pre-symbolized function list so that future
    cross-binary diffs (--diff-tag) work by function name without
    needing the original binary.  The binary is only required once —
    at save time.
    """
    total = data.meta.total_sec

    # Compute component and module breakdowns
    mod_groups = aggregate_by(data.by_flat, classify_module)
    comp_flat: dict[str, float] = {}
    for mod, stats in mod_groups.items():
        comp = classify_component(mod)
        comp_flat[comp] = comp_flat.get(comp, 0.0) + stats.flat_sec

    # Module breakdown (Lancet2 modules only)
    modules = {}
    for mod, stats in sorted_groups(mod_groups):
        if stats.flat_sec >= 0.1:
            modules[mod] = {"flat_sec": round(stats.flat_sec, 2),
                            "flat_pct": round(stats.flat_sec / total * 100, 1)}

    # Component breakdown
    components = {}
    for comp, sec in sorted(comp_flat.items(), key=lambda kv: -kv[1]):
        if sec >= 0.1:
            components[comp] = {"flat_sec": round(sec, 2),
                                "flat_pct": round(sec / total * 100, 1)}

    # Internal vs external
    lancet = sum(v for k, v in comp_flat.items() if k.startswith("Lancet2"))
    mm2 = comp_flat.get("minimap2", 0.0)
    htslib = comp_flat.get("HTSlib", 0.0)
    spoa = comp_flat.get("SPOA (MSA)", 0.0)
    other_sec = total - lancet - mm2 - htslib - spoa

    # Top 10 functions (kept for backward compatibility)
    top_10 = [
        {"name": e.name, "flat_sec": round(e.flat_sec, 2), "flat_pct": round(e.flat_pct, 1),
         "cum_sec": round(e.cum_sec, 2), "cum_pct": round(e.cum_pct, 1)}
        for e in data.by_flat[:10]
    ]

    # ALL profiled functions — the pre-symbolized snapshot used by --diff-tag.
    # Stored once at save time; the binary is never needed again for this entry.
    functions = [
        {"name": e.name, "flat_sec": round(e.flat_sec, 2), "flat_pct": round(e.flat_pct, 1),
         "cum_sec": round(e.cum_sec, 2), "cum_pct": round(e.cum_pct, 1)}
        for e in data.by_flat
    ]

    entry = {
        "tag": tag,
        "version": get_git_version(),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "profile_file": profile_path,
        "binary": binary_path or "",
        "total_sec": round(total, 2),
        "top_10": top_10,
        "functions": functions,
        "components": components,
        "modules": modules,
        "internal_vs_external": {
            "lancet2_own": {"flat_sec": round(lancet, 2), "flat_pct": round(lancet / total * 100, 1)},
            "minimap2": {"flat_sec": round(mm2, 2), "flat_pct": round(mm2 / total * 100, 1)},
            "htslib": {"flat_sec": round(htslib, 2), "flat_pct": round(htslib / total * 100, 1)},
            "spoa": {"flat_sec": round(spoa, 2), "flat_pct": round(spoa / total * 100, 1)},
            "other": {"flat_sec": round(other_sec, 2), "flat_pct": round(other_sec / total * 100, 1)},
        },
    }

    HISTORY_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(HISTORY_FILE, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, separators=(",", ":")) + "\n")

    console.print(f"[green]✓[/] Saved profile summary with tag [bold]{tag}[/] "
                  f"({len(functions)} functions) to {HISTORY_FILE.relative_to(REPO_ROOT)}")


def _load_history() -> list[dict]:
    """Load all entries from the profile history file."""
    if not HISTORY_FILE.exists():
        return []
    entries = []
    with open(HISTORY_FILE, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    return entries


def _find_history_entry(tag: str) -> dict | None:
    """Find a history entry by tag name (case-insensitive)."""
    for entry in _load_history():
        if entry.get("tag", "").lower() == tag.lower():
            return entry
    return None


def _compute_function_deltas(
    base_funcs: list[dict], new_funcs: list[dict],
) -> list[dict]:
    """Compute per-function deltas between two pre-symbolized function lists.

    Each input item has keys: name, flat_sec, flat_pct, cum_sec, cum_pct.
    Returns a list of delta dicts sorted by absolute flat_sec delta descending.
    """
    base_lookup = {f["name"]: f for f in base_funcs}
    new_lookup = {f["name"]: f for f in new_funcs}
    all_names = set(base_lookup) | set(new_lookup)

    deltas = []
    for name in all_names:
        b = base_lookup.get(name)
        n = new_lookup.get(name)
        b_flat = b["flat_sec"] if b else 0.0
        n_flat = n["flat_sec"] if n else 0.0
        d_flat = n_flat - b_flat
        d_pct = d_flat / b_flat * 100 if b_flat > 0 else (
            float("inf") if d_flat > 0 else (float("-inf") if d_flat < 0 else 0)
        )
        b_cum = b["cum_sec"] if b else 0.0
        n_cum = n["cum_sec"] if n else 0.0

        deltas.append({
            "name": name,
            "base_flat": b_flat, "new_flat": n_flat,
            "delta_flat": d_flat, "delta_pct": d_pct,
            "base_cum": b_cum, "new_cum": n_cum,
            "status": "new" if not b else "removed" if not n else (
                "regressed" if d_flat > 0 else "improved" if d_flat < 0 else "unchanged"
            ),
        })

    deltas.sort(key=lambda d: abs(d["delta_flat"]), reverse=True)
    return deltas


def _compute_component_deltas(base_entry: dict, new_entry: dict) -> list[dict]:
    """Compute per-component deltas between two history entries."""
    base_comps = base_entry.get("components", {})
    new_comps = new_entry.get("components", {})
    all_comp_names = sorted(set(base_comps) | set(new_comps))

    deltas = []
    for comp in all_comp_names:
        b = base_comps.get(comp, {})
        n = new_comps.get(comp, {})
        b_sec = b.get("flat_sec", 0.0)
        n_sec = n.get("flat_sec", 0.0)
        d_sec = n_sec - b_sec
        d_pct = d_sec / b_sec * 100 if b_sec > 0 else (
            float("inf") if d_sec > 0 else 0
        )
        deltas.append({
            "name": comp, "base_sec": b_sec, "new_sec": n_sec,
            "delta_sec": d_sec, "delta_pct": d_pct,
        })
    deltas.sort(key=lambda d: abs(d["delta_sec"]), reverse=True)
    return deltas


def render_tag_diff(base_tag: str, new_tag: str) -> None:
    """Render a cross-binary profile diff using pre-symbolized history data."""
    base_entry = _find_history_entry(base_tag)
    new_entry = _find_history_entry(new_tag)

    if not base_entry:
        console.print(f"[red]ERROR:[/] No history entry found for tag '{base_tag}'")
        console.print("Available tags:", ", ".join(
            e["tag"] for e in _load_history()) or "(none)")
        return
    if not new_entry:
        console.print(f"[red]ERROR:[/] No history entry found for tag '{new_tag}'")
        console.print("Available tags:", ", ".join(
            e["tag"] for e in _load_history()) or "(none)")
        return

    base_funcs = base_entry.get("functions", base_entry.get("top_10", []))
    new_funcs = new_entry.get("functions", new_entry.get("top_10", []))

    if "functions" not in base_entry:
        console.print(f"[yellow]WARNING:[/] Tag '{base_tag}' has no full function data "
                      f"(saved before --diff-tag support). Using top_10 only. "
                      f"Re-run --save-summary with the original binary for full data.")
    if "functions" not in new_entry:
        console.print(f"[yellow]WARNING:[/] Tag '{new_tag}' has no full function data. "
                      f"Re-run --save-summary with the original binary for full data.")

    # ── Overview ──────────────────────────────────────────────────────────────
    base_total = base_entry["total_sec"]
    new_total = new_entry["total_sec"]
    delta = new_total - base_total
    delta_pct = delta / base_total * 100 if base_total > 0 else 0

    verdict_style = "red bold" if delta > 0 else "green bold"
    verdict = "REGRESSION ▲" if delta > 0 else "IMPROVEMENT ▼"

    overview = Table(show_header=False, box=None, padding=(0, 2))
    overview.add_column("Key", style="bold")
    overview.add_column("Value")
    overview.add_row("Base", f"{base_tag} ({base_entry.get('version', {}).get('full', '?')})")
    overview.add_row("New", f"{new_tag} ({new_entry.get('version', {}).get('full', '?')})")
    overview.add_row("Base total", f"{base_total:.2f}s ({base_total / 60:.1f} min)")
    overview.add_row("New total", f"{new_total:.2f}s ({new_total / 60:.1f} min)")
    overview.add_row("Delta", Text(f"{delta:+.2f}s ({delta_pct:+.1f}%)", style=verdict_style))
    overview.add_row("Verdict", Text(verdict, style=verdict_style))
    console.print(Panel(overview, title="[bold]Cross-Binary Profile Diff[/]", border_style="cyan"))

    # ── Component deltas ─────────────────────────────────────────────────────
    comp_deltas = _compute_component_deltas(base_entry, new_entry)
    comp_table = Table(title="Component-Level Changes", title_style="bold",
                       border_style="dim", expand=True)
    comp_table.add_column("Component", no_wrap=True)
    comp_table.add_column("Base (s)", justify="right", width=10)
    comp_table.add_column("New (s)", justify="right", width=10)
    comp_table.add_column("Δ (s)", justify="right", width=10)
    comp_table.add_column("Δ%", justify="right", width=8)

    for cd in comp_deltas:
        if abs(cd["delta_sec"]) < 0.1:
            continue
        style = "red" if cd["delta_sec"] > 0 else "green"
        d_pct_str = f"{cd['delta_pct']:+.1f}%" if abs(cd["delta_pct"]) < 1e6 else "new"
        comp_table.add_row(
            cd["name"], f"{cd['base_sec']:.2f}", f"{cd['new_sec']:.2f}",
            Text(f"{cd['delta_sec']:+.2f}", style=style),
            Text(d_pct_str, style=style),
        )
    console.print(comp_table)

    # ── Function deltas ──────────────────────────────────────────────────────
    func_deltas = _compute_function_deltas(base_funcs, new_funcs)

    func_table = Table(title="Top Function Changes by Self-Time Delta", title_style="bold",
                       border_style="dim", expand=True)
    func_table.add_column("Function", no_wrap=True)
    func_table.add_column("Base (s)", justify="right", width=10)
    func_table.add_column("New (s)", justify="right", width=10)
    func_table.add_column("Δ (s)", justify="right", width=10)
    func_table.add_column("Δ%", justify="right", width=8)
    func_table.add_column("Status", width=14)

    STATUS_ICONS = {
        "regressed": ("red", "🔴 regressed"),
        "improved": ("green", "🟢 improved"),
        "new": ("yellow", "🆕 new"),
        "removed": ("dim", "⛔ removed"),
        "unchanged": ("dim", "— unchanged"),
    }

    regressions: list[str] = []
    improvements: list[str] = []
    for fd in func_deltas[:40]:
        if abs(fd["delta_flat"]) < 0.01:
            continue
        style, status_text = STATUS_ICONS.get(fd["status"], ("dim", fd["status"]))
        d_pct_str = (
            f"{fd['delta_pct']:+.1f}%" if abs(fd["delta_pct"]) < 1e6
            else "new" if fd["status"] == "new" else "removed"
        )
        func_table.add_row(
            fd["name"], f"{fd['base_flat']:.2f}", f"{fd['new_flat']:.2f}",
            Text(f"{fd['delta_flat']:+.2f}", style=style),
            Text(d_pct_str, style=style),
            Text(status_text, style=style),
        )
        if fd["delta_flat"] > 1.0 and fd["status"] in ("regressed", "new"):
            regressions.append(f"{fd['name']}: {fd['delta_flat']:+.2f}s ({d_pct_str})")
        elif fd["delta_flat"] < -1.0 and fd["status"] in ("improved", "removed"):
            improvements.append(f"{fd['name']}: {fd['delta_flat']:+.2f}s ({d_pct_str})")

    console.print(func_table)

    # ── Summary panels ───────────────────────────────────────────────────────
    if improvements:
        imp_text = "\n".join(f"  • {r}" for r in improvements[:15])
        console.print(Panel(
            f"[green bold]Top improvements (>1s self-time reduction):[/]\n{imp_text}",
            title="[green bold]✓ Improvements[/]", border_style="green",
        ))
    if regressions:
        reg_text = "\n".join(f"  • {r}" for r in regressions[:15])
        console.print(Panel(
            f"[red bold]Regressions (>1s self-time increase):[/]\n{reg_text}",
            title="[red bold]⚠ Regressions Detected[/]", border_style="red",
        ))


def render_history() -> None:
    """Show a trend table from the profile history file."""
    if not HISTORY_FILE.exists():
        console.print(f"[yellow]No history file found at {HISTORY_FILE.relative_to(REPO_ROOT)}[/]")
        console.print("Run with --save-summary <tag> after analyzing a profile.")
        return

    entries = _load_history()
    if not entries:
        console.print("[yellow]History file is empty.[/]")
        return

    table = Table(title="Profile History — Performance Trends", title_style="bold",
                  border_style="dim", expand=True)
    table.add_column("Tag", style="bold")
    table.add_column("Version", style="dim")
    table.add_column("Date", width=12)
    table.add_column("Total", justify="right", width=8)
    table.add_column("minimap2", justify="right", width=9)
    table.add_column("graph", justify="right", width=9)
    table.add_column("genotype", justify="right", width=9)
    table.add_column("SPOA", justify="right", width=9)

    for entry in entries:
        version = entry.get("version", {}).get("full", "?")
        ts = entry.get("timestamp", "")[:10]
        total_min = f"{entry['total_sec'] / 60:.1f}m"
        comps = entry.get("components", {})

        mm2_pct = comps.get("minimap2", {}).get("flat_pct", 0)
        graph_pct = comps.get("Lancet2 (graph assembly)", {}).get("flat_pct", 0)
        geno_pct = comps.get("Lancet2 (genotyping)", {}).get("flat_pct", 0)
        spoa_pct = comps.get("SPOA (MSA)", {}).get("flat_pct", 0)

        table.add_row(
            entry["tag"], version, ts, total_min,
            f"{mm2_pct:.1f}%", f"{graph_pct:.1f}%", f"{geno_pct:.1f}%", f"{spoa_pct:.1f}%",
        )

    console.print(table)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# HTML Report Generation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def _fix_pprof_svg(raw_svg: str) -> str:
    """Post-process pprof's SVG output for correct inline HTML embedding.

    The SVG generated by graphviz (via pprof -svg) has two problems:
    1. width="100%" height="100%" with no viewBox — the browser cannot compute
       intrinsic height, so it collapses to CSS default (150px).
    2. XML preamble and DOCTYPE that conflict when embedded inside HTML5.

    Fix: strip the preamble, remove the viewport scale(0.5,0.5) transform,
    extract the graph bounding polygon for exact dimensions, and set a viewBox
    at origin (0,0) — because graph0's translate() maps the polygon into
    positive coordinate space.
    """
    # Strip XML preamble and DOCTYPE for inline embedding
    svg = re.sub(r'<\?xml[^?]*\?>\s*', '', raw_svg)
    svg = re.sub(r'<!DOCTYPE[^>]*>\s*', '', svg)
    svg = re.sub(r'<!--[^-]*-->\s*', '', svg, count=1)  # strip graphviz comment

    # Remove the viewport scale(0.5,0.5) — we use the full coordinate space
    # in the viewBox and let the browser/CSS handle display scaling.
    svg = re.sub(
        r'transform="scale\([^)]+\)\s*translate\([^)]+\)"',
        'transform="translate(0,0)"',
        svg,
    )

    # Extract graph bounding box from the first polygon (graphviz background rect).
    # The polygon lives in graph0's local coordinate space (y-up, negative values).
    # graph0's transform="translate(4, 1276.25)" maps everything to positive space.
    # So the viewBox origin is (0, 0) and the size is the polygon's width × height.
    poly_match = re.search(r'<polygon[^>]*points="([^"]+)"', svg)

    if poly_match:
        coords = [tuple(map(float, p.split(',')))
                  for p in poly_match.group(1).strip().split()]
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        # Width and height of the graph in its native coordinate space
        vb_w = int(max(xs) - min(xs)) + 80  # padding for text labels
        vb_h = int(max(ys) - min(ys)) + 80

        # viewBox at (0,0) since graph0's translate maps content to positive space
        svg = re.sub(
            r'<svg\s+width="100%"\s+height="100%"',
            f'<svg viewBox="0 0 {vb_w} {vb_h}" '
            f'width="100%" '
            f'style="max-width:{vb_w}px; height:auto; aspect-ratio:{vb_w}/{vb_h}"',
            svg,
        )

    return svg


def generate_html_report(
    data: ProfileData,
    pprof_bin: str,
    profile_str: str,
    binary: str | None,
    output_path: str,
    *,
    diff_data: tuple[ProfileData, ProfileData, ProfileData] | None = None,
) -> None:
    """Generate a self-contained HTML report using the Jinja2 template."""
    from jinja2 import Environment, FileSystemLoader

    template_dir = Path(__file__).parent
    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=True)
    template = env.get_template("profile_report.html.j2")

    total = data.meta.total_sec

    # Compute component breakdown for the report
    mod_groups = aggregate_by(data.by_flat, classify_module)
    comp_flat: dict[str, float] = {}
    comp_subs: dict[str, dict[str, float]] = {}
    for mod, stats in mod_groups.items():
        comp = classify_component(mod)
        comp_flat[comp] = comp_flat.get(comp, 0.0) + stats.flat_sec
        comp_subs.setdefault(comp, {})[mod] = stats.flat_sec

    components = []
    for comp, sec in sorted(comp_flat.items(), key=lambda kv: -kv[1]):
        pct = sec / total * 100
        if pct < 0.01:
            continue
        subs = [{"name": m, "sec": round(s, 2), "pct": round(s / total * 100, 1)}
                for m, s in sorted(comp_subs[comp].items(), key=lambda kv: -kv[1])
                if s >= 1.0]
        components.append({"name": comp, "sec": round(sec, 2), "pct": round(pct, 1),
                           "style": _COMPONENT_STYLES.get(comp, "white"), "subs": subs})

    # Generate SVG callgraph
    svg = ""
    try:
        raw_svg = invoke_pprof(pprof_bin, profile_str, binary,
                               ["-svg", "-nodecount=20"])
        svg = _fix_pprof_svg(raw_svg)
    except SystemExit:
        pass  # SVG generation failed, skip it

    # Build functions list
    functions = [
        {"name": e.name, "flat_sec": round(e.flat_sec, 2), "flat_pct": round(e.flat_pct, 1),
         "cum_sec": round(e.cum_sec, 2), "cum_pct": round(e.cum_pct, 1)}
        for e in data.by_flat[:50]
    ]

    report = {
        "meta": {
            "binary": data.meta.filename,
            "profile_type": data.meta.profile_type,
            "total_sec": round(total, 2),
            "total_min": round(total / 60, 1),
            "shown_nodes": data.meta.shown_nodes,
            "total_nodes": data.meta.total_nodes,
        },
        "functions": functions,
        "components": components,
        "callgraph_svg": svg,
        "version": get_git_version(),
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
    }

    html = template.render(report=report)
    Path(output_path).write_text(html, encoding="utf-8")
    console.print(f"[green]✓[/] HTML report saved to [bold]{output_path}[/]")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CLI Entry Point
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze gperftools CPU profile data via pprof.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=VIEW_HELP,
    )
    parser.add_argument("profile", nargs="?", default=None,
                        help="Path to gperftools CPU profile binary file")
    parser.add_argument("--binary", default=None,
                        help="Path to the profiled binary (auto-detected if omitted)")
    parser.add_argument("--view", default="all", choices=list(VIEWS),
                        help="Which report view to display (default: all)")
    parser.add_argument("--top", type=int, default=30, dest="top_n",
                        help="Number of top functions to show (default: 30)")
    parser.add_argument("--focus", default=None,
                        help="Regex to focus the tree view on specific functions")
    parser.add_argument("--nodecount", type=int, default=200,
                        help="Max nodes for pprof to consider (default: 200)")
    parser.add_argument("--diff-base", default=None, dest="diff_base",
                        help="Base profile for same-binary diff comparison")
    parser.add_argument("--diff-tag", nargs=2, default=None, dest="diff_tags",
                        metavar=("BASE_TAG", "NEW_TAG"),
                        help="Cross-binary diff via pre-symbolized history entries")
    parser.add_argument("--html", default=None,
                        help="Output path for self-contained HTML report")
    parser.add_argument("--list", default=None, dest="list_pattern",
                        help="Regex for annotated source view (pprof -list)")
    parser.add_argument("--save-summary", default=None, dest="save_summary",
                        help="Tag name to save profile summary to history")
    parser.add_argument("--history", action="store_true",
                        help="Show profile history trend table")
    return parser


def main() -> int:
    args = build_parser().parse_args()

    # ── History mode (no profile needed) ─────────────────────────────────────
    if args.history:
        render_history()
        return 0

    # ── Cross-binary diff mode (no profile needed) ───────────────────────────
    if args.diff_tags:
        base_tag, new_tag = args.diff_tags
        console.print(Panel("[bold]Lancet2 Cross-Binary Profile Diff[/]", border_style="cyan"))
        render_tag_diff(base_tag, new_tag)
        return 0

    # ── Validate inputs ──────────────────────────────────────────────────────
    if args.profile is None:
        console.print("[red]ERROR:[/] Profile file is required (unless using --history)")
        return 1

    profile = resolve_profile(args.profile)
    if not profile.exists():
        searched = [f"  • {args.profile} (literal)"]
        for bd in _BUILD_DIRS:
            searched.append(f"  • {REPO_ROOT / bd / args.profile}")
        console.print(
            f"[red]ERROR:[/] Profile file not found: {args.profile}\n"
            f"Searched:\n" + "\n".join(searched)
        )
        return 1
    profile_str = str(profile)

    binary = args.binary or auto_detect_binary(profile_str)
    if binary and not Path(binary).exists():
        console.print(f"[yellow]WARNING:[/] Binary not found at {binary}, "
                      f"proceeding without symbolization")
        binary = None

    pprof_bin = find_pprof()

    # ── Annotated source view (standalone mode) ──────────────────────────────
    if args.list_pattern:
        render_list(pprof_bin, profile_str, binary, args.list_pattern)
        return 0

    # ── Diff mode ────────────────────────────────────────────────────────────
    if args.diff_base:
        base_profile = resolve_profile(args.diff_base)
        if not base_profile.exists():
            console.print(f"[red]ERROR:[/] Base profile not found: {args.diff_base}")
            return 1

        console.print(Panel("[bold]Lancet2 CPU Profile Diff Report[/]", border_style="cyan"))

        base, new, diff = load_diff(pprof_bin, str(base_profile), profile_str,
                                    binary, args.nodecount)
        render_diff(base, new, diff)

        if args.html:
            generate_html_report(new, pprof_bin, profile_str, binary, args.html,
                                 diff_data=(base, new, diff))
        return 0

    # ── Load and parse ───────────────────────────────────────────────────────
    data = load_profile(pprof_bin, profile_str, binary, args.nodecount)

    # ── Render ───────────────────────────────────────────────────────────────
    console.print()
    console.print(Panel("[bold]Lancet2 CPU Profile Analysis Report[/]", border_style="cyan"))

    show = args.view
    show_all = show == "all"

    if show_all or show == "overview":
        render_overview(data)
    if show_all or show == "top":
        render_top(data, args.top_n)
    if show_all or show == "modules":
        render_modules(data)
    if show_all or show == "components":
        render_components(data)
    if show_all or show == "tree":
        render_tree(pprof_bin, profile_str, binary, args.focus)
    if show_all or show == "hotpaths":
        render_hotpaths(pprof_bin, profile_str, binary, data)
    # 'lines' is NOT included in 'all' — must be explicitly requested
    if show == "lines":
        render_lines(pprof_bin, profile_str, binary, args.nodecount)

    # ── Save summary ─────────────────────────────────────────────────────────
    if args.save_summary:
        save_history_entry(args.save_summary, data, profile_str, binary)

    # ── HTML report ──────────────────────────────────────────────────────────
    if args.html:
        generate_html_report(data, pprof_bin, profile_str, binary, args.html)

    # ── Footer ───────────────────────────────────────────────────────────────
    console.print()
    footer = Table(show_header=False, box=None, padding=(0, 1))
    footer.add_column("Key", style="dim")
    footer.add_column("Value", style="dim")
    footer.add_row("Profile", profile_str)
    if binary:
        footer.add_row("Binary", binary)
    footer.add_row("pprof", pprof_bin)
    console.print(Panel(footer, border_style="cyan"))
    console.print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
