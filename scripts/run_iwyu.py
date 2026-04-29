#!/usr/bin/env python3
"""Run include-what-you-use (IWYU) on Lancet2 source files.

This script wraps iwyu_tool.py in two modes:

  Check mode (CI gate):   pixi run iwyu-check   →  ./scripts/run_iwyu.py
  Fix mode (developer):   pixi run iwyu-fix     →  ./scripts/run_iwyu.py --fix

Architecture
------------
IWYU analyzes every translation unit in compile_commands.json and reports
missing or unused #include directives. The output is parsed into structured
FileViolation objects, printed as a color-coded summary, and (in fix mode)
piped to fix_includes.py for automatic patching.

Bracket-Style Problem
---------------------
CMake declares third-party dependencies with SYSTEM (via FetchContent), which
adds their include paths with -isystem. IWYU interprets -isystem paths as
"system" headers and emits angle brackets for them:

    should add:    #include <absl/hash/hash.h>
    should remove: #include "absl/hash/hash.h"

Lancet2 convention reserves angle brackets exclusively for C/C++ standard
library headers. All third-party headers use quotes. This script handles the
mismatch in two ways:

  1. Check mode parser:  cancel_bracket_style_pairs() detects add/remove pairs
     that differ only in bracket style for known third-party prefixes and
     removes them from the violation report. This prevents CI false positives.

  2. Fix mode normalizer: normalize_include_style() runs after fix_includes.py
     and rewrites any angle-bracket third-party includes back to quotes.

Third-Party Prefix List
-----------------------
_THIRD_PARTY_PREFIXES must stay in sync with cmake/dependencies.cmake. When
adding a new FetchContent dependency, add its include prefix here.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BUILD_DIR = "cmake-build-release"

# IWYU mapping files resolved via PIXI_PROJECT_ROOT at runtime.
LANCET_MAPPING = "cmake/iwyu/lancet.imp"

# IWYU analysis flags for Lancet2's C++20 codebase.
IWYU_FLAGS = [
    "--error",
    "--no_fwd_decls",
    "--cxx17ns",
    "--quoted_includes_first",
    "--max_line_length=100",
]

# Third-party include prefixes whose headers IWYU emits with angle brackets
# (because CMake marks them with -isystem). Used by both the check-mode parser
# and the fix-mode normalizer to detect/correct the bracket-style mismatch.
# Derived from cmake/dependencies.cmake — keep in sync when adding new deps.
_THIRD_PARTY_PREFIXES = (
    "absl/",
    "spdlog/",
    "spoa/",
    "CLI/",
    "gperftools/",
    "htslib/",
    "benchmark/",
    "catch_amalgamated",
    "blockingconcurrentqueue",
    "concurrentqueue",
    "mimalloc",
)


# ══════════════════════════════════════════════════════════════════════════════
# Bracket-Style Normalization (fix mode)
#
# After fix_includes.py applies IWYU's recommendations, some third-party
# includes will have angle brackets (e.g. <absl/hash/hash.h>). This pass
# rewrites them to quotes ("absl/hash/hash.h") per project convention.
# ══════════════════════════════════════════════════════════════════════════════

# Matches: #include <absl/...>, #include <spdlog/...>, etc.
_RE_ANGLE_THIRD_PARTY = re.compile(
    r'^(\s*#include\s+)<((?:' + '|'.join(re.escape(p) for p in _THIRD_PARTY_PREFIXES) + r')[^>]*)>',
    re.MULTILINE,
)


def normalize_include_style(repo_root: Path) -> int:
    """Rewrite angle-bracket third-party includes to quoted includes.

    Scans src/, tests/, and benchmarks/ for .cpp and .h files.
    Returns the number of files modified.
    """
    modified = 0
    for src_dir in ("src", "tests", "benchmarks"):
        search_dir = repo_root / src_dir
        if not search_dir.exists():
            continue
        for pattern in ("**/*.cpp", "**/*.h"):
            for filepath in search_dir.glob(pattern):
                text = filepath.read_text()
                new_text = _RE_ANGLE_THIRD_PARTY.sub(r'\1"\2"', text)
                if new_text != text:
                    filepath.write_text(new_text)
                    modified += 1
    return modified


# ══════════════════════════════════════════════════════════════════════════════
# IWYU Output Parser
#
# IWYU emits text blocks like:
#
#   /path/to/file.cpp should add these lines:
#   #include <absl/hash/hash.h>  // for Hash
#
#   /path/to/file.cpp should remove these lines:
#   - #include "absl/hash/hash.h"  // lines 5-5
#
#   The full include-list for /path/to/file.cpp:
#   ...
#   ---
#
# The parser collects "should add" and "should remove" directives into
# FileViolation objects, filters to repository-local files, and cancels
# bracket-style-only pairs before reporting.
# ══════════════════════════════════════════════════════════════════════════════

_RE_SHOULD_ADD = re.compile(r"^(.*) should add these lines:$")
_RE_SHOULD_REMOVE = re.compile(r"^(.*) should remove these lines:$")
_RE_INCLUDE = re.compile(r"^[-#]?\s*#include\s+[\"<]")
_RE_HAS_CORRECT = re.compile(r"^\((.*) has correct #includes/fwd-decls\)$")
_RE_FULL_LIST = re.compile(r"^The full include-list for")
_RE_SEPARATOR = re.compile(r"^---$")

# Extracts the bare header path: #include <absl/foo.h> → absl/foo.h
_RE_HEADER_PATH = re.compile(r'#include\s+["<]([^">]+)[">]')


def _header_path(include_line: str) -> str:
    """Extract bare header path from an #include directive.

    Example: '#include <absl/hash/hash.h>  // for Hash' → 'absl/hash/hash.h'
    """
    match = _RE_HEADER_PATH.search(include_line)
    return match.group(1) if match else include_line


def _is_third_party_header(path: str) -> bool:
    """Return True if the header path belongs to a known third-party library."""
    return any(
        path.startswith(prefix) or path == prefix.rstrip("/")
        for prefix in _THIRD_PARTY_PREFIXES
    )


@dataclass
class FileViolation:
    """IWYU findings for a single source file.

    Attributes:
        path:      Repository-relative file path (e.g. "src/lancet/main.cpp").
        to_add:    Include directives IWYU says are missing.
        to_remove: Include directives IWYU says are unused.
    """

    path: str
    to_add: list[str] = field(default_factory=list)
    to_remove: list[str] = field(default_factory=list)

    def cancel_bracket_style_pairs(self) -> None:
        """Remove add/remove pairs that differ only in bracket style.

        When a file has `#include "absl/foo.h"` (quotes) but IWYU expects
        `#include <absl/foo.h>` (angle brackets), IWYU reports both:
          should add:    #include <absl/foo.h>
          should remove: #include "absl/foo.h"

        The header is present — only the bracket style differs. This method
        detects such pairs for known third-party headers and removes both
        entries, leaving only genuine missing/unused violations.
        """
        add_paths = {_header_path(line) for line in self.to_add}
        remove_paths = {_header_path(line) for line in self.to_remove}
        bracket_only = {
            p for p in (add_paths & remove_paths) if _is_third_party_header(p)
        }
        if not bracket_only:
            return
        self.to_add = [
            line for line in self.to_add if _header_path(line) not in bracket_only
        ]
        self.to_remove = [
            line for line in self.to_remove if _header_path(line) not in bracket_only
        ]

    @property
    def has_violations(self) -> bool:
        return bool(self.to_add or self.to_remove)


def parse_iwyu_output(output: str, repo_root: Path) -> list[FileViolation]:
    """Parse IWYU text output into structured FileViolation objects.

    Only files under repo_root are included (third-party dep output is
    ignored). Bracket-style-only pairs are cancelled before returning.
    """
    violations: dict[str, FileViolation] = {}
    current_file: str | None = None
    section: str | None = None  # "add", "remove", or None
    repo_prefix = str(repo_root) + "/"

    for line in output.splitlines():
        line = line.rstrip()

        # Skip correct-include confirmations and full-list/separator blocks
        if _RE_HAS_CORRECT.match(line):
            continue
        if _RE_FULL_LIST.match(line) or _RE_SEPARATOR.match(line):
            section = None
            continue

        # "should add" section header
        match = _RE_SHOULD_ADD.match(line)
        if match:
            filepath = match.group(1)
            if not filepath.startswith(repo_prefix):
                current_file = None
                section = None
                continue
            rel = filepath[len(repo_prefix):]
            current_file = rel
            if rel not in violations:
                violations[rel] = FileViolation(path=rel)
            section = "add"
            continue

        # "should remove" section header
        match = _RE_SHOULD_REMOVE.match(line)
        if match:
            filepath = match.group(1)
            if not filepath.startswith(repo_prefix):
                current_file = None
                section = None
                continue
            rel = filepath[len(repo_prefix):]
            current_file = rel
            if rel not in violations:
                violations[rel] = FileViolation(path=rel)
            section = "remove"
            continue

        # Empty line or mapping error resets the section
        if not line.strip() or line.startswith("Cannot open mapping file"):
            section = None
            continue

        # Collect #include lines within the current add/remove section
        if current_file and section and _RE_INCLUDE.match(line):
            entry = violations[current_file]
            cleaned = line.lstrip("- ").strip()
            if section == "add":
                entry.to_add.append(cleaned)
            elif section == "remove":
                entry.to_remove.append(cleaned)

    # Cancel bracket-style-only pairs before reporting
    for entry in violations.values():
        entry.cancel_bracket_style_pairs()

    return [v for v in violations.values() if v.has_violations]


# ══════════════════════════════════════════════════════════════════════════════
# Summary Display
# ══════════════════════════════════════════════════════════════════════════════

_YELLOW = "\033[33m"
_RED = "\033[31m"
_GREEN = "\033[32m"
_CYAN = "\033[36m"
_BOLD = "\033[1m"
_RESET = "\033[0m"


def _color(text: str, code: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"{code}{text}{_RESET}"


def print_summary(violations: list[FileViolation], total_files: int) -> None:
    """Print a concise, color-coded summary of IWYU findings."""
    clean = total_files - len(violations)
    total_adds = sum(len(v.to_add) for v in violations)
    total_removes = sum(len(v.to_remove) for v in violations)

    print()
    print(_color("─" * 80, _CYAN))
    print(_color("  IWYU Summary", _BOLD))
    print(_color("─" * 80, _CYAN))
    print(f"  Files analyzed:    {total_files}")
    print(f"  Files clean:       {_color(str(clean), _GREEN)}")
    print(f"  Files with issues: {_color(str(len(violations)), _RED)}")
    print(f"  Missing includes:  {total_adds}")
    print(f"  Unused includes:   {total_removes}")
    print(_color("─" * 80, _CYAN))

    if not violations:
        print(_color("\n  ✓ All files have correct includes.\n", _GREEN))
        return

    violations.sort(key=lambda v: v.path)
    for v in violations:
        print(f"\n  {_color(v.path, _BOLD)}")
        for inc in v.to_add:
            print(f"    {_color('+', _GREEN)} {inc}")
        for inc in v.to_remove:
            print(f"    {_color('-', _RED)} {inc}")

    print()
    print(
        _color(
            f"  ✗ {len(violations)} file(s) have include issues. "
            f"Fix the includes above to pass this check.",
            _YELLOW,
        )
    )
    print()


# ══════════════════════════════════════════════════════════════════════════════
# IWYU Command Builder
# ══════════════════════════════════════════════════════════════════════════════

def build_iwyu_cmd(build_dir: Path) -> list[str]:
    """Build the iwyu_tool.py command with mapping files and analysis flags."""
    project_root = os.environ.get("PIXI_PROJECT_ROOT", str(REPO_ROOT))
    lancet_map = str(Path(project_root) / LANCET_MAPPING)

    cmd: list[str] = [
        "iwyu_tool.py",
        "-p", str(build_dir),
        "--jobs", "0",
        "src/", "tests/", "benchmarks/",
    ]

    xiwyu: list[str] = []
    for flag in IWYU_FLAGS:
        xiwyu.extend(["-Xiwyu", flag])
    xiwyu.extend(["-Xiwyu", f"--mapping_file={lancet_map}"])

    cmd.extend(["--", *xiwyu])
    return cmd


# ══════════════════════════════════════════════════════════════════════════════
# Fix Mode: IWYU → fix_includes.py → normalize bracket style
#
# Pipeline:
#   1. Run IWYU analysis (up to 2 passes — headers may only become fixable
#      after their .cpp files are cleaned in pass 1).
#   2. Pipe IWYU output to fix_includes.py to apply add/remove edits.
#   3. Run normalize_include_style() to convert any angle-bracket third-party
#      includes back to quoted style.
#   4. Developer runs `pixi run fmt-fix` afterwards to reorder include blocks.
# ══════════════════════════════════════════════════════════════════════════════

def run_fix(build_dir: Path) -> int:
    """Run IWYU and apply fixes. Returns 0 on success, 1 on error."""
    iwyu_cmd = build_iwyu_cmd(build_dir)

    fix_cmd: list[str] = [
        "fix_includes.py",
        "--nocomments",
        "--noreorder",
        "--nosafe_headers",
        "--nokeep_iwyu_namespace_format",
        "--only_re=src/lancet|tests/|benchmarks/",
    ]

    for pass_num in (1, 2):
        print(f"==> IWYU fix pass {pass_num}: analyzing includes...")
        iwyu_result = subprocess.run(
            iwyu_cmd, cwd=REPO_ROOT,
            capture_output=True, text=True,
        )
        iwyu_output = iwyu_result.stdout + iwyu_result.stderr

        violations = parse_iwyu_output(iwyu_output, REPO_ROOT)
        if not violations:
            print(f"    No issues found in pass {pass_num}.")
            break

        total_adds = sum(len(v.to_add) for v in violations)
        total_removes = sum(len(v.to_remove) for v in violations)
        print(
            f"    Found {len(violations)} file(s) with issues "
            f"(+{total_adds} missing, -{total_removes} unused)"
        )

        print(f"==> IWYU fix pass {pass_num}: applying fixes...")
        fix_result = subprocess.run(
            fix_cmd, cwd=REPO_ROOT,
            input=iwyu_output, capture_output=True, text=True,
        )
        if fix_result.returncode != 0:
            print(f"ERROR: fix_includes.py failed:\n{fix_result.stderr}", file=sys.stderr)
            return 1

        edited = fix_result.stdout.count("Fixing #includes in")
        skipped = fix_result.stdout.count("no contentful changes")
        print(f"    Edited {edited} file(s), skipped {skipped} file(s)")

    # Post-fix: rewrite <third-party/...> → "third-party/..." per convention
    normalized = normalize_include_style(REPO_ROOT)
    if normalized > 0:
        print(f"==> Normalized include style in {normalized} file(s)")

    print("==> Done. Run 'pixi run fmt-fix' to reorder includes, "
          "then 'pixi run iwyu-check' to verify.")
    return 0


# ══════════════════════════════════════════════════════════════════════════════
# Check Mode: IWYU analysis only (CI gate)
#
# Runs IWYU, parses the output, cancels bracket-style false positives, and
# prints a summary. Returns 0 if all files are clean, 1 if violations exist.
# ══════════════════════════════════════════════════════════════════════════════

def run_check(build_dir: Path) -> int:
    """Run IWYU in check-only mode. Returns 1 if any violations are found."""
    cmd = build_iwyu_cmd(build_dir)

    print("==> Running IWYU on Lancet sources...")
    print(f"    Build dir: {build_dir.relative_to(REPO_ROOT)}")
    print(f"    Mapping:   {LANCET_MAPPING}")
    print()

    result = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    output = result.stdout + result.stderr

    correct_count = len(_RE_HAS_CORRECT.findall(output))
    violations = parse_iwyu_output(output, REPO_ROOT)
    total_files = correct_count + len(violations)

    print_summary(violations, total_files)

    return 1 if violations else 0


# ══════════════════════════════════════════════════════════════════════════════
# Entry Point
# ══════════════════════════════════════════════════════════════════════════════

def ensure_pixi() -> None:
    """Install pixi if it is not already on PATH."""
    if shutil.which("pixi") is not None:
        return
    print("pixi not found, installing...")
    subprocess.run(
        ["bash", "-c", "curl -fsSL https://pixi.sh/install.sh | bash"],
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run include-what-you-use on Lancet2 source files.",
        epilog="Check mode (default) is the CI gate. Fix mode applies changes in-place.",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Apply IWYU fixes in-place (default: check-only)",
    )
    parser.add_argument(
        "--build-dir",
        default=DEFAULT_BUILD_DIR,
        help=f"Path to CMake build directory (default: {DEFAULT_BUILD_DIR})",
    )
    args = parser.parse_args()

    build_dir = REPO_ROOT / args.build_dir
    compile_db = build_dir / "compile_commands.json"

    if not compile_db.exists():
        print(
            f"ERROR: {compile_db.relative_to(REPO_ROOT)} not found.\n"
            f"Run: pixi run build",
            file=sys.stderr,
        )
        return 1

    ensure_pixi()

    if args.fix:
        return run_fix(build_dir)
    return run_check(build_dir)


if __name__ == "__main__":
    raise SystemExit(main())
