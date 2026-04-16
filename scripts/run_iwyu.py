#!/usr/bin/env python3
"""Run include-what-you-use on Lancet2 source files.

Usage:
    ./scripts/run_iwyu.py                              # check-only (CI gate)
    ./scripts/run_iwyu.py --fix                        # fix includes in-place
    ./scripts/run_iwyu.py --build-dir build            # custom build dir

Prerequisites:
    - Build directory with compile_commands.json
    - pixi environment with include-what-you-use and iwyu_tool.py
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


def ensure_pixi() -> None:
    """Install pixi if it is not already on PATH."""
    if shutil.which("pixi") is not None:
        return
    print("pixi not found, installing...")
    subprocess.run(
        ["bash", "-c", "curl -fsSL https://pixi.sh/install.sh | bash"],
        check=True,
    )


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BUILD_DIR = "cmake-build-release"

# IWYU mapping files to load (resolved at runtime via PIXI_PROJECT_ROOT).
LANCET_MAPPING = "cmake/iwyu/lancet.imp"
BOOST_MAPPING = ".pixi/envs/default/share/include-what-you-use/boost-all.imp"

# IWYU flags tailored for Lancet2's C++20, template-heavy codebase.
IWYU_FLAGS = [
    "--error",
    "--no_fwd_decls",
    "--cxx17ns",
    "--quoted_includes_first",
    "--max_line_length=100",
]


@dataclass
class FileViolation:
    """IWYU findings for a single source file."""

    path: str
    to_add: list[str] = field(default_factory=list)
    to_remove: list[str] = field(default_factory=list)

    @property
    def has_violations(self) -> bool:
        return bool(self.to_add or self.to_remove)


# ──────────────────────────────────────────────────────────────────────────────
# IWYU output parser
# ──────────────────────────────────────────────────────────────────────────────
_RE_SHOULD_ADD = re.compile(r"^(.*) should add these lines:$")
_RE_SHOULD_REMOVE = re.compile(r"^(.*) should remove these lines:$")
_RE_INCLUDE = re.compile(r"^[-#]?\s*#include\s+[\"<]")
_RE_HAS_CORRECT = re.compile(r"^\((.*) has correct #includes/fwd-decls\)$")
_RE_FULL_LIST = re.compile(r"^The full include-list for")
_RE_SEPARATOR = re.compile(r"^---$")


def parse_iwyu_output(output: str, repo_root: Path) -> list[FileViolation]:
    """Parse IWYU output into structured file violations.

    Filters to only files under the repository root (skips third-party deps).
    """
    violations: dict[str, FileViolation] = {}
    current_file: str | None = None
    section: str | None = None  # "add", "remove", or None
    repo_prefix = str(repo_root) + "/"

    for line in output.splitlines():
        line = line.rstrip()

        # Skip lines that aren't about our repository files
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
            rel = filepath[len(repo_prefix) :]
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
            rel = filepath[len(repo_prefix) :]
            current_file = rel
            if rel not in violations:
                violations[rel] = FileViolation(path=rel)
            section = "remove"
            continue

        # Empty line or mapping error — reset section
        if not line.strip() or line.startswith("Cannot open mapping file"):
            section = None
            continue

        # Collect include lines in the current section
        if current_file and section and _RE_INCLUDE.match(line):
            entry = violations[current_file]
            cleaned = line.lstrip("- ").strip()
            if section == "add":
                entry.to_add.append(cleaned)
            elif section == "remove":
                entry.to_remove.append(cleaned)

    return [v for v in violations.values() if v.has_violations]


# ──────────────────────────────────────────────────────────────────────────────
# Summary printing
# ──────────────────────────────────────────────────────────────────────────────
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
    """Print a concise, human-readable summary of IWYU findings."""
    clean = total_files - len(violations)
    total_adds = sum(len(v.to_add) for v in violations)
    total_removes = sum(len(v.to_remove) for v in violations)

    print()
    print(_color("─" * 80, _CYAN))
    print(_color("  IWYU Summary", _BOLD))
    print(_color("─" * 80, _CYAN))
    print(f"  Files analyzed:  {total_files}")
    print(f"  Files clean:     {_color(str(clean), _GREEN)}")
    print(f"  Files with issues: {_color(str(len(violations)), _RED)}")
    print(f"  Missing includes:  {total_adds}")
    print(f"  Unused includes:   {total_removes}")
    print(_color("─" * 80, _CYAN))

    if not violations:
        print(_color("\n  ✓ All files have correct includes.\n", _GREEN))
        return

    # Group by module for readability
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


# ──────────────────────────────────────────────────────────────────────────────
# IWYU command builder
# ──────────────────────────────────────────────────────────────────────────────
def build_iwyu_cmd(build_dir: Path) -> list[str]:
    """Build the iwyu_tool.py command with all mapping files and flags."""
    project_root = os.environ.get("PIXI_PROJECT_ROOT", str(REPO_ROOT))
    lancet_map = str(Path(project_root) / LANCET_MAPPING)
    boost_map = str(Path(project_root) / BOOST_MAPPING)

    cmd: list[str] = [
        "iwyu_tool.py",
        "-p", str(build_dir),
        "--jobs", "0",
        "src/",
    ]

    xiwyu: list[str] = []
    for flag in IWYU_FLAGS:
        xiwyu.extend(["-Xiwyu", flag])
    xiwyu.extend(["-Xiwyu", f"--mapping_file={lancet_map}"])
    if Path(boost_map).exists():
        xiwyu.extend(["-Xiwyu", f"--mapping_file={boost_map}"])

    cmd.extend(["--", *xiwyu])
    return cmd


# ──────────────────────────────────────────────────────────────────────────────
# Fix mode: IWYU → fix_includes.py
# ──────────────────────────────────────────────────────────────────────────────
def run_fix(build_dir: Path) -> int:
    """Run IWYU and apply fixes with fix_includes.py.

    Runs two passes of iwyu → fix_includes to handle headers that only become
    fixable after .cpp files are cleaned in the first pass. Does NOT reformat;
    use clang-format (pixi run fmt-fix) afterwards to reorder include blocks.
    """
    iwyu_cmd = build_iwyu_cmd(build_dir)

    fix_cmd: list[str] = [
        "fix_includes.py",
        "--nocomments",
        "--noreorder",
        "--nosafe_headers",
        "--nokeep_iwyu_namespace_format",
        "--only_re=src/lancet",
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

    print("==> Done. Run 'pixi run fmt-fix' to reorder includes, then 'pixi run iwyu-check' to verify.")
    return 0

# ──────────────────────────────────────────────────────────────────────────────
# Check mode: IWYU analysis only (CI gate)
# ──────────────────────────────────────────────────────────────────────────────
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

    if violations:
        return 1
    return 0


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run include-what-you-use on Lancet2 source files."
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Apply IWYU fixes and reformat (default: check-only)",
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
