#!/usr/bin/env python3
"""Combine all Lancet2 documentation into a single Markdown file.

Dynamically reads the navigation structure from mkdocs.yml so the output
stays in sync even as pages are added, removed, or reordered.

Usage:
    python3 scripts/export_docs.py                          # writes to stdout
    python3 scripts/export_docs.py > Lancet2_Docs.md        # redirect to file
    python3 scripts/export_docs.py -o Lancet2_Docs.md       # write to file directly
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = ROOT / "docs"
MKDOCS_YML = ROOT / "mkdocs.yml"


def get_version_info() -> dict[str, str]:
    """Extract git version metadata matching the docker build script format.

    Format: VERSION_TAG-BRANCH-COMMIT[-dirty]
    Example: v0.9.7-main-a3b4c5d6e7
    """
    def git(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, stderr=subprocess.DEVNULL
            ).decode().strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "unknown"

    version_tag = git("describe", "--abbrev=0", "--tags")
    # Keep only major.minor.patch if tag has extra segments
    version_tag = ".".join(version_tag.split(".")[:3])

    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    commit = git("rev-parse", "--short=10", "--verify", "HEAD")

    try:
        subprocess.check_call(
            ["git", "diff", "--quiet"], cwd=ROOT, stderr=subprocess.DEVNULL
        )
        dirty = ""
    except (subprocess.CalledProcessError, FileNotFoundError):
        dirty = "-dirty"

    build_tag = f"{version_tag}-{branch}-{commit}{dirty}"
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    return {"build_tag": build_tag, "timestamp": timestamp, "commit": commit}


# ---------------------------------------------------------------------------
# Nav entry: (title, filepath | None, depth)
#   filepath=None means a section header with no content file.
# ---------------------------------------------------------------------------
NavEntry = tuple[str, str | None, int]


def parse_nav(nav_list: list, depth: int = 1) -> list[NavEntry]:
    """Recursively parse the mkdocs.yml nav list into flat (title, path, depth) entries."""
    entries: list[NavEntry] = []
    for item in nav_list:
        if isinstance(item, str):
            # Bare string: "file.md" with no title — use filename as title
            entries.append((Path(item).stem.replace("_", " ").title(), item, depth))
        elif isinstance(item, dict):
            for title, value in item.items():
                if isinstance(value, str):
                    # "Title: path.md"
                    entries.append((title, value, depth))
                elif isinstance(value, list):
                    # "Section Title:" followed by children
                    entries.append((title, None, depth))
                    entries.extend(parse_nav(value, depth + 1))
    return entries


def strip_frontmatter(content: str) -> str:
    """Remove YAML front matter (---...---) from markdown content."""
    return re.sub(r"^---\s*\n.*?\n---\s*\n", "", content, count=1, flags=re.DOTALL)


def strip_image_lines(content: str) -> str:
    """Remove image embed lines since they won't render outside the site."""
    return re.sub(r"^!\[.*?\]\(.*?\)\s*$", "", content, flags=re.MULTILINE)


def bump_headings(content: str, level: int) -> str:
    """Bump all markdown headings by `level` levels so page content nests under the TOC."""
    def replacer(m: re.Match) -> str:
        return "#" * (len(m.group(1)) + level) + m.group(2)
    return re.sub(r"^(#{1,6})(.*)", replacer, content, flags=re.MULTILINE)


def make_anchor(title: str) -> str:
    """Convert a title to a GitHub-style markdown anchor."""
    anchor = title.lower()
    anchor = re.sub(r"[^a-z0-9 -]", "", anchor)
    anchor = anchor.strip().replace(" ", "-")
    anchor = re.sub(r"-+", "-", anchor)
    return anchor


def build_document(entries: list[NavEntry]) -> str:
    """Assemble the combined markdown document from nav entries."""
    parts: list[str] = []
    version = get_version_info()

    # Title page
    parts.append("# Lancet2 — Complete Documentation\n")
    parts.append(f"> **Version:** `{version['build_tag']}`\n")
    parts.append(f"> **Generated:** {version['timestamp']}\n")
    parts.append("> **Source:** <https://github.com/nygenome/Lancet2>\n")

    # Table of contents
    parts.append("\n---\n\n## Table of Contents\n")
    for title, filepath, depth in entries:
        indent = "  " * (depth - 1)
        if filepath is None:
            parts.append(f"{indent}- **{title}**")
        else:
            parts.append(f"{indent}- [{title}](#{make_anchor(title)})")
    parts.append("\n---\n")

    # Content
    for title, filepath, depth in entries:
        if filepath is None:
            parts.append(f"\n---\n\n# {title}\n")
            continue

        src = DOCS_DIR / filepath
        if not src.exists():
            parts.append(f"\n## {title}\n\n> *File not found: {filepath}*\n")
            continue

        raw = src.read_text(encoding="utf-8")
        content = strip_frontmatter(raw)
        content = strip_image_lines(content)
        # Strip the first H1 from the file — we use our own section header
        content = re.sub(r"^#\s+.*\n", "", content, count=1)
        content = content.strip()

        if depth == 1:
            parts.append(f"\n---\n\n# {title}\n")
        else:
            parts.append(f"\n---\n\n## {title}\n")

        bumped = bump_headings(content, depth)
        parts.append(bumped)
        parts.append("")

    return "\n".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export Lancet2 docs to a single Markdown file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file path (default: stdout)")
    args = parser.parse_args()

    if not MKDOCS_YML.exists():
        sys.exit(f"mkdocs.yml not found at {MKDOCS_YML}")

    with open(MKDOCS_YML, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    nav_list = config.get("nav")
    if not nav_list:
        sys.exit("No 'nav' key found in mkdocs.yml")

    entries = parse_nav(nav_list)
    document = build_document(entries)

    if args.output:
        output_path = Path(args.output)
        output_path.write_text(document, encoding="utf-8")
        size_kb = output_path.stat().st_size / 1024
        print(f"✓ Written {output_path.name} ({size_kb:.0f} KB, {len(entries)} sections)", file=sys.stderr)
    else:
        sys.stdout.write(document)


if __name__ == "__main__":
    main()
