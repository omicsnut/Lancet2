#!/usr/bin/env bash
# Lancet2 PostToolUse observability hook: log_post_tool_use
#
# Silent file logging on PostToolUse for Edit/Write/MultiEdit/Read.
# Records, for each tool call, which path-scoped rules' globs would
# have matched the file_path — i.e., which rules' bodies were
# (probably) loaded into context as a result of this call. The
# matching is a substring + glob check that mirrors how Claude Code's
# memory system scopes rule loading.
#
# Why this matters: the InstructionsLoaded observability story has two
# halves —
#
#   log_session_start.sh records what was AVAILABLE (full rule set).
#   log_post_tool_use.sh records what FIRED (rules whose globs hit).
#
# Across many sessions, the ratio of "available but never fired"
# rules tells us which path-scoped rules are dead weight. The
# /audit-bundle quarterly review surfaces these for pruning.
#
# All output goes to .claude/observability/rule_usage.jsonl as
# append-only JSONL. Same gitignored directory as session start logs.
#
# This hook fires on every tool invocation, so it must be cheap — the
# whole script runs in well under 50ms (one Python invocation, no
# subprocess fan-out). The matching logic is conservative: if we
# can't determine the file_path for any reason, we silently no-op
# rather than logging a partial record.

set -e

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"
OBS_DIR="${PROJECT_DIR}/.claude/observability"
LOG_FILE="${OBS_DIR}/rule_usage.jsonl"

mkdir -p "${OBS_DIR}"

# Capture the JSON payload from stdin to a temp file so the heredoc-quoted
# Python script below can reliably read it (the bash heredoc syntax otherwise
# competes with the python3 - stdin-read for the same FD).
PAYLOAD_FILE=$(mktemp)
trap "rm -f '${PAYLOAD_FILE}'" EXIT
cat > "${PAYLOAD_FILE}"

python3 - "${PROJECT_DIR}" "${LOG_FILE}" "${PAYLOAD_FILE}" <<'PY'
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

project_dir = Path(sys.argv[1])
log_file = Path(sys.argv[2])
payload_file = Path(sys.argv[3])

try:
    payload = json.loads(payload_file.read_text())
except (json.JSONDecodeError, ValueError, OSError):
    sys.exit(0)

tool_name = payload.get("tool_name", "")
tool_input = payload.get("tool_input", {})

# Only Edit/Write/MultiEdit/Read affect rule loading. Other tools
# (Bash, Grep, Glob) don't currently trigger path-scoped rule loads.
if tool_name not in {"Edit", "Write", "MultiEdit", "Read"}:
    sys.exit(0)

file_path = tool_input.get("file_path", "")
if not file_path:
    sys.exit(0)

# Make file_path relative to project root for glob matching.
try:
    rel_path = str(Path(file_path).resolve().relative_to(project_dir.resolve()))
except (ValueError, OSError):
    rel_path = file_path

# Read each rule's frontmatter and extract its `paths:` glob.
rules_dir = project_dir / ".claude" / "rules"
matched_rules = []
if rules_dir.is_dir():
    for rule_file in sorted(rules_dir.glob("*.md")):
        if rule_file.name == "README.md":
            continue
        text = rule_file.read_text(errors="replace")
        m = re.match(r"^---\n(.*?)\n---", text, re.DOTALL)
        if not m:
            continue
        frontmatter = m.group(1)
        # Extract paths: list (YAML-ish; we don't need full YAML parsing)
        paths_section = re.search(r"^paths:\n((?:\s+-\s+.*\n?)+)", frontmatter, re.MULTILINE)
        if not paths_section:
            continue
        # Each line is `  - "pattern"` or `  - pattern`
        for line in paths_section.group(1).splitlines():
            line = line.strip()
            if not line.startswith("-"):
                continue
            pattern = line[1:].strip().strip('"').strip("'")
            # Convert glob to regex. ** matches anything (incl /), * matches one segment.
            regex_str = re.escape(pattern)
            regex_str = regex_str.replace(r"\*\*", ".*").replace(r"\*", "[^/]*")
            if re.match(regex_str + r"$", rel_path) or re.match(regex_str, rel_path):
                matched_rules.append(rule_file.name)
                break

# Skip writing if nothing matched — keeps the log focused on
# actual rule fires.
if not matched_rules:
    sys.exit(0)

session_id = os.environ.get("CLAUDE_SESSION_ID", "unknown")
record = {
    "event": "rule_match",
    "session_id": session_id,
    "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "tool_name": tool_name,
    "file_path": rel_path,
    "matched_rules": matched_rules,
}

with log_file.open("a") as f:
    f.write(json.dumps(record) + "\n")
PY

exit 0
