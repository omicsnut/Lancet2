#!/usr/bin/env bash
# Lancet2 SessionStart observability hook: log_session_start
#
# Silent file logging at session start. Captures three pieces of state
# for the quarterly /audit-bundle review:
#
#   1. Which path-scoped rules under .claude/rules/ are AVAILABLE
#      (their descriptions and globs). The actual loading happens
#      lazily on glob match; this hook just records what *could*
#      load this session. Pairs with log_post_tool_use.sh which
#      records what *did* load.
#
#   2. Whether AGENTS.md @import resolution succeeded — i.e., whether
#      the @AGENTS.md import in CLAUDE.md actually pulled the canonical
#      content. A failed import would silently drop the project memory,
#      which would be hard to detect without this log.
#
#   3. Approximate initial context size — the line count of AGENTS.md
#      plus the descriptions of all available rules and skills. This
#      is the "always-loaded" baseline; pairing this number against
#      observed session-end token usage helps calibrate cost-model.md's
#      claims.
#
# All output goes to .claude/observability/sessions.jsonl as
# append-only JSONL. The directory is gitignored so logs don't pollute
# git status. The /audit-bundle command parses this file during
# quarterly review.
#
# Why a separate hook from session_start.sh: that script does
# user-facing reminders (chr1/chr4 test fixtures, layer direction).
# This hook is silent observability with no console output. Mixing
# the two would either spam the user with debug info or hide
# observability data behind user-output gating.

set -e

# ── Resolve project root and observability directory ─────────────────────
PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"
OBS_DIR="${PROJECT_DIR}/.claude/observability"
LOG_FILE="${OBS_DIR}/sessions.jsonl"

# Create the dir on first session. Idempotent.
mkdir -p "${OBS_DIR}"

# ── Gather session metadata ──────────────────────────────────────────────
SESSION_ID="${CLAUDE_SESSION_ID:-unknown}"
TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

# AGENTS.md presence + line count (canonical content size).
if [ -f "${PROJECT_DIR}/AGENTS.md" ]; then
    AGENTS_LINES=$(wc -l < "${PROJECT_DIR}/AGENTS.md" | tr -d ' ')
else
    AGENTS_LINES=0
fi

# CLAUDE.md presence + import-target check. Verify the wrapper
# imports AGENTS.md (the @AGENTS.md line); if not, the session is
# running without the canonical project memory.
CLAUDE_IMPORT_OK=false
if [ -f "${PROJECT_DIR}/CLAUDE.md" ]; then
    if grep -q "^@AGENTS\.md" "${PROJECT_DIR}/CLAUDE.md"; then
        CLAUDE_IMPORT_OK=true
    fi
fi

# Available rules: list each .claude/rules/<layer>.md (excluding README).
AVAILABLE_RULES=()
if [ -d "${PROJECT_DIR}/.claude/rules" ]; then
    for f in "${PROJECT_DIR}"/.claude/rules/*.md; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        [ "$base" = "README.md" ] && continue
        AVAILABLE_RULES+=("$base")
    done
fi

# Available skills count.
SKILL_COUNT=0
if [ -d "${PROJECT_DIR}/.claude/skills" ]; then
    SKILL_COUNT=$(find "${PROJECT_DIR}/.claude/skills" -name SKILL.md -type f | wc -l | tr -d ' ')
fi

# ── Emit a JSONL record ──────────────────────────────────────────────────
# Use Python to safely encode the record. jq is not guaranteed.
python3 - "$SESSION_ID" "$TIMESTAMP" "$AGENTS_LINES" "$CLAUDE_IMPORT_OK" \
        "$SKILL_COUNT" "${AVAILABLE_RULES[@]}" <<'PY' >> "${LOG_FILE}"
import json, sys
session_id, ts, agents_lines, import_ok, skill_count, *rules = sys.argv[1:]
record = {
    "event": "session_start",
    "session_id": session_id,
    "timestamp": ts,
    "agents_md_lines": int(agents_lines),
    "claude_md_import_ok": import_ok == "true",
    "available_rules": rules,
    "available_skills_count": int(skill_count),
}
print(json.dumps(record))
PY

# Always exit 0 — observability never blocks.
exit 0
