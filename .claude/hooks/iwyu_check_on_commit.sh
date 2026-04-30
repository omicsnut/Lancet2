#!/usr/bin/env bash
# Lancet2 PreToolUse hook: iwyu_check_on_commit
#
# Hard-blocks `git commit` invocations when IWYU finds include-hygiene
# violations, mirroring CI's iwyu-check gate but firing locally before
# the bad commit lands.
#
# Why a hook rather than a slash command: the gate must be enforced even
# when Claude has been told to ignore it. Hooks override persuasion;
# AGENTS.md cannot.
#
# Why PreToolUse on git commit specifically: IWYU is a whole-tree analysis
# that runs the clang frontend across compile_commands.json — minutes per
# invocation, not seconds. Firing per-edit would tax the feedback loop
# enough that the user (and Claude) would set CLAUDE_LANCET_SKIP_*
# overrides to bypass it. Firing only when Claude is about to commit
# means the cost lands at exactly the moment the user has clearly
# decided "this work is ready to gate," matching CI's lifecycle moment.
#
# Why this is separate from validate_commit_message and pre_commit_summary:
# they validate the message and summarize the diff, respectively. This
# hook validates source quality. Three distinct concerns, three hooks,
# each easier to maintain in isolation.
#
# Architecture: the command-parsing fast-path uses the same inline-python
# JSON-read trick as pre_commit_summary.sh so we don't pay a heavy
# parser cost on every Bash invocation. The actual iwyu-check runs only
# when we have confirmed this is a `git commit` worth gating.

set -e

# ── Read JSON payload from stdin and extract the command ──────────────────
payload=$(cat)
command=$(printf '%s' "$payload" | python3 -c '
import json, sys
try:
    data = json.load(sys.stdin)
    print(data.get("tool_input", {}).get("command", ""))
except Exception:
    pass
' 2>/dev/null)

# ── Fast path: only fire on actual git commit invocations ─────────────────
# Skip non-commit bash, plus `git commit --amend` reusing a prior message
# (no source change happens in that case; matches validate_commit_message.py's
# skip behavior). Compound commands like `cd /repo && git commit -m ...`
# match because the substring `git commit` appears in $command.
case "$command" in
    *"git commit"*) ;;
    *) exit 0 ;;
esac

# ── --amend without -m / --message reuses the prior subject; skip ─────────
# A re-amend that ALSO restages source changes is rare in solo workflows
# and is caught by the next deliberate commit. The conservative skip here
# matches validate_commit_message.py's pattern.
if printf '%s' "$command" | grep -qE -- '--amend' \
   && ! printf '%s' "$command" | grep -qE -- '(-m\b|-am\b|--message\b|--message=)'; then
    exit 0
fi

# ── Run IWYU and gate the commit on its exit code ─────────────────────────
echo "─── iwyu_check_on_commit ──────────────────────────────────"
echo "Running pixi run iwyu-check before commit ..."
echo ""

iwyu_log=$(mktemp)
pixi run --quiet iwyu-check > "$iwyu_log" 2>&1
iwyu_rc=$?

if [ $iwyu_rc -ne 0 ]; then
    cat "$iwyu_log" >&2
    echo "" >&2
    echo "────────────────────────────────────────────────────────────" >&2
    echo "IWYU found include-hygiene violations. Commit blocked." >&2
    echo "Apply fixes with:" >&2
    echo "    pixi run iwyu-fix" >&2
    echo "and re-stage the modified files before committing." >&2
    echo "────────────────────────────────────────────────────────────" >&2
    rm -f "$iwyu_log"
    exit 2
fi

echo "✓ iwyu-check ok"
echo "────────────────────────────────────────────────────────────"
rm -f "$iwyu_log"
exit 0
