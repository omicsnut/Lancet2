---
description: Run the full Lancet2 validation suite on demand (build-debug + test + iwyu-fix + lint-check). Same as the Stop hook, but invokable mid-session.
allowed-tools: Bash
---

# /check — manual full validation

Run the full validation suite without waiting for the Stop hook to fire. This is the same sequence the `stop_validate.sh` hook runs: `pixi run build-debug`, then `pixi run test`, then `pixi run iwyu-fix`, then `pixi run lint-check`. It takes five to fifteen minutes and catches everything CI would catch.

Use this command in three situations. First, after completing a meaningful chunk of work mid-session, when you want validation feedback before continuing rather than at the end. Second, when you have set `CLAUDE_LANCET_SKIP_STOP_VALIDATE=1` to bypass the Stop hook for a fast-iteration session and want to re-enable validation on demand. Third, before invoking `fresh-reviewer` or pushing to a remote, as a sanity check that the change is in a green state before review.

The output format is staged: each of the four steps reports success or failure individually, with a tail of the relevant log on failure. A failed build halts the sequence (no point in running tests, IWYU, or lint against a broken build). A failed test, IWYU, or lint step reports its own output but does not halt; all three signals are useful, and you want them all in one /check run.

The third step, `iwyu-fix`, mutates files when it finds violations — running IWYU's `--fix` and then re-applying clang-format. If IWYU rewrites any include block, those changes need to be staged and committed alongside the rest of the work. The trailing lint-check sees the post-IWYU state, so its clang-tidy diagnostics reflect what would actually go to CI. The pre-commit hook `iwyu_check_on_commit.sh` runs `iwyu-check` (read-only) at commit time as a final gate; if /check passes locally, the pre-commit hook will too.

## What this command does

```bash
echo "─── /check validation ─────────────────────────────────────"

# Tempfiles for each phase. After all four steps complete, clean up
# the tempfiles. Each `rm` will be validated by block_dangerous_bash
# (paths are under /tmp/, so they pass) and Claude Code will surface
# a one-click approval prompt before the actual deletion runs. The
# extra prompt is acceptable: /check is a deliberate user-initiated
# command and the cleanup keeps /tmp tidy across many runs.

build_log=$(mktemp)
echo "1/4 build-debug ..."
pixi run --quiet build-debug > "$build_log" 2>&1
build_rc=$?
if [ $build_rc -ne 0 ]; then
  echo "❌ build-debug failed:"
  tail -40 "$build_log"
  exit 1
fi
echo "✓ build-debug ok"

test_log=$(mktemp)
echo "2/4 test ..."
pixi run --quiet test > "$test_log" 2>&1
test_rc=$?
if [ $test_rc -ne 0 ]; then
  echo "❌ tests failed:"
  tail -40 "$test_log"
else
  passed=$(grep -oE '[0-9]+ assertion' "$test_log" | head -1)
  echo "✓ tests ok ($passed)"
fi

iwyu_log=$(mktemp)
echo "3/4 iwyu-fix ..."
pixi run --quiet iwyu-fix > "$iwyu_log" 2>&1
iwyu_rc=$?
if [ $iwyu_rc -ne 0 ]; then
  echo "❌ iwyu-fix failed:"
  tail -40 "$iwyu_log"
elif [ -n "$(git status --porcelain -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx')" ]; then
  echo "⚠ iwyu-fix rewrote includes — review and stage:"
  git status --short -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx'
else
  echo "✓ iwyu-fix ok (no changes)"
fi

lint_log=$(mktemp)
echo "4/4 lint-check ..."
pixi run --quiet lint-check > "$lint_log" 2>&1
lint_rc=$?
if [ $lint_rc -ne 0 ]; then
  echo "❌ clang-tidy violations:"
  tail -30 "$lint_log"
else
  echo "✓ lint-check ok"
fi

if [ -n "$(git status --porcelain)" ]; then
  echo ""
  echo "⚠ uncommitted changes:"
  git status --short
fi

# Clean up the tempfiles. Each rm targets a single absolute path under
# /tmp/, so block_dangerous_bash will validate and Claude Code will
# prompt for approval. -f silences "no such file" if a phase aborted
# before its log was created.
rm -f "$build_log"
rm -f "$test_log"
rm -f "$iwyu_log"
rm -f "$lint_log"

echo "────────────────────────────────────────────────────────────"
```

## When you want narrower test runs

`/check` runs all four CI gates against the full test set. It's the right tool for "is this change merge-ready?" but it's overkill for "did my one change break the test I just wrote?". For narrower invocations, drive the test binary directly:

```bash
# Just the test you're iterating on
./cmake-build-debug/tests/TestLancet2 "test name"

# Just one layer
./cmake-build-debug/tests/TestLancet2 "[caller]"

# Reproduce a random-order failure (seed printed in the failing /check output)
./cmake-build-debug/tests/TestLancet2 --rng-seed <reported-seed>
```

The `add-cpp-test` skill has the full Catch2 invocation reference (tag filtering, sections, generators, reporters, sharding, seed reproduction). Use `/check` for the validation gauntlet; use the binary directly for everything else.
