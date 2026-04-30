#!/bin/bash
# Lancet2 Stop hook: stop_validate
# Runs the full validation suite when Claude declares a session done:
#   1. pixi run build-debug    (compiles in cmake-build-debug)
#   2. pixi run test           (runs Catch2 unit tests)
#   3. pixi run iwyu-fix       (IWYU --fix, then clang-format --fix)
#   4. pixi run lint-check     (clang-tidy, fails on warnings)
# This takes 5-15 minutes on a real codebase. The user accepted that
# trade-off in exchange for catching everything CI would catch.
#
# To skip validation for a fast iteration session, set:
#   export CLAUDE_LANCET_SKIP_STOP_VALIDATE=1
# in your shell before invoking `claude`. The /check slash command runs
# the same validation on demand without needing this flag.

set +e

if [ -n "$CLAUDE_LANCET_SKIP_STOP_VALIDATE" ]; then
  echo "ℹ Stop validation skipped (CLAUDE_LANCET_SKIP_STOP_VALIDATE is set)."
  if [ -n "$(git status --porcelain)" ]; then
    echo "⚠ uncommitted changes:"
    git status --short
  fi
  exit 0
fi

echo "─── Stop validation ────────────────────────────────────────"
echo "Running build-debug + test + iwyu-fix + lint-check via pixi."
echo "This may take 5-15 minutes. To skip in future sessions:"
echo "  export CLAUDE_LANCET_SKIP_STOP_VALIDATE=1"
echo ""

# Build first; test depends on it. The test target depends on build-debug
# anyway via pixi task dependencies, but running it explicitly gives clearer
# output if it fails.
echo "1/4 build-debug ..."
build_log=$(mktemp)
pixi run --quiet build-debug > "$build_log" 2>&1
build_rc=$?
if [ $build_rc -ne 0 ]; then
  echo "❌ build-debug failed:"
  tail -40 "$build_log"
  rm -f "$build_log"
  exit 0   # exit 0 because we are reporting, not blocking
fi
echo "✓ build-debug ok"
rm -f "$build_log"

echo "2/4 test ..."
test_log=$(mktemp)
pixi run --quiet test > "$test_log" 2>&1
test_rc=$?
if [ $test_rc -ne 0 ]; then
  echo "❌ tests failed:"
  tail -40 "$test_log"
else
  passed=$(grep -oE '[0-9]+ assertion' "$test_log" | head -1)
  echo "✓ tests ok ($passed)"
fi
rm -f "$test_log"

# iwyu-fix mutates files when it finds violations. The downstream
# lint-check then sees the post-IWYU state, which is what CI would.
# A clean iwyu-fix run leaves the tree untouched.
echo "3/4 iwyu-fix ..."
iwyu_log=$(mktemp)
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
rm -f "$iwyu_log"

echo "4/4 lint-check ..."
lint_log=$(mktemp)
pixi run --quiet lint-check > "$lint_log" 2>&1
lint_rc=$?
if [ $lint_rc -ne 0 ]; then
  echo "❌ clang-tidy violations:"
  tail -30 "$lint_log"
else
  echo "✓ lint-check ok"
fi
rm -f "$lint_log"

if [ -n "$(git status --porcelain)" ]; then
  echo ""
  echo "⚠ uncommitted changes:"
  git status --short
fi

echo "────────────────────────────────────────────────────────────"
exit 0
