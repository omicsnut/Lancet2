#!/bin/bash
# Lancet2 PostToolUse hook: post_edit_format
# Runs clang-format on changed C/C++ files via the project's pixi task.
# Delegating to pixi guarantees the tooling version matches CI exactly.
#
# CLAUDE_FILE_PATHS is space-separated; it contains every file Claude just
# edited. We filter to C/C++ extensions and pass to fmt-fix.

set +e

formatted=()
for f in $CLAUDE_FILE_PATHS; do
  case "$f" in
    *.cpp|*.cc|*.cxx|*.h|*.hpp|*.hxx)
      formatted+=("$f")
      ;;
  esac
done

if [ ${#formatted[@]} -eq 0 ]; then
  exit 0
fi

# pixi run fmt-fix runs scripts/run_clang_format.py --fix on the whole tree.
# Running the project task is more reliable than calling clang-format directly
# because it picks up the .clang-tidy-required clang-tools 22.1.4 from the
# pixi env. Output is suppressed unless something interesting happens.
pixi run --quiet fmt-fix 2>&1 | tail -10

# Always exit 0 — formatter failures should not halt Claude.
exit 0
