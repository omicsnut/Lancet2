#!/bin/bash
# Lancet2 SessionStart hook
# Prints orientation: branch, head commit, dirty status, build readiness.
# Must complete in under 1 second.

set +e

echo "─── Lancet2 session ────────────────────────────────────────"

branch=$(git branch --show-current 2>/dev/null)
head=$(git rev-parse --short HEAD 2>/dev/null)
subject=$(git log -1 --pretty=%s 2>/dev/null)

if [ -n "$branch" ]; then
  echo "branch:  $branch"
fi
if [ -n "$head" ]; then
  echo "head:    $head — $subject"
fi

if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
  echo "status:  ⚠ uncommitted changes (run: git status)"
fi

if [ ! -f cmake-build-debug/compile_commands.json ]; then
  echo "build:   ⚠ cmake-build-debug/compile_commands.json missing"
  echo "         run once: pixi run configure-debug"
else
  echo "build:   cmake-build-debug ready"
fi

# Check for project-significant recent changes that the agent should know about.
if [ -f CHANGELOG.md ]; then
  latest_version=$(grep -m1 -oE '\[v[0-9]+\.[0-9]+\.[0-9]+\]' CHANGELOG.md 2>/dev/null | head -1)
  if [ -n "$latest_version" ]; then
    echo "version: latest tag in CHANGELOG.md is $latest_version"
  fi
fi

echo "────────────────────────────────────────────────────────────"
exit 0
