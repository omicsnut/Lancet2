#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function ensure_deps() {
  if ! command -v git-chglog &> /dev/null; then
    echo "ERROR: git-chglog not found on PATH."
    echo "Run: pixi run ensure-git-chglog"
    exit 1
  fi
}

function main() {
  local PROJECT_DIR=$(git rev-parse --show-toplevel)

  echo "Updating changelog..."
  git-chglog --config "${PROJECT_DIR}/.chglog/config.yml" -o "${PROJECT_DIR}"/CHANGELOG.md
  
  git add "${PROJECT_DIR}"/CHANGELOG.md 
  git commit -vsm "chore: Update release notes in changelog"
  git push
  
  echo "Changelog updated and pushed successfully."
}

ensure_deps && main "$@"
