#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: validate_commit_message

Intercepts `git commit -m "..."` (and `-am`, `--message`, `--message=`) bash
commands and validates the message against the rules in .claude/commit-style.json.

The rules are grounded in the project's actual .chglog/config.yml. The chglog
header pattern is `^(\\w*)\\:\\s(.*)$` with pattern_maps [Type, Subject] — type
and subject only, no scope group, exactly one space after the colon. The filter
list is exactly [feat, fix, perf, chore]; any other type silently fails to
appear in the generated CHANGELOG, so this hook rejects it pre-commit.

Hard blocks (exit 2):
  - subject does not match the chglog header pattern (no colon, wrong spacing,
    or scope syntax which chglog does not support)
  - type is not in allowed_types (would not appear in CHANGELOG)
  - type is not lowercase (chglog filter is case-sensitive)
  - subject too long, too short, ends with period, or starts with uppercase
  - type requires a body (feat/fix/perf) but none was given

Soft warns (printed to stderr, exit 0):
  - subject's first word suggests non-imperative mood
  - body line exceeds the configured wrap limit
  - fixup!/squash! markers (should be rebased before push)

Skips entirely (exit 0 silently):
  - any non-commit bash command
  - `git commit` without -m (opens $EDITOR, out of scope; use a real
    prepare-commit-msg git hook for that surface)
  - `git commit --amend` reusing the prior message
  - merge and revert commits matching configured prefixes
"""
import json
import os
import re
import shlex
import sys
from pathlib import Path

# Heuristic: words at the start of a subject that suggest non-imperative mood.
# Intentionally narrow to keep false-positive rate low.
NON_IMPERATIVE_HINTS = {
    "added", "adds", "adding",
    "fixed", "fixes", "fixing",
    "updated", "updates", "updating",
    "changed", "changes", "changing",
    "removed", "removes", "removing",
    "refactored", "refactoring",
    "improved", "improves", "improving",
    "implemented", "implementing",
    "created", "creates", "creating",
    "deleted", "deletes", "deleting",
    "renamed", "renaming",
    "moved", "moving",
}


def load_config(project_dir: str) -> dict | None:
    """Load commit-style config from project root, then user home."""
    candidates = [
        Path(project_dir) / ".claude" / "commit-style.json",
        Path.home() / ".claude" / "commit-style.json",
    ]
    for path in candidates:
        if path.is_file():
            try:
                with path.open() as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError):
                # Don't block on a broken config; let the user see the error elsewhere.
                return None
    return None


def extract_commit_messages(command: str) -> list[str] | None:
    """
    Parse a bash command string and return the list of -m / --message values
    if it is a `git commit` invocation. Returns None if not a git commit at all,
    or [] if it is `git commit` without -m (opens $EDITOR — out of scope).

    Multiple -m flags concatenate as separate paragraphs (git's documented
    behavior). The first -m is the subject; subsequent -m flags become the body.
    Compound shell commands like `cd /repo && git commit -m "..."` are handled
    by scanning for `git ... commit` and stopping at shell separators.
    """
    try:
        tokens = shlex.split(command, posix=True)
    except ValueError:
        # Unclosed quote or similar — let the actual git invocation surface it.
        return None
    if not tokens:
        return None

    messages: list[str] = []
    found_commit = False
    SHELL_SEP = ("&&", "||", ";", "|", "&")

    i = 0
    while i < len(tokens):
        if tokens[i] == "git":
            j = i + 1
            # Walk forward looking for `commit`, allowing `git -C dir commit` etc.
            while j < len(tokens) and tokens[j] not in SHELL_SEP:
                if tokens[j] == "commit":
                    found_commit = True
                    k = j + 1
                    while k < len(tokens) and tokens[k] not in SHELL_SEP:
                        ct = tokens[k]
                        if ct in ("-m", "--message"):
                            if k + 1 < len(tokens):
                                messages.append(tokens[k + 1])
                                k += 2
                                continue
                        elif ct.startswith("--message="):
                            messages.append(ct[len("--message="):])
                        elif ct in ("-am", "-ma"):
                            # Combined -a + -m
                            if k + 1 < len(tokens):
                                messages.append(tokens[k + 1])
                                k += 2
                                continue
                        k += 1
                    break
                j += 1
        i += 1

    if not found_commit:
        return None
    return messages


def is_auto_generated(subject: str, allowed_prefixes: list[str]) -> bool:
    return any(subject.startswith(p) for p in allowed_prefixes)


def is_fixup_or_squash(subject: str) -> bool:
    return subject.startswith("fixup!") or subject.startswith("squash!")


def validate_message(messages: list[str], cfg: dict, command: str = "") -> tuple[list[str], list[str]]:
    """
    Validate the parsed commit messages against the config.
    Returns (hard_errors, soft_warnings). Hard errors block; soft warnings print only.
    """
    hard: list[str] = []
    soft: list[str] = []

    if not messages:
        # No -m means git will open $EDITOR; out of scope for this hook.
        return hard, soft

    subject = messages[0].strip()
    body_paragraphs = [m.strip() for m in messages[1:] if m.strip()]

    # Skip auto-generated merge/revert commits entirely.
    auto_prefixes = cfg.get("auto_generated_subject_prefixes_allowed", [])
    if is_auto_generated(subject, auto_prefixes):
        return hard, soft

    # Soft-warn on fixup!/squash! commits.
    if is_fixup_or_squash(subject):
        if cfg.get("fixup_squash_handling", "soft_warn") == "soft_warn":
            soft.append(
                "Subject is a fixup!/squash! marker. These should be rebased away "
                "(`git rebase -i --autosquash`) before push."
            )
        return hard, soft

    # ── Header parsing — match chglog's pattern exactly ──────────────────
    pattern_str = cfg.get("header_pattern", r"^(?P<type>\w+): (?P<subject>.+)$")
    try:
        header_re = re.compile(pattern_str)
    except re.error as e:
        # Misconfigured pattern; fail open with a stderr note.
        print(f"⚠ commit-style: invalid header_pattern in config: {e}", file=sys.stderr)
        return hard, soft

    m = header_re.match(subject)
    if not m:
        hard.append(
            f"Subject does not match the project header pattern from .chglog/config.yml.\n"
            f"  Pattern: {pattern_str}\n"
            f"  Got:     {subject!r}\n"
            f"  Required form: '<type>: <subject>' with exactly one space after the colon.\n"
            f"  Note: scopes like 'fix(caller): ...' are NOT supported by .chglog/config.yml.\n"
            f"  Examples of correct subjects:\n"
            + "\n".join(f"    {ex}" for ex in cfg.get("examples_good", [])[:3])
        )
        return hard, soft

    # Pull captures by group name when available, falling back to position.
    # The default config uses named groups; a user-overridden pattern from
    # .chglog/config.yml ported verbatim might use numbered groups instead.
    try:
        ctype = m.group("type")
        subj_text = m.group("subject")
    except IndexError:
        groups = m.groups()
        ctype = groups[0] if len(groups) >= 1 else ""
        subj_text = groups[1] if len(groups) >= 2 else ""

    # ── Type validation ──────────────────────────────────────────────────
    allowed_types = cfg.get("allowed_types", [])
    if ctype != ctype.lower():
        hard.append(
            f"Type must be lowercase. Got '{ctype}'. The chglog filter is "
            f"case-sensitive, so '{ctype}' would silently fail to be grouped "
            f"into the changelog."
        )
    elif allowed_types and ctype not in allowed_types:
        # Build a helpful message that explains why this matters.
        sections = cfg.get("type_to_changelog_section", {})
        section_hints = ", ".join(
            f"{t}={sections.get(t, '?')!r}" for t in allowed_types
        )
        hard.append(
            f"Type '{ctype}' is not in .chglog/config.yml's filter list.\n"
            f"  Allowed types and their changelog sections: {section_hints}.\n"
            f"  A '{ctype}:' commit would parse but never appear in the generated "
            f"CHANGELOG.md. If this is a refactor, use 'chore:'. If it is a doc "
            f"or test or build change that genuinely should not appear in the "
            f"changelog, use 'chore:' or extend allowed_types in commit-style.json."
        )

    # ── Subject text validation ──────────────────────────────────────────
    max_len = cfg.get("max_subject_length", 72)
    if len(subject) > max_len:
        hard.append(
            f"Subject line is {len(subject)} characters; limit is {max_len}. "
            f"Shorten the subject and move detail to the body."
        )

    min_len = cfg.get("min_subject_length", 10)
    if len(subj_text) < min_len:
        hard.append(
            f"Subject text is {len(subj_text)} characters after the type; "
            f"minimum is {min_len}. Be specific about what changed."
        )

    if cfg.get("subject_must_not_end_with_period", True) and subj_text.endswith("."):
        hard.append("Subject must not end with a period.")

    if cfg.get("subject_must_be_lowercase_first_letter", False) and subj_text:
        first = subj_text[0]
        if first.isupper():
            hard.append(
                f"Subject text must start with a lowercase letter. Got '{first}' "
                f"({subj_text[:30]!r}...)."
            )

    # Imperative-mood heuristic.
    imperative_setting = cfg.get("subject_imperative_mood", "soft_warn")
    if imperative_setting in ("soft_warn", "hard_block") and subj_text:
        words = subj_text.split()
        first_word = words[0].lower() if words else ""
        if first_word in NON_IMPERATIVE_HINTS:
            msg = (
                f"Subject's first word '{first_word}' suggests non-imperative mood. "
                f"Conventional commits use imperative ('add' not 'added', "
                f"'fix' not 'fixes')."
            )
            if imperative_setting == "hard_block":
                hard.append(msg)
            else:
                soft.append(msg)

    # ── Body validation ───────────────────────────────────────────────────
    # Per docs_dev/style/cpp_style.md § Git commit messages, trivial commits
    # (typos, exec-bit, formatting, dep bumps) ship with no body; substantial
    # commits use a two-section body.
    #
    # We classify "substantial" using two signals, both regardless of type:
    #   1. Subject pattern — matches one of trivial_commit_subject_patterns.
    #   2. Diff size — insertions + deletions in the staged diff.
    #
    # Combinations:
    #   trivial + small diff   + no body → silent pass
    #   trivial + large diff   + no body → soft-warn (subject claims trivial,
    #                                       but diff is large; double-check)
    #   non-trivial + small    + no body → soft-warn (consider a body)
    #   non-trivial + large    + no body → HARD BLOCK
    #   any        + body      + any     → run body-shape checks
    #
    # If the diff size is unknowable (no git, not in a repo, etc.), we fail
    # open: treat the commit as if it were small. Better to occasionally let
    # a too-big-no-body slip through than to block legitimate commits in
    # weird repo states.

    is_trivial = _matches_trivial_pattern(subject, cfg)
    threshold = cfg.get("small_diff_line_threshold", 30)

    if not body_paragraphs:
        diff_size = _get_diff_size(command)
        is_small = diff_size is None or diff_size <= threshold

        if is_trivial:
            if not is_small:
                soft.append(
                    f"Subject matches a trivial pattern ('{subject}') but the staged "
                    f"diff is {diff_size} lines (threshold {threshold}). If this is "
                    f"actually a substantive change, add a body explaining what and why."
                )
            # else: silent pass — typical trivial commit
        else:
            # Subject is not trivial; body is expected unless diff is small.
            if is_small:
                size_phrase = (
                    "diff size unavailable" if diff_size is None
                    else f"only {diff_size} lines changed"
                )
                soft.append(
                    f"Subject does not match a trivial pattern and {size_phrase}. "
                    f"Consider adding a body explaining what changed and why; even "
                    f"a one-paragraph context note helps reviewers and future-you."
                )
            else:
                hard.append(
                    f"Substantive change ({diff_size} lines) without a body.\n"
                    f"  Subject does not match a trivial pattern (typos, exec-bit, format,\n"
                    f"  dep bumps), and the staged diff exceeds the {threshold}-line threshold\n"
                    f"  beyond which the subject alone is unlikely to convey what changed.\n"
                    f"  Add a body with a context paragraph and file-list bullets:\n"
                    f"    git commit -m \"{subject}\" \\\n"
                    f"      -m \"Context paragraph explaining why...\" \\\n"
                    f"      -m \"- src/lancet/<layer>/file.cpp: what changed there\"\n"
                    f"  If the change really is trivial, use a subject pattern that\n"
                    f"  matches (e.g. 'chore: format ...', 'chore: bump ...', 'fix: typo ...')."
                )

    # Body line-length check applies regardless of shape.
    max_body_line = cfg.get("max_body_line_length", 100)
    for paragraph in body_paragraphs:
        for line in paragraph.split("\n"):
            if len(line) > max_body_line:
                soft.append(
                    f"Body line exceeds {max_body_line} chars: "
                    f"{line[:60]}... ({len(line)} chars)"
                )
                break  # one warning per paragraph is enough

    # Body-shape check fires only when a body is present and is non-trivial.
    if body_paragraphs and not is_trivial:
        shape_warnings = _check_body_shape(body_paragraphs, cfg)
        if cfg.get("body_shape_handling", "soft_warn") == "hard_block":
            hard.extend(shape_warnings)
        else:
            soft.extend(shape_warnings)

    return hard, soft


def _matches_trivial_pattern(subject: str, cfg: dict) -> bool:
    """True if subject matches any configured trivial-commit pattern."""
    patterns = cfg.get("trivial_commit_subject_patterns", [])
    for pat in patterns:
        try:
            if re.match(pat, subject, re.IGNORECASE):
                return True
        except re.error:
            continue  # skip malformed pattern silently
    return False


def _get_diff_size(command: str) -> int | None:
    """
    Return the staged diff size as insertions + deletions, suitable for the
    body-required threshold check.

    For `git commit -m`, this is `git diff --cached --shortstat` — the changes
    already in the index. For `git commit -am` (or -ma, --all), the `-a` flag
    has not yet executed at hook time; we add `git diff --shortstat` (working
    tree vs index) so the count reflects what `-a` would stage. Untracked
    files are not included in either case (they are not picked up by `-a`).

    Returns None if the diff size cannot be determined (no git binary, not
    in a repo, command timed out). Callers should treat None as "unknown,
    fail open" — a missing measurement should not block a commit.
    """
    import subprocess

    is_all_flag = bool(re.search(r"(?:^|\s)(?:-a|-am|-ma|--all)(?:\s|$)", command))

    def _shortstat_to_lines(stdout: str) -> int:
        """Parse output like ' 3 files changed, 42 insertions(+), 5 deletions(-)'."""
        if not stdout.strip():
            return 0
        ins = 0
        dels = 0
        m = re.search(r"(\d+)\s+insertions?\(\+\)", stdout)
        if m:
            ins = int(m.group(1))
        m = re.search(r"(\d+)\s+deletions?\(-\)", stdout)
        if m:
            dels = int(m.group(1))
        return ins + dels

    try:
        cached = subprocess.run(
            ["git", "diff", "--cached", "--shortstat"],
            capture_output=True, text=True, timeout=2, check=False,
        )
        if cached.returncode != 0:
            return None
        total = _shortstat_to_lines(cached.stdout)

        if is_all_flag:
            unstaged = subprocess.run(
                ["git", "diff", "--shortstat"],
                capture_output=True, text=True, timeout=2, check=False,
            )
            if unstaged.returncode == 0:
                total += _shortstat_to_lines(unstaged.stdout)
            # If unstaged probe fails, return what we have (cached only) rather than None.
        return total
    except (subprocess.SubprocessError, OSError, FileNotFoundError):
        return None


def _check_body_shape(body_paragraphs: list[str], cfg: dict) -> list[str]:
    """
    Enforce the two-section body shape from docs_dev/style/cpp_style.md
    § Git commit messages.

    The body is interpreted as:
      [context paragraph(s)] [optional language sub-header + file-list bullets]+

    What we check:
      - If file-list bullets are present at all, each must match the bullet pattern.
      - If a language sub-header line is present, the lines after it (until the
        next sub-header or end) should be file-list bullets.
      - If the body has 3+ paragraphs and zero file-list bullets, recommend
        adding the file list — at that volume the bullets aid review.

    What we do NOT check:
      - Whether context paragraphs use backticks (style preference, not
        structural).
      - Whether every changed file is mentioned in the file list (we don't
        have the diff).
      - Whether the brace shorthand `name.{h,cpp}` is used (the bullet
        pattern allows it but doesn't require it).
    """
    warnings: list[str] = []

    bullet_pat_str = cfg.get(
        "file_list_bullet_pattern",
        r"^- (?P<path>[\w./{}, *-]+):\s+(?P<change>.+)$",
    )
    subheader_pat_str = cfg.get(
        "language_subheader_pattern",
        r"^(?:C\+\+|Python|Shell|Bash|CMake|YAML|JSON|Markdown|Rust|Go) changes:$",
    )
    try:
        bullet_re = re.compile(bullet_pat_str)
        subheader_re = re.compile(subheader_pat_str)
    except re.error:
        return warnings  # fail open on misconfigured patterns

    # Flatten paragraphs into lines for sequential analysis.
    all_lines: list[str] = []
    for p in body_paragraphs:
        for line in p.split("\n"):
            stripped = line.rstrip()
            if stripped:
                all_lines.append(stripped)

    bullet_lines: list[str] = []
    malformed_bullets: list[str] = []
    in_file_list_block = False  # True after we see a sub-header or first bullet

    for line in all_lines:
        if subheader_re.match(line):
            in_file_list_block = True
            continue

        # A line starting with "- " in the file-list region: must be a valid bullet.
        if line.startswith("- "):
            in_file_list_block = True
            if bullet_re.match(line):
                bullet_lines.append(line)
            else:
                malformed_bullets.append(line)
        elif in_file_list_block:
            # In a file-list region but not a bullet — usually a continuation
            # line or a back-to-prose. Don't flag; the writer may have added
            # a closing remark. Reset the region tracker so a subsequent
            # rogue "- ..." line is still caught above.
            pass

    if malformed_bullets:
        warnings.append(
            "File-list bullets do not match the expected `- path: change` format:\n"
            + "\n".join(f"    {b[:80]}" for b in malformed_bullets[:3])
            + "\n  Expected form: `- src/lancet/<layer>/<file>: <what changed>` "
            "or `- name.{h,cpp}: ...` for header/source pairs."
        )

    if len(body_paragraphs) >= 3 and not bullet_lines:
        warnings.append(
            "Body has 3+ paragraphs but no file-list bullets. "
            "docs_dev/style/cpp_style.md § Git commit messages recommends "
            "a `- path: change` section for substantial commits to aid "
            "review and bisection."
        )

    return warnings


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        # Hook protocol changed or input malformed; fail open.
        return 0

    command = payload.get("tool_input", {}).get("command", "") or ""
    if not command:
        return 0

    messages = extract_commit_messages(command)
    if messages is None:
        # Not a `git commit` command at all.
        return 0
    if not messages:
        # `git commit` with no -m — opens $EDITOR. Out of scope.
        return 0

    project_dir = os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())
    cfg = load_config(project_dir)
    if cfg is None:
        # No config; nothing to enforce. Fail open.
        return 0

    hard, soft = validate_message(messages, cfg, command=command)

    for w in soft:
        print(f"⚠ commit-style: {w}", file=sys.stderr)

    if hard:
        print(
            "BLOCKED: commit message does not satisfy Lancet2 commit-style rules.",
            file=sys.stderr,
        )
        for e in hard:
            print(f"  - {e}", file=sys.stderr)
        print("", file=sys.stderr)
        print(
            "Rules live in .claude/commit-style.json (grounded in .chglog/config.yml).",
            file=sys.stderr,
        )
        print(
            "Run /commit to have Claude compose a conformant message from the staged diff.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
