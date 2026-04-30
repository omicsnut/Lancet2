# Hooks

Hooks are deterministic. They fire every time, can block actions before they execute (`PreToolUse` with exit code 2), and run outside Claude's context window so they cost zero tokens. Skills and `AGENTS.md` are probabilistic — Claude can ignore or forget them. Hooks cannot be ignored.

This is the entire reason hooks exist in this bundle: there are rules that must be enforced even when Claude has been told to ignore them, even when the rule is buried under context, even when Claude is focused on something else. The twelve hooks here are the rules that meet that bar.

## What's here

### Blocking hooks (PreToolUse, exit 2 on violation)

- **`block_protected_paths.py`** — refuses edits to `pixi.toml`, `cmake/`, `.github/`, `data/`, `pixi.lock`, `Dockerfile`, `cmake-build-*/`, `_deps/`, `.pixi/`, `CHANGELOG.md`. Substring match. The blocklist is conservative; `notes/scratch/` and `notes/<feature>/` are NOT blocked.
- **`block_dangerous_bash.py`** — validates bash commands before execution. Three layers: (1) **blanket-banned patterns** (`git push`, `git reset --hard origin/...`, `dd if=`, `mkfs.<fs>`, the fork-bomb literal, `chmod -R 777 /`, `chown -R`, `xargs rm`) — these have no safe variant; the maintainer handles them by hand outside Claude Code. `xargs rm` is in this list because its paths come from stdin and cannot be statically validated. (2) **path-validated deletions** (`rm`, `rmdir`, `shred`, `find … -delete`, `find … -exec rm`) — the hook parses out the path arguments and verifies each is absolute, contains no shell metacharacters that would re-expand at execution time (`$`, backtick, `*`, `?`, `[`, `]`, `{`, `}`), resolves cleanly (no `..` escape, no symlink to outside the allowlist), is not the bare allowlist root, and lands under `/tmp/`, `/var/tmp/`, `/scratch/`, OR inside `$CLAUDE_PROJECT_DIR` AND is gitignored (verified via `git check-ignore`). Any path that fails validation blocks the whole command. (3) **user-approval prompt** (downstream of the hook) — `settings.json` does NOT list `Bash(rm:*)` / `Bash(rmdir:*)` / `Bash(shred:*)` in either allow or deny, so Claude Code prompts the user per invocation. The hook prints a verdict line ("safe: under /tmp/", "safe: gitignored under project root") to stderr so the user-approval prompt surfaces the categorization the hook applied. Together this forms a two-key gate: the hook validates path safety AND the user approves each delete by hand. Word-boundary regex matching on the blanket-banned layer avoids false positives on words containing `rm` (firmware, storm, warm).
- **`validate_layer_direction.py`** — enforces the six-layer dependency rule from CMakeLists.txt: base → hts → cbdg → caller → core → cli. Blocks any include that goes against the chain.
- **`validate_naming.py`** — catches `mPascalCase` violations on members, `using namespace std` in headers, bare `assert()`, `std::format` and `std::print` (the project uses fmtlib via spdlog), and bare `// NOLINT` (any inline same-line form, with or without a check name; only scoped `NOLINTNEXTLINE` and `NOLINTBEGIN`/`NOLINTEND` are permitted).
- **`validate_commit_message.py`** — validates `git commit -m` messages against `commit-style.json` (which is grounded in `.chglog/config.yml` and `docs_dev/style/cpp_style.md` § Git commit messages). Hard-blocks invalid types, scopes, malformed subjects, and substantive changes (above the configured `small_diff_line_threshold`) without a body. Soft-warns on heuristic concerns.
- **`iwyu_check_on_commit.sh`** — gates `git commit` invocations on `pixi run iwyu-check`. Hard-blocks if IWYU finds include-hygiene violations (the same gate CI runs, just earlier). Skips `--amend` reusing a prior message. Prints full IWYU output plus a remediation hint pointing at `pixi run iwyu-fix`.

### Informational hooks (always exit 0; print to stderr)

- **`pre_commit_summary.sh`** — prints a one-paragraph summary of the staged change (file count, lines, layers touched, kinds of files) before a commit executes. Heuristic warnings for multi-layer-without-docs, source-without-tests, and VCF schema changes. Never blocks.
- **`session_start.sh`** — prints branch and status on session start. Useful orientation for solo-developer context-switching between branches.
- **`post_edit_format.sh`** — runs `pixi run fmt-fix` on edited files. Cheap, deterministic, removes a common source of CI friction.
- **`stop_validate.sh`** — runs the lint/test/build sweep on session stop. The 5-10 minute runtime is acceptable when the session was productive.

## Why hooks vs skills vs AGENTS.md

A rule is hookable when (a) it has zero false positives — every match is a real violation — and (b) you want it enforced regardless of what Claude is currently focused on.

- **Layer-direction violations:** hookable (the dependency chain is in CMakeLists, every upward include is wrong).
- **Naming violations** like bare `assert()` or `using namespace std` in headers: hookable (the rules are crisp, the source enforces them via clang-tidy).
- **Comment language quality:** not hookable (judgment calls about whether a word is filler in this specific context).
- **Choice of which algorithm to use:** not hookable (judgment about clarity vs. cleverness).

A rule that has any judgment call belongs in a skill or in AGENTS.md, not a hook. Hooks that try to encode judgment produce false positives that train Claude (and you) to ignore the hook output.

## Why use hooks at all when AGENTS.md could state the rule

AGENTS.md is read at session start (loaded by Claude Code through the CLAUDE.md wrapper) and at the start of each message, but it's not enforced — Claude can be persuaded out of any rule there by a strong enough context signal ("but in this case I should ignore the rule because…"). Hooks override the persuasion.

For the layer rule specifically: AGENTS.md explains the architecture so Claude understands the design; the hook prevents accidental violations even when Claude has correctly understood the design but slipped while writing.

## Why use hooks at all when CI catches violations

Feedback latency. CI catches the violation after you push; the hook catches it before the edit is even applied. The faster the feedback, the smaller the rework cost.

## Hook calibration philosophy

Every Lancet2 hook is calibrated to a near-zero-false-positive bar. The reason is empirical: if a hook produces false positives at any meaningful rate, both the human and Claude learn to ignore it. Once an enforcement mechanism is being routinely overridden, it has negative value — it adds noise to the session's stderr without preventing anything.

The bar is operationalized as three rules, applied to every hook in this directory:

1. **Hard blocks (exit 2) require near-zero false positives.** A pattern earns hard-block status only if (a) the violation is explicit and unambiguous in the source text, (b) reviewers have flagged it repeatedly in the past, and (c) the cost of an erroneous block is trivially recoverable (the user adds an inline comment override and re-runs).

2. **Soft warnings (exit 0 with stderr) are the right answer when the pattern is high-signal but ambiguous in some real cases.** Soft warnings MUST be acknowledged: when Claude sees a soft warning, the next assistant turn should briefly reference what was warned about (one phrase is enough — "noted, the m-prefix warning is on a struct field, intentional"). The acknowledgment keeps the warning's signal alive instead of letting it become wallpaper.

3. **Calibration is validated quarterly.** The `/audit-bundle` slash command reviews each hook's stderr output from the prior quarter (logged via the InstructionsLoaded observability hooks) and checks that hard-block invocations were genuine and soft warnings were either fixed or acknowledged. A hook with a >5% false-positive rate is downgraded to soft warning or removed.

The same doctrine is documented in `validate_naming.py`'s docstring; keep the two references in sync. The philosophy is project-wide, not per-hook.

When you propose a new check, classify it against this bar before choosing the exit code. If you're not sure whether the pattern is FP-near-zero, default to soft warning — promoting later is easy; demoting a hard block after Claude has been trained to ignore it is much harder.

## The hook ceiling

The hook count is technical debt. Every hook has to be maintained as the codebase changes (the layer-direction hook needs updates when layer rules change; the naming hook needs updates when conventions evolve; the observability hooks need updates if Claude Code's hook protocol changes; etc.). The current count of twelve is at the upper end of what's sustainable for a solo-maintainer project. Adding a thirteenth hook should require articulating why the existing twelve aren't enough — and should ideally pair with retiring one that's no longer earning its keep.

The twelve break down as: three deterministic enforcement (block_protected_paths, validate_layer_direction, validate_naming), two informational (pre_commit_summary, validate_commit_message), one bash-safety (block_dangerous_bash), one Stop-time check (stop_validate), one SessionStart user-facing reminder (session_start), two observability (log_session_start, log_post_tool_use), one PostToolUse formatter (post_edit_format), and one pre-commit gate (iwyu_check_on_commit).

Specific anti-patterns to avoid:

- **A hook that duplicates a skill.** The hook fires deterministically; the skill is suggestive. If both fire, the skill's effort is wasted. Pick one.
- **A hook whose rule changes frequently.** Frequent edits to a hook's regex or path list means the rule has soft edges that the regex can't capture. Convert to a skill (where Claude can apply judgment) or remove.
- **A hook with a long allowlist of exceptions.** Exceptions are how false positives leak back in. If the hook needs more than a handful of exceptions, the rule isn't actually deterministic.

## Maintenance lifecycle

### Adding a new hook

A new hook earns its place when:

1. The rule is deterministic — every match is a genuine violation.
2. The rule needs to be enforced rather than suggested. A skill or AGENTS.md note isn't enough.
3. CI doesn't catch it tightly enough. The feedback-latency improvement justifies the maintenance cost.

Procedure: write the script (Python or bash, both supported). Convention: the script reads the JSON payload from stdin, exits 0 to allow, exits 2 to block. Print clear error messages to stderr — Claude reads stderr and uses it to explain the violation to the user. Wire the hook into `settings.json` under the appropriate event matcher (`PreToolUse` for blocking, `PostToolUse`/`Stop`/`SessionStart` for informational). Add to this README's "What's here" section. Write at least a smoke test (for non-trivial logic, multiple cases) — the commit-message validator was 28-case smoke-tested when written.

The hook event lifecycle is documented in Anthropic's hooks reference.[^hooks-docs] The most useful events are `PreToolUse` (can block; happens before the tool runs), `PostToolUse` (informational; after the tool runs), `SessionStart`/`Stop` (session lifecycle), `UserPromptSubmit` (can modify or block prompts), `PreCompact` (before context compaction).

[^hooks-docs]: <https://code.claude.com/docs/en/hooks>

### Deleting a hook

Delete a hook when:

- It has not blocked anything in the last quarter (it's not catching violations because the rule is internalized, or because the violation never happens).
- Its rule has been superseded by a clang-tidy check, a CMake check, or another deterministic enforcement.
- It has accumulated false positives or exceptions to the point where you regularly bypass it.

Procedure: remove the script, remove the entry from `settings.json`, update this README. Git history preserves the script if you ever need it back.

### Refactoring a hook

Three patterns are common:

**Tightening the regex.** A hook that fires false positives needs its match logic constrained. The fix is usually a more specific pattern (e.g., a path prefix rather than a substring, or a word-boundary regex rather than a substring). Resist the urge to add an allowlist of exceptions; that's a smell that the rule isn't deterministic.

**Splitting one hook into two.** When a hook's logic has grown to handle two distinct rules, splitting often makes both clearer. The original `validate_layer_direction.py` was once `validate_includes.py` covering both layer direction and unrelated include hygiene; splitting into separate hooks made each one easier to reason about.

**Promoting a soft-warn to a hard-block.** When a soft-warn check fires and is consistently right, promote to hard-block. The validate_commit_message hook's "body required for substantive changes" check went through this path — it started as soft-warn, was right enough often enough that hard-blocking became safe.

**Adding informational output.** Many hooks benefit from richer stderr output without changing their block/allow logic. The pre_commit_summary hook is informational-only; some validation hooks could include more diagnostic context to help Claude explain violations to the user.

### Reviewing the hook set

Quarterly. Walk all hooks:

- Has this hook fired in the last quarter? Did it catch real violations or false positives?
- Is the rule still deterministic? Have any judgment calls crept into the regex?
- Does the error message still make sense? Does Claude know how to explain it to the user?
- Could this hook be deleted (CI now catches it; clang-tidy enforces it; the rule was internalized)?

The grounding for what each hook checks should match the actual source. The layer-direction hook's chain is grounded in CMakeLists.txt; if the layer chain changes, the hook needs updating. The naming hook's rules are grounded in `.clang-tidy`; if those change, the hook needs updating. Periodic re-grounding catches drift.

### Retiring a hook

Retire a hook that has been on the deletion candidate list for two quarters running. Hooks are technical debt; the bias should be toward fewer hooks rather than more.

The exception is the pre_commit_summary hook and similar informational hooks — they don't enforce anything, so they don't carry the false-positive risk. Informational hooks are cheap to keep around if they're useful at all.

## Cost model

Hooks cost zero context tokens. They run outside the model's context window entirely. The cost they do have is maintenance — every hook is code that has to be kept in sync with the rules it enforces. See `../cost-model.md` for how this fits into the overall mechanism cost picture.

## Recent changes (this directory)

This section records changes to the hook set — additions, deletions,
and behaviour changes (new rules enforced, thresholds adjusted, exit
codes changed). Cosmetic edits do not belong here. Bundle-wide
reorganizations are recorded in the top-level `README.md` instead.
Entries accumulate from production cutover onward.
