---
name: documentation-sync
description: Use when the user changes code that has corresponding documentation, or changes documentation that has corresponding code, and the three layers (in-source comments, dev docs in docs_dev/, user docs in docs/) need to stay in sync. Trigger on "I just changed X, did I miss the docs?", "the architecture guide says Y but the code does Z", "update the VCF schema docs", "add this new flag to the wgs_analysis guide", "the docstring on Foo doesn't match the implementation", "we renamed FOO to BAR; what else needs updating?", "the probe_tracking subsystem doc is out of date". Walks the multi-direction sync (code ↔ docs_dev ↔ docs), surfaces drift, proposes coordinated edits, and gates the change behind an explicit consistency check. The authoritative design for the synchronization rule lives in docs_dev/style/sync_and_verification.md; this skill is the procedural playbook that applies that rule. Does NOT cover docs_dev/investigations/ (immutable by design — never updated when the codebase changes), the agent-memory files, the bundle's own .claude/ READMEs, or AGENTS.md.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Documentation sync on Lancet2

Lancet2 has three layers of documentation that must stay synchronized with the source code. The split is by audience:

1. **In-source documentation** — header comments, function-level docstrings, inline rationale in `.cpp` files. Audience: a contributor reading the source. Voice: "what does this code do, why is it designed this way, and what would break if it were done differently."

2. **Developer documentation under `docs_dev/`** — deep-dives on subsystems, runbooks for recurring workflows, architecture decision records, postmortems. Audience: a contributor working on the codebase but not currently reading the specific source file. Voice varies by subdirectory: explanatory in `subsystems/`, procedural in `workflows/`, decisional in `architecture/`, narrative in `investigations/`.

3. **User-facing documentation under `docs/`** — the rendered MkDocs site. Audience: someone running Lancet2, not contributing to it. Voice: task-oriented — "how do I run a tumor-normal pipeline," "what does this output field mean," "which CLI flag controls this behavior."

The synchronization rule is mechanical: **every claim about the codebase made in any of these three layers must remain consistent with the actual code as the codebase evolves.** When source code changes, the relevant comments, dev docs, and user docs must change in the same commit. Splitting the code change and the doc change across commits is how documentation goes stale — the code commit lands, the doc commit gets forgotten, and a future contributor reads documentation that confidently misleads.

The authoritative source for the rule's exact wording, the file globs to grep, and the pre-commit verification checklist is `docs_dev/style/sync_and_verification.md`. Read it before applying this skill the first time; it carries the project-wide rule with the project owner's exact language. This skill is the procedural playbook that applies that rule in a specific session.

The exception that crosses the rule is `docs_dev/investigations/`. Investigations capture a moment in time — what was understood, what was tried, what worked and what didn't. They are intentionally NOT maintained as the codebase evolves. If understanding later changes, write a new investigation that supersedes the old one rather than editing the old one in place. The exception is the entire reason `investigations/` is its own subdirectory.

## Step 1 — Identify the direction(s) of sync

Three possible starting directions:

- **Code-to-docs.** The user just changed code (added a function, renamed a CLI flag, changed a default, modified a VCF field, restructured a subsystem). Find every place in `docs_dev/` and `docs/` that describes the changed code; those places are now potentially stale. Walk forward in both layers.
- **docs_dev-to-code-and-user-docs.** The user updated a developer doc (e.g., fixed a probe-tracking subsystem description, landed an ADR that shifts a design contract). Find code comments that describe the same thing and user-facing docs that depend on it. Code-comment drift means a future contributor reading the source gets a different story than one reading `docs_dev/`. User-doc drift means the user-facing site disagrees with the dev docs.
- **user-docs-to-code-and-dev-docs.** The user updated the user-facing site (e.g., clarified a CLI flag's behavior, corrected a description of what a FORMAT field means). Find code comments that describe the flag's implementation and dev docs that describe its design history.

Most cases are code-to-docs. The other two directions usually start when someone catches a documentation error and propagates the fix back through the layers.

## Step 2 — Run the cross-layer audit

Use grep across all three layers for the renamed term, the deprecated workflow, the changed default, the new field. The audit pattern from `docs_dev/style/sync_and_verification.md`:

```bash
# Replace OLD_NAME with the term being renamed/changed
grep -rn "OLD_NAME" \
    src/**/*.h src/**/*.cpp \
    docs/**/*.md \
    docs_dev/**/*.md \
    tests/**/*.cpp
```

The list above covers the typical sync surface. Specific high-value paths to always include for VCF-schema-adjacent changes:

- **`src/lancet/cli/vcf_header_builder.cpp`** — VCF FORMAT/INFO/FILTER definitions (the canonical declaration site).
- **`src/lancet/caller/variant_call.{h,cpp}`** — FORMAT field-name comments and per-record emission code.

For CLI flag changes, the surface narrows: `src/lancet/cli/cli_interface.cpp`, the `pipeline_runner.cpp` consumer, `tests/cli/`, the user-facing flag references in `docs/index.md` plus the relevant `docs/guides/`, and any `docs_dev/architecture/` ADR that decided the flag's semantics.

For subsystem changes (a new pipeline stage, a refactored component): the implementing layer (`src/lancet/<layer>/`), `docs_dev/subsystems/<topic>.md` if it exists, `docs_dev/workflows/<topic>.md` if applicable, and `docs/guides/<topic>.md`.

For algorithm changes (a new heuristic, a tuned constant, a different statistical model): the implementing source files, possibly `docs_dev/architecture/` if the change has cross-cutting consequences, and the user-facing `docs/guides/<topic>.md` describing what the user observes.

## Step 3 — Surface every divergence before editing

Before applying any change, list every location found and present them to the user:

```
Code change: renamed Genotyper::EvaluateAllele to Genotyper::ScoreAllele.

Drift found:
  Code comments:
    - src/lancet/caller/feature_emitter.h:27 — references EvaluateAllele as score source
    - src/lancet/caller/genotyper.cpp:84 — comment block describes EvaluateAllele step
  docs_dev:
    - docs_dev/subsystems/probe_tracking.md:142 — refers to "EvaluateAllele step" in cascade narrative
    - docs_dev/architecture/0003-genotyping-redesign.md:51 — cites EvaluateAllele in the design rationale
  docs:
    - docs/guides/variant_discovery_genotyping.md:42 — refers to "EvaluateAllele step"
    - docs/guides/scoring_somatic_variants.md:18 — references EvaluateAllele in pseudocode

Proposed updates: rename in all six locations.

Confirm to proceed (yes/no/refine)?
```

Wait for confirmation. The same "show before applying" discipline as `schema-migration` — the user catches mistakes (or scope misjudgments) before edits propagate. If the user says "refine," accept their corrections and re-present.

## Step 4 — Apply the coordinated edit set

A single coordinated change set across all three layers wherever drift was found. Run the relevant verification after each major edit:

- For `docs/` changes: `pixi run docs-build` (mkdocs `--strict` mode; broken cross-links fail the build).
- For `docs_dev/` changes: there is no separate build step — `docs_dev/` is not rendered. The verification is grep-driven: after the edit, re-run the Step 2 audit pattern with the OLD name and confirm zero matches outside `docs_dev/investigations/`.
- For code comment edits: the next compilation pass surfaces any syntactic problem (highly unlikely for comment-only edits), but a comment-only edit cannot break compilation. Review is the only verification.

## Step 5 — Verify mkdocs strict-build

```bash
pixi run docs-build
```

This is `mkdocs build --strict`. A clean exit confirms `docs/` is internally consistent (no broken cross-links, no missing images, no malformed front matter). It does NOT confirm `docs/` is consistent with `docs_dev/` or with code — that's what the cross-layer audit above is for.

## Step 6 — Note the doc updates in the commit

The Lancet2 chglog filter accepts only `feat`, `fix`, `perf`, and `chore` as commit types; `docs:`, `refactor:`, `build:`, etc. parse but are silently dropped from `CHANGELOG.md`. The chglog header pattern also does not support scopes — `refactor(caller): ...` and `feat(scripts): ...` are dropped for the same reason.

For documentation-only commits, use `chore:`. For mixed code+docs changes, the primary code commit type (`feat:`, `fix:`, `perf:`, `chore:`) wins; mention the doc updates in the body:

```
chore: rename EvaluateAllele to ScoreAllele

Function name now matches the operation (computing a score, not
evaluating an allele's existence). No behavior change.

Doc updates:
  Code comments: feature_emitter.h, genotyper.cpp
  docs_dev:      subsystems/probe_tracking.md, architecture/0003-...
  docs:          guides/variant_discovery_genotyping.md,
                 guides/scoring_somatic_variants.md
```

If the doc changes are the primary commit (a guide rewrite, a new subsystem deepdive, a new ADR), use `chore:` and mention any code-comment updates in the body.

## When NOT to use this skill

Do not use this skill for:

- **`docs_dev/investigations/`.** Investigations are immutable by design. They are NOT updated when the codebase changes. If understanding later changes, write a new investigation that supersedes the old one — do NOT edit the old one in place. This is the only path-scoped exception to the synchronization rule.
- **AGENTS.md, CLAUDE.md, or `.claude/` READMEs.** Those are bundle-internal documentation with their own update conventions; changes go through the `/audit-bundle` slash command's quarterly review.
- **Agent-memory files.** Those are append-only logs the agents themselves maintain. Manual edits should be rare and don't follow this skill's procedure.
- **CHANGELOG.md.** Generated by `pixi run update-changelog` from commit messages; the right place to influence it is the commit body.
- **Vendored dependency documentation under `cmake-build-*/_deps/`.** Read-only; never edited.

## When the change is large enough to need a new doc

If the code change introduces a new subsystem, a new substantive feature, or fundamentally changes a workflow, the right move is often a new document rather than incremental edits to existing ones:

- New subsystem → new `docs_dev/subsystems/<topic>.md` (and possibly a paired user-facing guide at `docs/guides/<topic>.md`).
- New recurring workflow → new `docs_dev/workflows/<topic>.md`.
- Substantive design decision → new `docs_dev/architecture/NNNN-<slug>.md` (ADR). Use the `/arch-decision-record` slash command, which interviews the user and produces a draft following `docs_dev/architecture/README.md`'s writing standards.
- Postmortem or debug archaeology → new `docs_dev/investigations/<date>-<slug>.md`. Use the `/investigation` slash command.

Discuss the structure with the user before drafting; a new document is a longer-lived artifact than an inline comment and benefits from intentional placement.
