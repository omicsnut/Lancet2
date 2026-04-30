---
name: semantic-audit
description: Use when the user asks for a comment audit, a stale-docs sweep, a rigorous documentation sync, or to check whether comments and docs still say true things after a refactor. Use proactively on a quarterly cadence over one namespace at a time. Five-pass methodology: structural compliance, comment quality, doc accuracy, cross-file sync, terminology. Read-only walk of source, tests, and docs. Returns a backlog of findings ranked by severity. NOT for pre-merge review of a specific diff — that is fresh-reviewer's job.
allowed-tools: Read, Glob, Grep
---

# Semantic audit

A semantic audit is a file-by-file, line-by-line read of a target subdirectory's source, tests, and documentation to verify that comments and docs say true things, explain the right things, use the right words, are in the right place, and stay synchronized across code and docs. The methodology is derived from `docs_dev/style/sync_and_verification.md` (the synchronization rule and the pre-commit verification checklist).

## When to use this skill, and when not to

Use this skill when the user asks for a "rigorous doc check and sync", a "comment audit", a "stale doc sweep", or any variant of those phrases. Use it after a refactor that renamed a metric, a FORMAT field, or a CLI flag, to catch the resulting drift before the next release. Use it on a quarterly cadence over one namespace at a time as a maintenance practice.

Do NOT use this skill for pre-merge review of a specific diff — that is the `fresh-reviewer` subagent's job. The two are complementary: `fresh-reviewer` asks "did this diff introduce a new violation?" with priority on correctness, threading, and layer-direction; `semantic-audit` asks "what's the existing comment-and-docs backlog across this subdirectory?" with priority on comment quality and cross-file sync. A pre-existing filler word in code the diff doesn't touch is invisible to `fresh-reviewer` (correctly — it isn't what the writer changed) but is exactly what `semantic-audit` is for. For a large refactor PR, the right workflow is both: `fresh-reviewer` on the diff first, then this skill on the most-touched subdirectory second.

## Procedure

The audit is a single read-through per file in which all five passes happen simultaneously. Do not split the audit into five sequential reads of the same file; a skilled auditor catches all five categories of finding in one pass.

Process namespaces in dependency order: `src/lancet/base/` → `src/lancet/hts/` → `src/lancet/cbdg/` → `src/lancet/caller/` → `src/lancet/core/` → `src/lancet/cli/`. Process tests and benchmarks (only Pass 2 applies — they do not need algorithmic commentary). Then process `docs/guides/*.md` and `docs/reference.md` (only Passes 3 and 5 apply).

### Pass 1 — Structural compliance

For each `.h` and `.cpp` file, verify:

- **Include ordering**: the four-tier convention (main header → project → third-party → C++ stdlib → C stdlib) is respected, with blank lines between tiers.
- **Include style**: angle brackets only for C/C++ stdlib; quoted includes for everything else (Abseil, spdlog, htslib, minimap2, spoa, moodycamel, Catch2). Clang-format does not catch this.
- **`extern "C"` wrapping**: C-only third-party headers (htslib, minimap2) are wrapped with a brief comment noting they are POSIX/C headers.
- **Member layout**: structs and classes declare members in descending alignment (8B → 4B → 2B → 1B) with `// 8B`, `// 4B`, `// 2B`, `// 1B` size annotations. Constructor initializer lists match declaration order.
- **NOLINT discipline**: every suppression is `NOLINTNEXTLINE(check-name)` or `NOLINTBEGIN/END(check-name)` — never bare, never inline, always with a rationale.

### Pass 2 — Comment quality

For each comment block, ask four questions:

**Does it explain why, not what?** The code shows what; the comment must explain why this approach was chosen and what would break if it were done differently. `++mAligned; // increment aligned count` and `// Track the minimum base quality at an isolated read position natively.` are both violations: the first restates the code; the second describes the code's effect without rationale. The replacement is `// Record the weakest base quality in the variant region — the confidence in an indel observation is bounded by the least confident base (weakest-link).`

**Is every claim factually correct?** Read the code the comment describes. Comments commonly drift in these ways: comment says "returns X" but code returns Y after a refactor; comment says "uses parameter Z" but Z was renamed or removed; comment references a VCF FORMAT field name that was changed; comment says "O(N)" but the actual complexity is O(N log N) due to a sort.

**Is the math correct and annotated?** For any formula, verify each variable corresponds to a real program variable, check that representative values are correct (e.g. `Q=30 → 0.999`), and confirm the formula produces the stated result.

**Are invariants and boundary conditions stated?** For any function with non-obvious preconditions, the coordinate system (0-based vs 1-based, absolute vs relative), edge cases (empty input, zero coverage, single-element groups), and sentinel return values (`nullopt`, `0.0`, `SIZE_MAX`) should be explicit.

### Pass 3 — Language compliance

This pass is not a checklist exercise. The rules below catch common violations, but new filler words appear in future code that no table can anticipate. The core skill is **semantic reading of English**: reading each sentence and asking *does every word carry technical meaning, or is something just making the sentence sound more authoritative without adding information?*

**The deletion test is the primary tool.** Remove the suspect word or phrase. Re-read the sentence. If the technical meaning is unchanged, the word was filler. If the sentence loses a real distinction (compile-time vs runtime, contrasting with an alternative, quantifying a tradeoff), the word earns its place.

Filler words previously caught: `natively` (always delete; means "the code does what the code does"), `organically` (always delete; code is deterministic, not organic), `seamlessly` (always delete; never true for code readers debugging an issue), `fundamentally` (acceptable when describing a structural or topological difference; delete when it means "very" or "importantly"), `inherently` (acceptable for a mathematical property; delete when it means "obviously"), `intrinsically` (same test as `inherently`).

Filler adverbs previously caught: `mathematically` (keep when describing a proof or derivation step; delete when restating that code does arithmetic), `explicitly` (keep when contrasting with implicit behavior; delete when emphasizing that code does what code does), `structurally` (keep when describing data-structure topology; delete when modifying a verb for emphasis), `aggressively` (keep when quantifying a design tradeoff with a number; delete as a vague intensifier), `dynamically` (keep when distinguishing runtime from compile-time; delete when it means "it happens at runtime", which is obvious), `deliberately` (keep when contrasting with accidental behavior; delete when redundant with a `CRITICAL:` note), `unconditionally` (keep for guard clauses; delete when redundant with "always"), `affirmatively` (never keep; replace with specifics), `maliciously` (never keep; anthropomorphizes software).

**Watch for adverb clusters.** Two or more adverbs modifying the same clause is a strong signal that at least one is filler. "Aggressively and structurally prevents maliciously winning" — three adverbs, all filler.

**Watch for anthropomorphization.** Code does not "want", "try", "refuse", "fight", or "maliciously win". If a comment attributes intent to software, rephrase to describe the mechanism.

**Watch for hedging.** "Essentially", "basically", "more or less", "in a sense" — if the claim is true, state it directly. If it is an approximation, quantify the approximation.

**Unexplained jargon.** Apply the biologist test: would a biologist or clinician have to stop and look it up? If yes, either replace with plain language ("independent" not "orthogonal") or define on first occurrence with a brief parenthetical (`Bessel's correction (divides by n−1 instead of n to avoid underestimating spread from a sample)`).

### Pass 4 — Header vs implementation split

Headers (`.h`) document the public API: what each method does, what types it accepts and returns, what invariants it maintains. They should NOT contain algorithmic walkthroughs. Implementation files (`.cpp`) document the algorithm: step-by-step logic, formula derivations, ASCII diagrams, coordinate-system handling. They should NOT duplicate the API contract. Massive comment blocks in headers explaining algorithms are a Tier 2 finding; the explanation belongs in the `.cpp`.

### Pass 5 — Cross-file synchronization

For each concept the current file documents, verify:

- If the file mentions a VCF FORMAT/INFO/FILTER field, the field name matches `variant_call.h`, the registration in `vcf_header_builder.cpp`, and the user-facing description in `docs/guides/vcf_output.md`.
- If the file documents a metric (e.g., NPBQ, SB, RPCD), the formula and description match the implementation in `variant_support.cpp` and the corresponding web doc in `docs/guides/`.
- If the file mentions another file, function, or class by name, that name still exists.

The patterns that go stale most reliably: architecture diagram comments referencing FORMAT field lists, metric descriptions in tangentially related functions, `ReadEvidence` member comments referencing FORMAT field names, operating-mode tables in `vcf_output.md` referencing QUAL computation methods, cross-references in one doc page to a metric described in another.

## Findings classification

Every finding is classified into exactly one tier.

**Tier 1 — Convention violation, must fix.** A codified rule from `docs_dev/style/` (entry point: `docs_dev/style/README.md`) or `.clang-tidy` is broken: angle-bracket include on a third-party header, banned filler word with no technical meaning, members declared out of alignment order, inline `NOLINT`, etc.

**Tier 2 — Comment quality issue, should fix.** The comment is technically not rule-breaking but fails to serve its purpose: comment restates the code instead of explaining why, non-obvious logic has no rationale, formula in comment does not match the code, missing boundary-condition documentation, algorithmic walkthrough in a header rather than the `.cpp`.

**Tier 3 — Documentation sync gap.** A web documentation page or cross-file reference is stale or inconsistent: `docs/guides/architecture.md` describes a metric that was renamed in code, a VCF FORMAT field comment in `variant_call.h` uses a different formula than `variant_support.cpp`, a doc page uses filler jargon that the code comments have already been cleaned of.

## Report template

Structure the report per-tier, per-file. Each finding has:

```markdown
### [path/to/file.cpp:NNN] — <Rule Category>

**Current**:
<verbatim offending comment text>

**Rule**: <quote the specific convention rule being violated>

**Fix**:
<the exact replacement text>
```

End the report with an explicit list of files that passed the audit with no findings. The pass list proves the audit was exhaustive — every file in scope appears either in findings or in the pass list — and identifies gold-standard files that new contributors should study.

## Anti-patterns to avoid during audit

**Grep-only scanning.** Grep finds words, not meaning. `grep "natively"` catches the word but not "the DM model handles multi-allelic sites" (a sentence that restates the class docstring 20 lines above). Read the code.

**Checking comments in isolation.** A comment is correct only in the context of the code it describes. Reading a comment without reading the adjacent code cannot detect factual inaccuracy.

**Applying rules mechanically.** "fundamentally" is sometimes correct (structural topology) and sometimes filler (emphasis). "mathematically" is correct in a derivation step, filler when restating arithmetic. Every instance requires judgment.

**Reporting style preferences as violations.** The audit checks rules and factual accuracy, not personal writing style. "Subtracting `local_raw_score` carves a hole" is fine — it is a metaphor that aids understanding. "Subtracting `local_raw_score` natively carves an exact algebraic hole" is a violation — three filler words.

**Skipping exemplary files.** List them. The pass list proves completeness and identifies gold-standard files that new developers should study.

## Audit cadence

The audit fits into the project at three frequencies:

**Per-PR (writer self-audit, then reviewer).** The developer who writes the code applies Passes 1 through 4 to their own changes before opening the PR. The reviewer re-applies them during review — but the reviewer's pass is naturally diff-scoped, so the writer's self-audit is the more thorough application of this skill at PR time.

**Quarterly (full-scope sweep).** A full-namespace audit (one namespace at a time, in dependency order) catches drift that no per-PR audit can. The findings report is committed as a tracking artifact under `notes/audits/YYYY-QN/`.

**Post-refactor (focused sync).** Any rename, restructure, or metric change triggers a focused Pass-5 sync audit on every file that referenced the changed concept. This is the audit that most reliably catches stale comments, because the writer of the rename has the field name fresh in mind and can grep for it exhaustively.
