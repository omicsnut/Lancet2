# Synchronization and Verification

This document holds two things: the cross-document synchronization rule that keeps the codebase, code comments, dev docs, and website docs consistent with each other, and the pre-commit verification checklist that catches violations before they merge.

The sync rule applies project-wide — every audience and every document type is bound by it, except `docs_dev/investigations/` (which is immutable by design and intentionally not maintained as the codebase evolves).

## The synchronization rule

Every sentence in code comments, every sentence in website documentation, and every sentence in `docs_dev/` (excluding `investigations/`) must stay in sync with the actual logic in the codebase. Stale comments and stale docs are worse than no documentation: they confidently mislead.

When any refactor changes a field name, metric computation, or behavioural semantic, perform an exhaustive codebase-wide audit by grepping for the old names across:

- `src/**/*.h` and `src/**/*.cpp` — code comments and docstrings
- `docs/**/*.md` — website documentation and tables
- `docs_dev/**/*.md` — internal developer documentation
- `tests/**/*.cpp` — test names and assertions that may reference old field names
- VCF header definitions in `src/lancet/cli/vcf_header_builder.cpp`
- FORMAT string comments in `src/lancet/caller/variant_call.cpp` and `variant_call.h`

The rule is mechanical: if the new code emits or computes something different, the comments and docs that describe it must change in the same commit. Splitting the code change and the doc change across commits is how documentation goes stale — the code commit lands and the doc commit gets forgotten.

### Common patterns that go stale

These are the patterns most likely to drift if you don't audit explicitly:

- Architecture diagram comments referencing VCF FORMAT field lists.
- Metric descriptions in tangentially related functions (e.g., a comment about `mScore` in a function that doesn't compute it).
- `ReadEvidence` and similar struct member comments referencing VCF field names.
- Cross-references in one doc page to a metric described in another, when the metric is renamed.
- `mkdocs.yml` nav titles when a page is renamed but the nav entry is not.
- Pixi task names in CLAUDE.md or in `docs_dev/workflows/` runbooks when `pixi.toml` task names change.
- `.clang-tidy` complexity thresholds that drift away from the table in `cpp_style.md`.

## Pre-commit verification checklist

Run through this checklist before every commit that touches documentation or code that documentation describes. The list is short by design — items here are the ones that catch real, recurring failures.

### Website documentation

If you edited a file in `docs/`:

- [ ] Page opens with a 2–3 sentence orientation paragraph immediately under `# Title` (no `## Introduction` heading).
- [ ] Heading hierarchy is `## → ### → ####`, no skipped levels.
- [ ] Inline code uses backticks; key terms on first definition use bold (not both).
- [ ] Cross-references at the end (`* **Read more:**` or `* **CLI reference:**`).
- [ ] No unexplained jargon — apply the "would a clinician have to look this up?" test.
- [ ] No hedging language ("should", "is expected to") without an explicit reason.
- [ ] Tables left-aligned with `:------`.
- [ ] If new admonitions added: at most one or two per page, used for genuine warnings/notes/tips.
- [ ] Mathematical notation uses Unicode (×, −, ≥, →) rather than ASCII approximations.

See `website_docs.md` for the full set of rules.

### Code comments

If you edited code comments in `src/`:

- [ ] Every named constant has a comment explaining what it represents, why this value, and what changes if modified.
- [ ] Comment blocks above ~15 lines have been considered for migration to `docs_dev/subsystems/` or `docs/guides/` with a cross-reference.
- [ ] No restated code (`++count; // increment count`).
- [ ] No filler adverbs (drop "strictly", "fundamentally", "essentially").
- [ ] Member-variable size annotations (`// 8B`, `// 24B`) are present and correct on declaration lines.
- [ ] ASCII art and visual tables are wrapped in `// clang-format off` / `// clang-format on`.
- [ ] TODO/FIXME comments include what needs to change, why, and what blocks them — no bare `// TODO: fix this`.

See `code_comments.md` for the full set of rules.

### Code style and linting

If you edited C++ in `src/`:

- [ ] `pixi run fmt-check` is clean.
- [ ] `pixi run lint-check` is clean (or every NOLINT suppression is scoped, named, and has a rationale comment).
- [ ] Struct/class members are in 8B → 4B → 2B → 1B order, with constructor initializer lists matching declaration order.
- [ ] Quote-vs-angle-bracket include rule respected: stdlib uses `<>`, everything else uses `""`.
- [ ] If a function exceeds a complexity threshold (cognitive=35, statements=150, lines=200, params=6, nesting=4), it has been decomposed or has a justified scoped NOLINT.
- [ ] Manual loops have been considered for replacement with `<algorithm>`, `<ranges>`, or Abseil utilities — but only where it improves clarity.
- [ ] New NOLINT suppressions are scoped (`NOLINTNEXTLINE` or `NOLINTBEGIN`/`NOLINTEND`), named, with rationale, and disclosed in the commit message.

See `cpp_style.md` for the full set of rules.

### Cross-document synchronization

If your change touches a field name, metric, threshold, or behavioural semantic:

- [ ] Audited `src/**/*.{h,cpp}` for occurrences of the old name in comments and code.
- [ ] Audited `docs/**/*.md` and `docs_dev/**/*.md` (excluding `investigations/`) for occurrences in prose and tables.
- [ ] Audited `tests/**/*.cpp` for assertions that reference the old name.
- [ ] If the change touches VCF emission: VCF header line in `vcf_header_builder.cpp`, FORMAT comments in `variant_call.{h,cpp}`, the `vcf_output.md` annotations guide, and any related `docs/guides/` page are all updated.
- [ ] If the change touches a CLI flag: `cli_interface.cpp`, `docs/reference.md`, and any guide that mentions the flag are updated.
- [ ] If the change touches a pipeline stage name: code comments + architecture diagrams in source + the relevant `docs/guides/` page + any subsystem deep-dive in `docs_dev/subsystems/` are updated.

### Architecture and investigation documents

If you added a document to `architecture/` or `investigations/`:

- [ ] You read the README of the relevant subdirectory before drafting.
- [ ] The document follows the structure prescribed there.
- [ ] For architecture documents: a "last reviewed" date header is present.
- [ ] For investigation documents: a "captured-on" date header is present, the document explicitly states its immutability, and you intend to leave it unchanged after merge.
- [ ] The document is cross-referenced from any pre-existing document where readers will benefit from finding it.

## When violations are found post-merge

If you discover a sync violation in a document that has already merged — a comment that no longer matches the code, a stale claim in a guide — fix it in a dedicated `chore:` commit with a subject like `chore: sync <metric/field/etc> documentation across <files>`. Do not roll the fix into an unrelated change; sync fixes should be findable in the git history.

The exception, again, is `investigations/`. If an investigation contains a claim that is no longer true, do not edit the investigation. Write a new one that supersedes it, with a header that links back to the original.
