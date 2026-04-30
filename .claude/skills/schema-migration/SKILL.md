---
name: schema-migration
description: Use when the user proposes adding, removing, renaming, changing the cardinality of, or silently changing the semantics of a VCF FORMAT, INFO, or FILTER field, or when changing a user-facing CLI flag. Trigger on "add a FORMAT field", "rename SB to STRBIAS", "change AD from R to A", "deprecate this field", "the field still exists but means something different now", "rename the --tumor flag", "change the default for --num-threads". Walks an interview before applying any edit, classifies the change against the five-operation matrix (add / rename / cardinality-change / remove / silent-semantic-change), produces the coordinated edits across vcf_header_builder.cpp + reader/writer code + tests + docs, contributes a structured paragraph to the commit body, and invokes the vcf-validator subagent before finalizing. Does NOT cover internal struct fields or non-VCF data formats — those go through normal review.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Schema migration on Lancet2

Lancet2's VCF schema is a stable interface that downstream pipelines depend on. Most fields are read by automated tooling that doesn't tolerate silent renames, cardinality changes, or semantic shifts. This skill formalizes the procedure so that schema changes are explicit, coordinated across all the places that need to know, and reviewed by the vcf-validator subagent before merge.

The same discipline applies to user-facing CLI flags. Renaming `--tumor` or changing the default value of `--num-threads` is a schema change for the CLI surface; downstream wrapper scripts break in the same way as VCF parsers.

## The five operations

Every schema change falls into one of these:

1. **Add** — new FORMAT/INFO/FILTER field or new CLI flag. Lowest-risk; downstream consumers either start using the new field or keep ignoring it. Still needs documentation.
2. **Rename** — an existing field gets a new name with the same semantics. High-risk: every downstream parser breaks. Requires deprecation period and dual-write if the field is widely consumed.
3. **Cardinality change** — a FORMAT field's `Number=` attribute changes (e.g., `Number=R` to `Number=A`, `Number=1` to `Number=A`). Subtle: the new file parses successfully but downstream code that indexed by allele count gets wrong values. Highest-risk class because the breakage is silent.
4. **Remove** — a FORMAT/INFO/FILTER definition or CLI flag is deleted. Downstream code that referenced it now fails or gets default values. Requires deprecation period.
5. **Silent semantic change** — the field's name and cardinality stay the same but the underlying meaning changes (e.g., `SB` previously was Phred-scaled strand bias, now is log odds ratio). The most dangerous because no parser breaks; downstream interpretations are silently wrong. ALWAYS rename the field instead of changing semantics — even if the rename is just `SB` → `SB2`.

## Step 1 — Interview before editing

Do NOT begin editing on first mention of a schema change. Ask the user:

1. **Which operation is this?** Show the five-operation list above and ask them to pick. If they say "we're just changing what SB means," walk through why that's operation 5 and propose a rename instead.
2. **Which field, and what is the change?** Get the exact field name, the current `##FORMAT`/`##INFO`/`##FILTER` line in `src/lancet/cli/vcf_header_builder.cpp` (the FORMAT defs occupy a contiguous block beginning around line 35; confirm the exact span by reading the file rather than assuming), and the proposed new line.
3. **What downstream consumers will be affected?** Names of pipelines, scripts, or notebooks that parse this field. The user may not know all of them, but document what they do know in the commit body.
4. **For renames and removes: what is the deprecation strategy?** Either dual-write for N releases, or a hard cutover with a CHANGELOG warning. Hard cutover is rare and should be justified.
5. **For cardinality changes and silent-semantic changes: are you sure?** These are the highest-risk classes. Push back if the user hasn't fully reasoned about downstream impact.

If any answer is unclear, stop and ask. Don't proceed with the edits until the operation type and the rationale are explicit.

## Step 2 — Find every coordinated location

A FORMAT/INFO/FILTER change touches at least these places:

- `src/lancet/cli/vcf_header_builder.cpp` — the canonical declaration site. The `FORMAT_STR_HEADER` raw string is the source of truth for FORMAT and INFO definitions in single-sample mode; a separate `CASE_CTRL_INFO_HDR_LINES` raw string declares the case/control INFO fields (SHARED, CTRL, CASE) that emit only in case-control mode. Read the file before editing to see the exact line spans.
- `src/lancet/caller/` — the code that populates the field. For a FORMAT field, this is typically the genotyper or the per-sample feature emitter.
- `tests/caller/` and `tests/cli/` — tests that consume or assert against the field.
- `docs/guides/vcf_output.md` — the user-facing schema documentation.
- `docs/PROBE_TRACKING.md` — if the field is used by the probe pipeline (FORMAT fields like AD, DP, RMQ are referenced there).
- `.chglog/CHANGELOG.tpl.md` and the upcoming CHANGELOG entry — schema changes are a release note category.

For a CLI flag change, the locations are:

- `src/lancet/cli/cli_interface.cpp` and the relevant `pipeline_runner.cpp` — flag definition and handling.
- `src/lancet/cli/argument_validator.h` — if the flag has constraints.
- `tests/cli/` — tests that exercise the flag.
- `docs/index.md`, `docs/guides/wgs_analysis.md`, `docs/guides/targeted_analysis.md` — user-facing CLI usage.

Use `pixi run` `Glob` and `Grep` to find every occurrence of the field/flag name across the codebase. Build a list before editing; never edit one location at a time.

## Step 3 — Propose the coordinated edit set

Before applying any change, present the user with the full list of edits:

```
Schema migration: rename FORMAT/SB to FORMAT/SBLOR (strand bias log odds ratio)
Operation: 2 (rename)
Deprecation: dual-write SB and SBLOR for one release; remove SB in next.

Edits to apply:
  1. src/lancet/cli/vcf_header_builder.cpp:45 — change ##FORMAT line; add SBLOR with same Description.
  2. src/lancet/caller/feature_emitter.cpp:NNN — emit both SB and SBLOR (dual-write).
  3. tests/caller/feature_emitter_test.cpp:NNN — add SBLOR assertion alongside SB.
  4. tests/cli/vcf_output_test.cpp — add SBLOR header check.
  5. docs/guides/vcf_output.md — document SBLOR; mark SB as deprecated.
  6. CHANGELOG entry: "FORMAT/SB renamed to FORMAT/SBLOR. Both fields emitted this release; SB will be removed in 3.0.0."

Confirm to proceed (yes/no/refine)?
```

Wait for explicit confirmation. If the user says "refine," accept their corrections and re-present.

## Step 4 — Apply the edits

Apply the edit set as a single coordinated change. Run `pixi run lint-check` and `pixi run test` after each major edit to surface compilation or test failures early. The protected-paths hook will not block any of these (none are on the protected list).

## Step 5 — Invoke the vcf-validator subagent

Before finalizing, hand the change to the vcf-validator subagent:

```
Use the vcf-validator subagent to review this schema migration.

Operation: <add/rename/cardinality-change/remove/silent-semantic-change>
Field: <FORMAT/INFO/FILTER name and direction>
Edits applied: <list from Step 3>

Verify the schema invariants in vcf-validator's frame and surface any
breakage (Number= cardinality, Description completeness, Type validity,
multi-sample consistency, missing-value semantics).
```

The subagent reads the actual VCF output (use a small `pixi run test-asan` invocation to produce one quickly, or run the pipeline against `LANCET_TEST_*_REGION_SMALL`) and checks every invariant. Address findings before merge.

## Step 6 — Contribute to the commit body

Schema changes get a structured paragraph in the commit body:

```
feat: rename caller FORMAT/SB to FORMAT/SBLOR

Operation: rename (operation 2 of the schema-migration matrix).
Deprecation: dual-write SB and SBLOR this release; SB removed in 3.0.0.
Rationale: SB historically was Phred-scaled strand bias but the
implementation has been log-odds-ratio for several releases. Renaming
makes the meaning self-documenting and prevents future readers from
assuming the original Phred semantics.

Downstream impact: any consumer indexing on FORMAT/SB will continue
to work this release. Consumers should switch to FORMAT/SBLOR before
3.0.0.

vcf-validator: clean.
```

The structured form makes the CHANGELOG entry mechanical to generate (`scripts/update_changelog.sh` via git-chglog) and gives downstream maintainers the impact summary they need.

## When NOT to use this skill

Do not use this skill for changes to internal struct fields, internal class APIs, or non-VCF data formats (e.g., the probe pipeline's `probe_results.tsv` columns or any private TSV format). Those go through normal code review. Do not use this skill for purely additive Description changes (typo fixes in an existing `##FORMAT` line) — those are spelling fixes, not schema changes.

## When the change is operation 5 (silent semantic change)

Strongly resist this. Walk through with the user why a rename is the correct response — even a small one like `SB` → `SBNEW`. Silent semantic changes are how downstream pipelines break in production six months later. If the user insists, document in the commit body that this was a deliberate operation-5 with full awareness of the downstream risk, and add a CHANGELOG warning at the top of the entry.
