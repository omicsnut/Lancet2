# /audit-probe-pipeline — Drift check for the probe tracking infrastructure

Walk `docs_dev/subsystems/probe_tracking.md` against the actual probe-tracking source and scripts, report any drift in stage names, file paths, ownership chain, API surface, or the 27-stage cascade. The probe pipeline has many moving parts (8 C++ files, 2 scripts, 1 doc) and the bundle's claims about it must stay correct against all of them.

This is the third member of the bundle's drift-command set, alongside `/sync-cost-model` (drifts against Anthropic) and `/audit-bundle` + `/audit-vcf-schema` (drift against project source). Probe tracking gets its own command rather than being folded into `/audit-bundle` because the surface is wide enough that bundling would dilute focus.

## When to run

After any change to the probe-tracking source: `src/lancet/cbdg/probe_index.{h,cpp}`, `src/lancet/cbdg/probe_tracker.{h,cpp}`, `src/lancet/cbdg/probe_results_writer.{h,cpp}`, `src/lancet/core/probe_diagnostics.{h,cpp}`. After any change to the scripts: `scripts/truth_concordance.py`, `scripts/analyze_probe_results.py`. After any edit to `docs_dev/subsystems/probe_tracking.md`.

Also run as a quarterly maintenance check, paired with `/audit-bundle` and `/sync-cost-model`.

## Procedure

The audit walks five pairings. For each, read the source-of-truth files, read the bundle's claims, diff them, record findings with file:line on both sides. Do not fix inline; the report is the deliverable.

### Pairing 1 — Stage cascade (27 levels)

**Source of truth:** the `STAGE_ORDER` constant in `scripts/analyze_probe_results.py`.

Read the file and extract every entry from `STAGE_ORDER`. The expected categories (in attribution-priority order, deepest first):

1. Genotyper (3): `geno_reads_reassigned`, `geno_zero_alt_reads`, `geno_no_overlap`
2. Variant extraction / MSA (3): `msa_not_extracted`, `msa_shifted`, `msa_subsumed`
3. Path finding (2): `bfs_exhausted`, `no_path`
4. Pruning (6): `pruned_at_tips`, `pruned_at_compress2`, `pruned_at_lowcov2`, `pruned_at_compress1`, `pruned_at_lowcov1`, `pruned_at_build`
5. Graph construction (5): `graph_has_cycle`, `graph_too_complex`, `no_anchor`, `short_anchor`, `variant_in_anchor`
6. Not processed (7): `not_processed` plus 6 sub-stages (`not_processed:ref_all_n`, `not_processed:ref_repeat`, `not_processed:inactive`, `not_processed:low_coverage`, `not_processed:no_alt_haplotype`, `not_processed:other_variant_called`)
7. Survived (1): `survived`

Categories sum to 3+3+2+6+5+7+1 = 27 stages.

**Bundle locations to check:**

- `.claude/agents/probe-interpreter.md` — the stage→source mapping table.
- `.claude/skills/probe-tracking/SKILL.md` — any stage references.
- `docs_dev/subsystems/probe_tracking.md` — the canonical doc itself.

Report any stage name in any of these that does not appear in `STAGE_ORDER`, plus any stage in `STAGE_ORDER` that is not mentioned in `probe-interpreter.md` (a missing stage means the agent will not know what code is responsible for it).

### Pairing 2 — `not_processed` sub-stages

**Source of truth:** the StatusCode enum referenced in `src/lancet/core/variant_builder.cpp` plus the mapping table in `analyze_probe_results.py` (the section that maps Window status → not_processed sub-stage).

Expected sub-stages:

| Sub-stage | Window status |
|:----------|:--------------|
| `not_processed:ref_all_n` | `SKIPPED_NONLY_REF_BASES` |
| `not_processed:ref_repeat` | `SKIPPED_REF_REPEAT_SEEN` |
| `not_processed:inactive` | `SKIPPED_INACTIVE_REGION` |
| `not_processed:low_coverage` | `SKIPPED_ANCHOR_COVERAGE` |
| `not_processed:no_alt_haplotype` | `SKIPPED_NOASM_HAPLOTYPE` |
| `not_processed:other_variant_called` | `FOUND_GENOTYPED_VARIANT` |

Verify each sub-stage exists in both source files and that the mapping is consistent. New StatusCode values not surfaced as a sub-stage are findings.

**Bundle locations to check:** `.claude/skills/probe-tracking/SKILL.md`, `.claude/agents/probe-interpreter.md`, `docs_dev/subsystems/probe_tracking.md`.

### Pairing 3 — Probe API surface

**Source of truth:** the `Probe*` method declarations in `src/lancet/cbdg/probe_tracker.h` and the call sites in `src/lancet/cbdg/graph.cpp` plus `src/lancet/core/variant_builder.cpp`.

Expected probe call surface (from `docs_dev/subsystems/probe_tracking.md` Phase 3):

- `ProbeGenerateAndTag`, `ProbeCountInReads`
- `ProbeLogStatus(stage, ctx)` for stages BUILD, LOWCOV1, COMPRESS1, LOWCOV2, COMPRESS2, TIPS
- `ProbeSetNoAnchor`, `ProbeSetShortAnchor`, `ProbeCheckAnchorOverlap`
- `ProbeSetGraphCycle`, `ProbeSetGraphComplex`, `ProbeSetTraversalLimit`
- `ProbeCheckPaths`
- `CheckMsaExtraction`, `CheckGenotyperResult` (in `probe_diagnostics`)

Verify every method named in `docs_dev/subsystems/probe_tracking.md` exists in the header. Report any documented-but-missing methods, plus any methods present in the header but not documented (the latter is a documentation gap).

### Pairing 4 — TSV schema column count

**Source of truth:** the `WriteHeader` and `WriteRow` functions in `src/lancet/cbdg/probe_results_writer.cpp`.

`docs_dev/subsystems/probe_tracking.md` claims a 40-column TSV. Count the actual emitted columns. Report any drift.

The Python attribution side reads these columns by name; if column names change, both the C++ and the Python side must change in coordinated fashion. A name change on only one side is a silent failure mode.

**Bundle locations to check:** any prose claim about TSV column count or column names in `probe-tracking/SKILL.md`, `probe-interpreter.md`, or `docs_dev/subsystems/probe_tracking.md`.

### Pairing 5 — Script CLI shape

**Source of truth:** `scripts/truth_concordance.py` and `scripts/analyze_probe_results.py` argparse blocks.

For each script, extract every flag, its required-ness, and its default. Compare to the flag tables in:

- `.claude/commands/probe-concordance.md` — `truth_concordance.py` flags
- `.claude/commands/probe-analyze.md` — `analyze_probe_results.py` flags
- `.claude/commands/probe-run.md` — Lancet2 `--probe-variants` / `--probe-results` / `--verbose` flags (from `cli_interface.cpp`)
- `.claude/skills/probe-tracking/SKILL.md` — flag tables
- `docs_dev/subsystems/probe_tracking.md` — Step 1 / Step 2 / Step 3 flag tables

Report any flag claimed by the bundle that does not exist, any required flag that the bundle treats as optional (or vice versa), and any flag in the script that is not surfaced anywhere in the bundle (the latter is a documentation gap, not a correctness issue).

## Report format

```
## /audit-probe-pipeline drift report

### Summary
<N total findings across <K of 5> pairings>

### Pairing 1 — Stage cascade
Status: <CLEAN | DRIFT FOUND>
Findings:
  - bundle/<file>:<line> mentions stage `<name>`; not in STAGE_ORDER
  - STAGE_ORDER includes `<name>`; not mentioned in probe-interpreter.md
  - ...

[... pairings 2–5 ...]

### Suggested fixes
<for each finding, a proposed edit; one bullet per file>

### Manual review needed
<findings the script can't auto-resolve, e.g., a category restructure>
```

After producing the report, ask the user how to proceed: apply all auto-fixable findings, apply selectively, or report-only mode. Do not apply without explicit approval.

## When NOT to use this command

Do not use this command if the probe-tracking source files are mid-refactor (uncommitted partial changes). The audit needs the source to be parseable; a half-finished refactor produces noisy false-positive findings.

Do not use this command for runtime correctness checks. The audit catches documentation drift; it does not validate that the probe pipeline still produces correct attributions. For correctness, run `/probe-concordance` → `/probe-run` → `/probe-analyze` end-to-end on the known germline fixture and verify the report sanity.

Do not use this command for general bundle drift. The probe pipeline has its own command because its surface is wide; for everything else, `/audit-bundle` is the right tool.

## Maintenance

The five pairings here are the surfaces where bundle claims about probe tracking can drift. If a new surface emerges (e.g., a new helper script), add a sixth pairing.

The expected sub-stage list in pairing 2 will grow when new StatusCode values are added in `variant_builder.cpp`. Keep that table in sync with the actual enum.
