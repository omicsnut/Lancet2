# /audit-vcf-schema — Focused VCF FORMAT/INFO/FILTER drift check

Walk every `##FORMAT`, `##INFO`, and `##FILTER` line in `src/lancet/cli/vcf_header_builder.cpp` and verify that the bundle's claims about VCF schema match. This is the focused complement to `/audit-bundle` for the case where you want to check VCF schema specifically without running the full audit.

## When to run this

After every release that touches VCF emission. After any change to `vcf_header_builder.cpp`, `variant_call.{h,cpp}`, `variant_support.cpp`, `sample_format_data.{h,cpp}`, or `variant_extractor.{h,cpp}`. Before invoking `vcf-validator` on a PR that touches schema, as a sanity check that the bundle's invariants list is current.

The VCF schema is the surface that drifts most often (every release adds, modifies, or retires fields), so this command is intentionally fast and narrow. Run it freely.

## Procedure

### Phase 1 — Enumerate the current schema

Read `src/lancet/cli/vcf_header_builder.cpp` and extract every header line. Each line has the form:

```
##FORMAT=<ID=NAME,Number=N,Type=T,Description="...">
##INFO=<ID=NAME,Number=N,Type=T,Description="...">
##FILTER=<ID=NAME,Description="...">
```

Build a structured list:

```
FORMAT fields: [(ID, Number, Type, Description), ...]
INFO fields:   [(ID, Number, Type, Description), ...]
FILTER fields: [(ID, Description), ...]
```

If the source has been moved (e.g., header construction migrated to a new file), report that as the first finding and stop — the rest of the audit depends on having the right source file.

### Phase 2 — Compare against bundle claims

The bundle makes specific claims about VCF schema in the following files. Read each, extract the field-related claims, and compare against the structured list from phase 1.

**`.claude/agents/vcf-validator.md`** — the "Schema invariants you must enforce" section codifies version-agnostic rules (always-present fields, missing-value conventions, Number cardinality conventions, score and metric semantics, field-range conventions, per-ALT vs per-record semantics, multi-sample generalization). For each rule and any specific FORMAT IDs it names (e.g., AD, DP, GT, GQ, PL, ASMD, AHDD, BQCD, CMLOD, FLD, FSSE, HSE, MQCD, PANG, PDCV, PRAD, RPCD, SB, SDFC, SCA, NPBQ, RMQ), check:
- Is the field still present in the header? If not, the invariant is stale.
- Does the field's `Number` and `Type` match what the agent says? If not, flag.
- Does the missing-value behavior described still apply? (Hard to verify from header alone; flag for manual review if the field's emission code might have changed.)

**`.claude/agents/assembly-and-calling-expert.md`** — the "VCF FORMAT/INFO/FILTER semantics" bullet should now read as a delegation pointer to `vcf-validator` rather than carry its own invariants list. If the bullet has re-accumulated specific schema claims (e.g., naming specific fields and their Number/Type/missing-value semantics), that is a finding — those claims belong only in `vcf-validator.md` to prevent the dual-source drift the consolidation was meant to eliminate. Verify the bullet still ends with the delegation language.

**`AGENTS.md`** — the downstream-sync section names `variant_call.{h,cpp}` and `variant_support.cpp` as the field-definition and value-computation locations. Verify these files still exist and contain field-related code (a sanity check, not an exhaustive grep).

### Phase 3 — Identify new fields the bundle does not yet acknowledge

For every FORMAT/INFO/FILTER field present in the header but NOT mentioned anywhere in the bundle, list them as "new fields the bundle does not document." For some fields this is fine (the bundle does not need to enumerate every field), but if a field is non-obvious, it may earn a mention in `vcf-validator`'s invariants list.

The threshold for "earns a mention" is rough: fields with unusual cardinality (`Number=A`, `Number=R`, `Number=G`), fields with breaking-change-prone semantics (missing-value rules, per-sample vs per-site), and fields whose existence implies architectural assumptions (e.g., a per-ALT field implies the multi-ALT path is supported). Trivial Number=1 Type=Integer fields can be ignored.

### Phase 4 — Identify removed fields the bundle still references

For every field name mentioned in the bundle but NOT present in the header, flag as "removed field still referenced." This is the most actionable category — every such reference is wrong and should be removed.

### Phase 5 — Produce the report

Format:

```
## /audit-vcf-schema drift report

### Schema summary
- N FORMAT fields in header
- N INFO fields in header
- N FILTER fields in header

### Findings

#### Stale claims (field changed semantically)
- bundle/<file>:<line> says X about field FOO; header now says Y
- ...

#### Removed fields still referenced
- bundle/<file>:<line> mentions FOO; not in header
- ...

#### New fields not yet acknowledged (judgment call)
- BAR (Number=A, Type=Integer) added; bundle does not mention; suggest adding to vcf-validator if it's a non-trivial field.
- ...

### Suggested fixes
<for each stale or removed finding, the specific edit to apply>
```

After producing the report, ask the user which fixes to apply. Do not apply without explicit approval — schema claims are subtle and a wrong "fix" can make the bundle worse.

## When NOT to use this command

Do not use this command if `vcf_header_builder.cpp` has unsaved or uncommitted edits — the audit depends on a stable source. If you are in the middle of a schema change, finish the change first, commit, then run this.

Do not use this command to plan schema changes (i.e., "what should the new field's Number be?"). That is the `vcf-validator` subagent's job; it has the v4.5 spec.

Do not use this command as a substitute for `vcf-validator` review on a PR. The audit catches drift between the bundle and the source; it does not validate that a new field is well-formed against the v4.5 spec. Both are useful and they do different jobs.

## Maintenance

This command's pairing list (which bundle files to check, what kind of claims they make about VCF) should grow as the bundle accumulates more VCF-related claims. If a new agent or skill makes a specific field claim, add the file to phase 2.

The "earns a mention" judgment in phase 3 is the part most likely to be wrong. If the report consistently flags fields you consider trivial, tighten the threshold. If it consistently misses fields you care about, loosen it.
