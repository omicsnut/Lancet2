# /probe-analyze — Step 3 of probe tracking: attribution + interpretation handoff

Run `scripts/analyze_probe_results.py` to attribute each missed variant to a pipeline stage via the 27-level cascade, then hand off to the `probe-interpreter` subagent for focused recommendations grounded in the C++ source.

This command bridges the operational layer (running the script, producing the report) with the interpretive layer (what does the report mean, where do I dig in next). The script is the operational side; the subagent is the interpretive side. Together they replace the failure mode where Claude reads a report cold, focuses on the wrong section, and recommends generic next steps.

## When to use

Run this when you have a `probe_results.tsv` from step 2 (`/probe-run`) and you want stage attribution + a focused recommendation. Run it again with `--view <name>` for iterative single-section debugging without re-rendering the whole report.

For "I have a report but want to re-interpret it without re-running the script" — invoke the `probe-interpreter` subagent directly instead.

## Procedure

### Phase 1 — Verify inputs

If the user did not supply paths, default to the most recent `notes/probe-debug-<YYYY-MM-DD>/` directory. Check:

```bash
ls <outdir>/probe_results.tsv <outdir>/missed_variants.txt
```

`probe_results.tsv` and `missed_variants.txt` are required. If either is missing, surface the gap and suggest the prior step (`/probe-run` or `/probe-concordance` respectively).

Optional but high-value:

- `concordance_details.txt` — enriches §1 scorecard. If present in the same directory, pass it.
- `lancet2_debug.log` — enables sub-classification of `not_processed` into 6 sub-stages. If present in the same directory, pass it via `--log`.

If the user wants to skip a single section or render only one, gather `--view` via `AskUserQuestion` from the choices: `scorecard`, `funnel`, `survival`, `breakdown`, `genotyper`, `targets`, `deepdive`, `windows`, `all`. Default to `all`.

### Phase 2 — Run the script

```bash
pixi run -e hts-tools python3 scripts/analyze_probe_results.py \
    --probe-results <outdir>/probe_results.tsv \
    --missed-variants <outdir>/missed_variants.txt \
    [--concordance-details <outdir>/concordance_details.txt] \
    [--log <outdir>/lancet2_debug.log] \
    --output-dir <outdir> \
    --view <view>
```

The script writes:

- `probe_analysis_report.txt` — full rich-text report (8 sections)
- `probe_stage_attribution.txt` — TSV: one row per probe with final `lost_at_stage`, type, size
- `probe_survival_matrix.txt` — TSV: one row per `(probe, window, comp, k)` with 6 survival counts

It also prints the report to stdout. Capture both stdout (the rich report) and stderr (any progress / warnings).

### Phase 3 — Surface the top-line attribution

After the script completes, parse the §2 funnel from `probe_analysis_report.txt` and surface the top 3 stages with percentages in chat. Also surface the §1 scorecard's coverage gap; if it is nonzero, flag it as a concern (incomplete data).

The point of this surface is to give the user enough context to decide whether to invoke the interpreter or read the report themselves. Do not produce a summary that pre-empts the interpreter; just surface the dominant stage(s).

### Phase 4 — Hand off to probe-interpreter

End the response with a structured handoff:

> Step 3 complete. Report: `<outdir>/probe_analysis_report.txt`. Top attribution from §2 funnel: \[stage1 X%, stage2 Y%, stage3 Z%\]. Coverage gap: N (should be 0).
>
> To get a focused recommendation grounded in the C++ source — what code is responsible for the dominant stage, what to investigate first, sensitivity vs specificity trade-offs — invoke the `probe-interpreter` subagent.

Do not auto-invoke the subagent. The user may want to read the report themselves first, or pick a particular `--view` to drill in. The subagent invocation should be a deliberate "I want recommendations now" decision.

The handoff structure matters. Without it, Claude tends to either (a) produce a generic recommendation in the slash command itself (which the operational/interpretive split is designed to prevent), or (b) say "the report is at <path>, please read it" and stop (which isn't useful). The structured handoff threads the needle: enough surface to decide what to do next, no pre-empting the interpreter's deeper read.

## Iterative usage

For follow-up work without re-running the script:

```
/probe-analyze --view funnel    # render just §2
/probe-analyze --view targets   # render just §6
```

The script is fast; re-rendering specific views during iteration is cheap. The interpreter is the expensive part; invoke it once after you have the views you want to discuss, not after every script run.

For re-interpreting an existing report (after a Lancet2 source change, for example):

```
/probe-run                       # capture fresh probe_results.tsv
/probe-analyze --view all        # fresh attribution
[invoke probe-interpreter manually with reference to the new outputs]
```

## When NOT to use this command

Do not use this command without `probe_results.tsv` from step 2. The script does not produce useful output without it.

Do not use this command for one-off "what does stage X mean" questions when you don't have a report to interpret. The `probe-tracking` skill documents stages in the abstract; use that.

Do not use this command as a substitute for invoking the `probe-interpreter` subagent. The slash command does the script run + handoff; the subagent does the interpretive work. They compose; they don't replace each other.

## Maintenance

The 8-section view list and the 27-stage cascade reflect `scripts/analyze_probe_results.py`. If the script adds, removes, or renames sections or stages, update this command. The `/audit-probe-pipeline` slash command catches drift between the bundle's claims and the actual script.

The handoff prompt format is the leverage point. If the interpreter consistently fails to fire because the prompt doesn't surface enough context, tighten the surface; if it fires too eagerly without the user wanting it, soften the wording. Treat the prompt as a description-as-trigger style problem.
