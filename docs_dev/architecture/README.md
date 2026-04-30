# Architecture

This directory holds **cross-cutting design documentation** for Lancet2. Two formats live here:

- **Architecture Decision Records (ADRs)** — short, numbered documents that capture a single architectural decision, the alternatives considered, and the consequences. Use when you've made a decision that affects how the codebase is structured and a future contributor will need to know *why* the system is the way it is.
- **Architecture Overview Documents** — longer prose documents that explain the shape of a system across multiple subsystems. Use when describing a cross-cutting structural pattern (e.g., "the six-layer dependency rule") or a system-wide concern (e.g., "how Lancet2 manages parallelism").

Both formats share a quality bar that is stricter than the rest of `docs_dev/`. This document defines that bar.

## Lifecycle: evergreen

Architecture documents are **evergreen**. They describe the system as it is *now*, and they must be kept current as the system evolves. The synchronization rule (`../style/sync_and_verification.md`) applies in full.

This is the load-bearing distinction from `../investigations/`, which captures snapshots of past moments and is immutable. If you are documenting a past investigation, that goes in `investigations/` — the document records *what was true at that moment*. If you are documenting how the system *currently is and should remain*, that goes here — the document records the *intended* shape of the system.

When an architecture document goes out of date, fix it. Do not write a new one to supersede it.

## When to write an ADR vs an overview

| Situation | Format |
|:---|:---|
| You're choosing one option from several, and the choice has architectural consequences | ADR |
| Future contributors will ask "why did we do it this way?" and the answer needs context | ADR |
| The decision constrains future work in non-obvious ways | ADR |
| You're explaining the shape of a system that already exists, end-to-end | Overview |
| The explanation will exceed ~600 lines if done well | Overview |
| The document needs to walk through multiple subsystems' interactions | Overview |

A useful disambiguation: if the document is fundamentally answering *"why this and not that?"*, it's an ADR. If it's answering *"how does this system fit together?"*, it's an overview.

If you find yourself wanting to write an ADR longer than ~250 lines or an overview shorter than ~150, you've probably picked the wrong format. Reconsider before continuing.

# Part 1 — Architecture Decision Records (ADRs)

## What an ADR is

An ADR captures **one** architectural decision. It explains the situation that prompted the decision, the options considered, the choice made, and the consequences (both the constraints accepted and the doors closed).

The ADR format is well-known in software engineering. The variant Lancet2 uses follows Michael Nygard's original format with light modifications. The point of the format is consistency — a reader who has seen one project's ADRs can navigate any other project's ADRs immediately.

## File naming

ADRs are numbered sequentially and named by short kebab-case slug:

```
0001-pixi-over-conda.md
0002-six-layer-dependency-rule.md
0003-spoa-int16-simd-path.md
```

The number is zero-padded to four digits. Numbers are assigned in the order ADRs land in main, not when they're drafted. Find the highest existing number and increment by one. If two ADRs are drafted in parallel, the second to merge gets the next number; rebase if needed.

The slug captures the decision in 2–5 words. Match the verb tense of the section "Decision" — i.e. write as a present-tense statement of the chosen position, not as a question or as the prior state.

## Required sections

Every ADR has exactly these sections, in this order. Do not add or remove top-level sections.

```markdown
# ADR 0001: <slug as title — capitalized>

**Status:** <Proposed | Accepted | Superseded by ADR ####>
**Date:** YYYY-MM-DD
**Last reviewed:** YYYY-MM-DD

## Context

What forced the decision? What was true about the system, the constraints,
the team, the requirements? Two or three paragraphs. Concrete and verifiable
against the codebase or external constraints — not aspirational.

## Decision

What did we decide? State the chosen position as a single sentence at the top,
then expand. Do not bury the decision under the alternatives.

## Alternatives considered

What else was on the table? For each alternative: what it was, why it was
considered, why it lost. Be specific — "performance" is not a reason; "adds
~15% wall-clock cost on the standard somatic fixture" is.

## Consequences

What does this decision constrain? What follow-on work or non-obvious
trade-offs does it create? This section is where future contributors learn
why working around the decision is harder than working with it.
```

The four sections are non-negotiable because together they answer the four questions a future reader will ask: *what was the situation?* *what did we choose?* *what did we reject?* *what did we accept in exchange?*

## Voice and tone

- **Written in past tense for Context, present tense for Decision and Consequences.** The context describes what was; the decision is the position we now hold.
- **Specific, not abstract.** Cite file paths, function names, version numbers, measured performance, dependency requirements. *"Conda's environment-creation time exceeded 90s on cold runs"* is useful; *"Conda was slow"* is not.
- **Honest about trade-offs.** An ADR that lists no consequences is a red flag — every architectural decision closes some doors. Name them.
- **No marketing language.** *"This robust solution elegantly addresses the underlying complexity"* is noise. State the choice and the consequences.
- **No hedging.** ADRs are the historical record of decisions made. Once accepted, they speak with authority.

## Status field

- **Proposed** — drafted, not yet accepted. Open in a PR for review.
- **Accepted** — merged. The decision is in force.
- **Superseded by ADR ####** — replaced by a later decision. The original ADR stays in place; do not delete it. Add a header line at the top:

  ```markdown
  > **Note:** This decision was superseded by [ADR 0017](0017-new-decision.md).
  > Kept as historical record.
  ```

ADRs do not get a "Rejected" status. If a proposal is rejected, either delete the draft (it never landed) or capture the rejection in a new ADR that records the decision *not* to do the thing.

## Length and detail

| Section | Typical length |
|:---|:---|
| Context | 2–4 paragraphs (~150–300 words) |
| Decision | 1 sentence + 1 paragraph elaboration (~50–150 words) |
| Alternatives considered | 2–4 alternatives, each ~50–150 words |
| Consequences | 1 paragraph (~100–200 words) |

A complete ADR is typically 100–250 lines. Above 250 lines and you are probably writing an overview document instead, or trying to bundle multiple decisions into one ADR (split them).

## Cross-references

Every ADR that supersedes or refines an earlier one must link to it. Every ADR that constrains a subsystem documented in `../subsystems/` should link both ways.

## ADR template (copy-paste)

```markdown
# ADR ####: <Decision in 2–5 words>

**Status:** Proposed
**Date:** YYYY-MM-DD
**Last reviewed:** YYYY-MM-DD

## Context

<What was true that forced this decision? Two or three paragraphs.>

## Decision

<One-sentence statement of the chosen position.>

<Expansion paragraph: what does adopting this look like in concrete terms?>

## Alternatives considered

### <Alternative 1 name>

<What it was. Why it was considered. Why it lost — be specific.>

### <Alternative 2 name>

<Same shape.>

## Consequences

<What this decision constrains. What follow-on work it implies. What doors
it closes. What non-obvious trade-offs future contributors should know about.>
```

# Part 2 — Architecture Overview Documents

## What an overview document is

An overview describes the **shape of a system** that already exists in the codebase. It explains how multiple subsystems relate, what cross-cutting patterns govern behaviour across them, and what mental model a developer needs to navigate the source.

Overviews are not feature documentation (that's `subsystems/`) and not user-facing (that's `docs/guides/`). They sit at a higher altitude: a developer about to make a substantive change to multiple parts of the codebase reads an overview to orient before diving in.

## Examples of what would warrant an overview document

- **`pipeline_architecture.md`** — the end-to-end execution shape: how the pipeline runner constructs the executor, how windows are batched, how parallelism is structured, how output is serialized.
- **`layer_dependencies.md`** — the six-layer dependency rule, what layers may include from what, what enforcement exists, why the constraint matters.
- **`memory_model.md`** — how memory is allocated and freed across the pipeline, what allocator is used, where the hot allocation paths are.
- **`error_handling.md`** — how errors propagate, what's a fatal vs recoverable failure, what the user sees vs what's logged.
- **`testing_strategy.md`** — what's covered by unit tests, what's covered by E2E, what's covered by truth-set comparison via the probe pipeline.

These are illustrative — none of them need to exist before they're worth writing, and several may not be worth writing at all unless a real reader will benefit.

## File naming

Lowercase, underscore-separated, descriptive. No number prefix (overviews are not sequenced like ADRs are):

- `pipeline_architecture.md`
- `parallelism_model.md`
- `layer_dependencies.md`

## Structure

Overviews don't have a fixed section list the way ADRs do — the right structure depends on what's being explained. But every overview document should answer four questions, in roughly this order:

1. **What is this?** A 2–3 paragraph orientation: what system is being described, why it matters, who needs to read this document.
2. **What's the mental model?** The simplest correct picture of how the system fits together. ASCII diagrams, tables, or short prose. The reader should be able to navigate the source after this section even if they read nothing else.
3. **How does it actually work?** The detailed walkthrough. Reference specific files, classes, functions. Cite line ranges or function names. Show the code path for a typical case.
4. **What's non-obvious?** The traps, the historical reasons, the constraints that surprise newcomers. This is the section that will save future readers the most time.

## Voice and tone

- **Written for a developer about to make a substantive change.** They have the source open, they know C++, they don't know this system. Write to that person.
- **Specific over abstract.** Always cite file paths and function names. *"BuildHaplotypes (graph.cpp:348)"* is more useful than *"the haplotype builder"*.
- **Diagrams when they help.** ASCII diagrams for control flow, tables for stage-by-stage walkthroughs, mathematical notation for formulas. Use the same conventions as the code-comments style guide.
- **No marketing language, no hedging, no filler.** Same Core Principles as the rest of `docs_dev/`.
- **Cross-reference generously.** Overview documents are the natural place to surface "see ADR 0007 for why we chose this" or "subsystem X has a deep-dive at `../subsystems/x.md`."

## Length and detail

Overview documents typically run 200–800 lines. Below 200 they're usually trying to be ADRs (and would fit there better). Above 800 they're usually trying to be subsystem deep-dives (and should be moved to `../subsystems/`).

The right length is "as long as it needs to be, and no longer." A 600-line overview that earns every line is fine. A 400-line overview where 100 lines are filler is too long.

## Required header

Every overview document has these three lines as its first content under the title:

```markdown
**Type:** Overview
**Last reviewed:** YYYY-MM-DD
**Audience:** <one sentence describing who this is for>
```

The "Last reviewed" date is what distinguishes evergreen documents from snapshots — a reader can tell at a glance whether the document has been kept current. Update it when you substantively edit the document or confirm it still reflects reality.

## Overview template (copy-paste)

```markdown
# <Title>

**Type:** Overview
**Last reviewed:** YYYY-MM-DD
**Audience:** <one sentence>

<2–3 paragraph orientation: what is this, why it matters, what the reader will learn.>

## Mental model

<The simplest correct picture. Diagram, table, or short prose. The reader
should be able to navigate the source after this section alone.>

## How it actually works

<The walkthrough. Cite files, classes, functions. Show the code path for a
typical case.>

### <Sub-aspect 1>

<...>

### <Sub-aspect 2>

<...>

## Non-obvious behaviour

<The traps, the historical reasons, the constraints that surprise newcomers.>

## Cross-references

- ADRs that constrain this system: ADR ####, ADR ####
- Subsystem deep-dives that live below this overview: `../subsystems/<>.md`
- Public guides that explain user-visible behaviour: `../../docs/guides/<>.md`
```

# Existing documents

<!-- Update this list as architecture documents are added. -->

*No architecture documents have been written yet. This section will populate as ADRs and overviews are added.*

When the first document lands, add an entry here in this format:

- **ADR 0001:** [Pixi over Conda](0001-pixi-over-conda.md) — chose Pixi for environment management.
- **Overview:** [Pipeline architecture](pipeline_architecture.md) — end-to-end execution shape.

# Cross-references

- **Style hub:** `../style/README.md` — Core Principles and the four focused style documents apply here too.
- **Sync and verification:** `../style/sync_and_verification.md` — the synchronization rule applies in full to architecture documents (they are evergreen, like the style guide).
- **Investigations:** `../investigations/README.md` — for documenting what was learned during a past incident or debug session, which is a different lifecycle (immutable snapshots).
