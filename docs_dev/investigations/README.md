# Investigations

This directory holds **immutable snapshots of past debugging, postmortems, and performance investigations**. Three closely related document types live here:

- **Postmortems** — write-ups of incidents, regressions, or production-impacting bugs. *"What went wrong, why we didn't catch it earlier, what changed in response."*
- **Debug archaeology** — deep technical investigations of subtle bugs that took serious time to track down. *"Here is the bug, here is the path that led to finding it, here is the trap that future readers should know about."*
- **Performance investigations** — write-ups of "why was X slow, what we tried, what worked, what we learned." Often paired with the resulting `perf:` commit.

The three types share enough structure that one quality bar covers them all. This document defines that bar.

## Lifecycle: immutable

Investigations are **immutable once merged**. They are a snapshot of what was understood at a specific moment in time, captured before the understanding fades and while the artefacts (logs, profiles, repros) are still fresh.

Do not edit an investigation to reflect later understanding, or to "fix" claims that turned out to be wrong. If understanding changes, write a new investigation that supersedes the old one, with a header that links back to the original. The original stays as historical record.

This is the load-bearing distinction from `../architecture/`, which is evergreen and is kept current as the system evolves. If the document records *what was true at that moment*, it goes here. If it records *the intended shape of the system going forward*, that's an architecture document.

The synchronization rule (`../style/sync_and_verification.md`) **does not apply** to investigation documents. They are deliberately not maintained against the evolving codebase. A reader of an investigation must understand they are reading history, not current behaviour.

## When to write an investigation document

Not every bug fix or performance tweak needs an investigation document. Most don't. Write an investigation when at least two of the following are true:

- The investigation took meaningful time (more than a working day) and the path to the answer was non-obvious.
- The bug or perf issue was caused by a subtle interaction between subsystems, a violated invariant, or a misunderstood library behaviour — something a future contributor could plausibly hit again.
- A future maintainer working in the same area would benefit from knowing what was already tried and ruled out.
- The fix or finding is in production but the *reasoning* behind it is not obvious from the diff or the commit message.

If the bug was found in five minutes and the fix is one line, do not write an investigation. The commit message carries enough context.

If the issue is "user reported flaky behaviour, can't reproduce" and the resolution is "couldn't reproduce, closing the issue," do not write an investigation. There's nothing to capture.

The bar is: *will a future contributor stuck on something similar be glad this exists?*

## When to write each type

| Type | Trigger | Primary value |
|:---|:---|:---|
| Postmortem | Incident or regression with user-visible impact | Why we didn't catch it earlier; what process or test changed |
| Debug archaeology | A subtle bug that took serious effort to find | The path through the symptoms, the false leads, the eventual root cause |
| Performance investigation | "Why is X slow" that required real measurement work | What was measured, what was tried, what worked, why |

The distinctions can blur. Don't over-think which type fits — the structural template below works for all three. Pick the type that best describes the situation; you can mention "this began as a perf issue and became a debug archaeology" in the document itself.

## File naming

Investigations are named by date and short slug:

```
2026-04-15-msa-shifted-coordinate-bug.md
2026-Q1-genotyper-regression-postmortem.md
2025-09-23-extractor-cache-locality-perf.md
```

Two date conventions are accepted:

- **Specific date** (`YYYY-MM-DD`) when the investigation happened on or around a specific day.
- **Quarter** (`YYYY-Q#`) when the investigation spanned a longer period or the exact start date doesn't matter.

The slug is descriptive lowercase kebab-case. Keep it specific — `2026-04-15-bug.md` is useless; `2026-04-15-msa-shifted-coordinate-bug.md` lets a future grep find it.

Date prefixes serve double duty: they sort chronologically in `ls` and they signal the immutable-snapshot nature of the document. A document without a date in its filename is suspicious in this directory.

## Required structure

Every investigation document has these sections, in this order:

```markdown
# <Type prefix>: <Title>

**Type:** Postmortem | Debug Archaeology | Performance Investigation
**Captured on:** YYYY-MM-DD
**Status:** Final (or "Superseded by [later](path.md)")
**Author:** <name or git handle>

## Summary

One paragraph (~5 sentences). What was the problem? What was found? What was
done about it? A reader should be able to read this paragraph alone and
decide whether the rest of the document is relevant to them.

## What we observed

The symptoms — concrete and verifiable. Logs, error messages, performance
numbers, screenshots, test failures. State what was seen, not what we
suspected at the time.

## Investigation path

How we got from symptom to root cause. Be honest about the false leads —
"we initially suspected X because of Y, but ruled it out when Z" is more
useful than a clean retroactive narrative. Future readers stuck on similar
symptoms benefit from knowing what was tried and didn't work.

## Root cause

What was actually wrong. Cite specific code, specific commits if relevant,
specific assumptions that turned out to be false. Be concrete.

## Resolution

What changed? Link to the fix commit or the PR. If the fix was non-trivial,
explain the reasoning briefly — but do not duplicate the commit message.

## What we learned

Lessons that would benefit a future contributor. What invariant was violated?
What test gap allowed this through? What documentation was missing or wrong?
What was the surprising behaviour of a library or tool? Be specific enough
that someone hitting similar symptoms can recognize the pattern.

## What we did not change

Honest section: what *should* arguably have changed but didn't, and why. If
a class of bugs is now known to exist but we chose not to fix them all, say
so. If a refactor would prevent recurrence but is too expensive, say so.
This is the section that prevents readers from later asking "why didn't they
fix this properly?"
```

The seven sections are non-negotiable because together they answer what a future contributor will ask: *what happened, what did you see, how did you find it, what was actually wrong, what changed, what should I learn, and what's still not fixed?*

## Voice and tone

- **Written in past tense.** Investigations describe events that already happened. *"We observed that..."*, *"The first hypothesis was..."*, *"After eliminating..."*.
- **First-person plural ("we") is acceptable.** This is the one place in `docs_dev/` where "we" reads naturally — investigations are reports of what a person or team did, and the genre conventions allow it. Singular "I" is also acceptable for solo investigations. Avoid the awkward third-person "the investigator" formulation.
- **Honest, not retroactive.** Do not present the investigation as if the answer was obvious from the start. Document the false leads. *"We spent three days suspecting the genotyper before noticing the issue was in the MSA coordinate translation"* is the kind of detail that saves future readers time.
- **Specific to the moment.** Cite commit SHAs (the *current at-time-of-writing* SHA), file paths, line numbers, library versions, build flags. The codebase will move on; the investigation is anchored to its moment.
- **Concrete, no marketing language.** Investigations are reports. Write them as reports.
- **No filler, no hedging.** Same Core Principles apply. The audience is a developer in the same shoes you were in.

## Length and detail

Investigations typically run 150–600 lines. The "right" length depends on the complexity of the path from symptom to root cause. A bug found in five minutes with a one-line fix doesn't need an investigation document at all (see "When to write" above). A subtle bug that took a week to track down may earn 500 lines.

Be specific over comprehensive. A 200-line investigation that pinpoints the root cause precisely is worth ten times more than an 800-line investigation that hedges and gestures vaguely.

## Voice differences from `architecture/`

The voice contrast with architecture documents is intentional and worth surfacing:

| Aspect | Architecture (evergreen) | Investigations (immutable) |
|:---|:---|:---|
| Tense | Present ("the pipeline runs windows in parallel") | Past ("we observed that windows were not parallelizing") |
| Subject | "Lancet2", "the system" | "We", "I" |
| Authority | Authoritative — describes how things are | Reportorial — describes what happened |
| Maintenance | Updated as system evolves | Frozen at capture date |
| Audience | Developer about to change the system | Developer hitting similar symptoms later |

Both types share the underlying Core Principles (source-as-truth, density, why-not-just-what, quantify, address-the-real-question). The difference is what kind of question they answer.

## Investigation template (copy-paste)

```markdown
# <Type>: <Specific descriptive title>

**Type:** Postmortem | Debug Archaeology | Performance Investigation
**Captured on:** YYYY-MM-DD
**Status:** Final
**Author:** <name or git handle>

## Summary

<One paragraph, ~5 sentences. Problem, finding, action.>

## What we observed

<Concrete symptoms. Logs, errors, perf numbers, test failures.>

## Investigation path

<How we got from symptoms to root cause, including the false leads.>

## Root cause

<What was actually wrong. Cite code and assumptions concretely.>

## Resolution

<What changed. Link to the commit/PR. Explain non-obvious reasoning briefly.>

## What we learned

<Lessons for future contributors. Invariants violated, gaps in tests/docs.>

## What we did not change

<Honest section on what should arguably have changed but didn't, and why.>
```

## Superseded investigations

If a later investigation supersedes an earlier one, the earlier document gets a header line at the top:

```markdown
> **Note:** This investigation has been superseded by [2026-08-12-...md](path.md).
> Kept as historical record. The conclusions in §"Root cause" below were later
> revised in light of new information; see the superseding document for current
> understanding.
```

The earlier investigation is not edited beyond this header. The reasoning, the false leads, and the conclusions all remain visible — they are still useful as a record of what was understood at the time.

# Existing investigations

<!-- Update this list as investigation documents are added. Most-recent at top. -->

*No investigation documents have been written yet. This section will populate as investigations are captured.*

When the first investigation lands, add an entry here in this format:

- **2026-04-15** — [MSA shifted coordinate bug](2026-04-15-msa-shifted-coordinate-bug.md) (Debug Archaeology) — coordinate translation off-by-one in `variant_extractor.cpp` for windows containing homopolymer runs.
- **2026-Q1** — [Genotyper regression postmortem](2026-Q1-genotyper-regression-postmortem.md) (Postmortem) — silent regression in PBQ score weighting introduced in v2.8.0.

# Cross-references

- **Architecture:** `../architecture/README.md` — for documenting decisions and overviews about how the system *is* and *should remain*. Different lifecycle (evergreen vs immutable).
- **Style hub:** `../style/README.md` — Core Principles apply here too. Voice differences from website docs and architecture docs are documented above.
- **Sync rule does not apply.** The synchronization rule in `../style/sync_and_verification.md` deliberately excludes this directory. Investigations are not kept current with the evolving codebase by design.
