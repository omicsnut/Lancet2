# Subsystems

This directory holds **deep-dives on how a specific Lancet2 subsystem works at the source level**. Audience: a developer about to debug or extend that subsystem. Voice: explanatory and developer-facing — what the subsystem does, why it's structured this way, how to read its source without getting lost.

## What belongs here

A document earns its place in `subsystems/` when:

- It explains the internal mechanics of a single subsystem (the probe-tracking pipeline, the de Bruijn graph builder, the genotyper, the read-evidence engine).
- The audience is a developer about to read or modify the subsystem's source.
- The level of detail is too implementation-specific for `docs/guides/` (which is user-facing) but more focused than an `architecture/` overview (which spans multiple subsystems).
- The document will be kept current as the subsystem evolves — it's evergreen, like architecture documents and the style guide.

If the document is about *how the user interacts with a subsystem*, that's a public guide and goes in `docs/guides/` instead. If the document is about *how multiple subsystems fit together*, that's an architecture overview and goes in `../architecture/`.

## Voice and structure

Subsystem deep-dives don't have a fixed structure the way ADRs and investigations do — the right shape depends on what's being explained. But every subsystem document should answer:

1. **What is this subsystem?** A 2–3 paragraph orientation: where in the source it lives, what it produces or operates on, who calls it and who it calls.
2. **What's the conceptual model?** The mental picture a reader needs before they touch the source.
3. **What does the source actually do?** A walkthrough with file:line references. Cite specific functions, specific data structures, specific flow paths.
4. **What's non-obvious?** Invariants that aren't visible in the code, historical reasons for surprising choices, traps that caught earlier developers.

Apply the Core Principles from `../style/README.md`. Use the comment-style conventions from `../style/code_comments.md` for any included pseudocode or formula notation. Cite specific source paths — *"BuildHaplotypes (graph.cpp:348)"* not *"the haplotype builder"*.

The existing `probe_tracking.md` is the prototype. New subsystem documents should match its information density and citation discipline.

## Lifecycle

Subsystem documents are **evergreen**. The synchronization rule (`../style/sync_and_verification.md`) applies in full. When a subsystem changes substantively, the document changes in the same commit. Stale subsystem documentation is worse than missing subsystem documentation — readers trust it and follow it into a wall.

## File naming

Lowercase, underscore-separated, named for the subsystem:

- `probe_tracking.md`
- `genotyper.md`
- `cbdg_graph_builder.md`

# Existing subsystem documents

- **[Probe tracking pipeline](probe_tracking.md)** — the variant forensic pipeline (truth_concordance → Lancet2 with --probe-variants → analyze_probe_results), the 27-stage attribution cascade, the C++/Python contract.

# Cross-references

- **Architecture:** `../architecture/` — overview documents that span multiple subsystems live there.
- **Workflows:** `../workflows/` — runbooks for *doing* something with a subsystem (e.g., "how to capture and analyze a profile") live there.
- **Public guides:** `../../docs/guides/` — user-facing explanations of the same subsystems live there.
