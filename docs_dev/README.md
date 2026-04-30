# Lancet2 Developer Documentation

This directory holds documentation written for **developers working on Lancet2**, not for users running it. The user-facing site (rendered by MkDocs from `docs/`) explains what Lancet2 does and how to invoke it. Everything here is for the people who write the code, debug the code, decide how the code should change, and explain the code to future contributors.

The split is deliberate. A user reading the rendered site does not need to know how the de Bruijn graph is constructed at the source level, what padding rule applies to struct members, or what was learned from a 2025-Q3 sanitizer investigation. A developer reading source files needs all three. Mixing the audiences in one tree forces every document to compromise on voice and depth. Two trees, two voices.

For an overview of the .claude/ bundle and why these docs are split into .claude/ vs docs_dev/, see ai_bundle_overview.md.

## Layout

```
docs_dev/
├── README.md                ← this file
├── style/                   ← writing and code conventions
├── subsystems/              ← deep-dives on how subsystems work
├── workflows/               ← runbooks for recurring development tasks
├── architecture/            ← cross-cutting design and architecture decisions
└── investigations/          ← postmortems, debug archaeology, perf investigations
```

Each subdirectory has its own `README.md` that states what belongs there and what doesn't. Read the relevant subdirectory README before adding a document — placement is part of the quality bar.

## What belongs in each subdirectory

**`style/`** — the project's writing and code conventions, split by audience. Four documents: how to write user-facing website docs (`website_docs.md`), how to write inline code comments (`code_comments.md`), C++ formatting and lint rules (`cpp_style.md`), and the cross-document synchronization rule plus pre-commit verification checklist (`sync_and_verification.md`). Read-mostly, change rarely. The `style/README.md` is the index plus the small set of principles that apply across all four.

**`subsystems/`** — explanatory documents about how a specific Lancet2 subsystem actually works at the source level. The current example is `probe_tracking.md` (the variant forensic pipeline). The audience is a developer about to debug or extend that subsystem; the voice is "here is what this code does and why, with enough detail that you can read the source without getting lost." Different from the public `docs/guides/` because the public guides explain user-facing behavior, not internal mechanics.

**`workflows/`** — runbooks for recurring development tasks. The current example is `profiling.md` (capture, analyze, diff a CPU profile). The audience is a developer about to do task X for the first or fifth time; the voice is "do these steps in this order, here's what you'll see, here are the failure modes." Different from subsystems because workflows tell you *how to do something*, while subsystems tell you *how something is*.

**`architecture/`** — cross-cutting design documentation. Two formats live here: numbered Architecture Decision Records (ADRs) capturing single decisions with their context and consequences, and longer overview documents explaining the shape of a system across multiple subsystems. Architecture documents are evergreen — they are kept current as the system evolves, like the style guide. See `architecture/README.md` for the writing standards.

**`investigations/`** — postmortems, deep debug archaeology, and performance investigations. These document what was understood at a specific moment in time, often during or after a hard-to-find bug or a performance puzzle. Investigations are immutable once merged — they are a snapshot of the moment, not a living document. If understanding later changes, write a new investigation that supersedes the old one rather than editing in place. See `investigations/README.md` for the writing standards.

## The rule that crosses all of `docs_dev/`

**Source code is the single source of truth.** Every claim in any document here must be verifiable against the current codebase. When the codebase changes, the relevant evergreen documents change with it. When that synchronization slips, the document is wrong, no matter how well-written it is. The full synchronization audit procedure lives in `style/sync_and_verification.md`.

The exception is `investigations/`, which document a moment in time and are not maintained against the evolving codebase by design. That exception is the entire reason `investigations/` is its own subdirectory.

## What is NOT in `docs_dev/`

- **User-facing documentation.** That lives in `docs/` and ships to the rendered MkDocs site.
- **API reference for external consumers.** Lancet2 does not currently have a stable external API; if it ever does, a `docs/api/` directory would be the right home.
- **Generated artifacts.** Profile binaries, benchmark JSON outputs, debug logs, build outputs. These are gitignored or ephemeral.
- **Drafts and scratch.** Use `notes/scratch/` for in-progress thinking. A document graduates from `notes/scratch/` to `docs_dev/` when it has been written for a reader who is not the author.

## Cross-references

A document in `docs_dev/` may link to a public `docs/` page when user-facing context helps the developer reader. The public site never links into `docs_dev/` — public links would 404 on the rendered site. This is enforced by convention, not tooling.

Cross-references between subdirectories of `docs_dev/` are expected and welcome. The four `style/` documents reference each other; subsystem deep-dives reference their related workflows; investigations reference whatever architecture documents existed at the time the incident happened.

## How to contribute

If you are adding a document:

1. Read the README of the subdirectory you intend to write in. The README states what belongs there, what voice to use, and what the quality bar is. If your document does not fit any subdirectory, either propose a new one in a separate change or reconsider whether the document is needed.
2. If you are writing in `architecture/` or `investigations/`, follow the writing standards in those subdirectory READMEs precisely. Both have stricter quality bars than the rest of `docs_dev/`.
3. If your document changes how the codebase should be read (style, conventions), update the relevant `style/` document instead of writing a new sibling.
4. Cross-reference adjacent documents where it helps the reader. Every link earns its place by saving the reader a search.

If you are editing an existing document:

1. Confirm it is not in `investigations/` — those are immutable. To revise an investigation's content, write a new investigation that supersedes it.
2. If the edit is substantive (changes a claim, not just typo fixes), check whether other documents make the same claim. The synchronization rule applies; stale claims in any document are bugs.
3. Run the pre-commit verification checklist in `style/sync_and_verification.md`.

## How this directory relates to the bundle

The Claude Code configuration bundle (`.claude/` and `CLAUDE.md`) references documents here. When a bundle agent or skill points at `docs_dev/style/cpp_style.md` for naming conventions, or at `docs_dev/subsystems/probe_tracking.md` for the probe pipeline architecture, it is treating these documents as authoritative. The bundle's `/audit-bundle` command walks several of these references for drift; see the bundle's documentation for the full list.

If a future restructure of `docs_dev/` changes paths, the bundle's path references must update in lockstep. Run `/audit-bundle` after any such restructure.
