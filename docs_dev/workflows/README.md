# Workflows

This directory holds **runbooks for recurring development tasks**. Audience: a developer about to do task X for the first or fifth time. Voice: procedural — the steps, in order, with the failure modes and the verification checks that confirm each step worked.

## What belongs here

A document earns its place in `workflows/` when:

- It walks through a repeatable development task end-to-end (CPU profiling, sanitizer iteration, capturing a benchmark baseline, cutting a release).
- The reader's question is *"how do I do X?"* — not *"what is X?"* (subsystems) and not *"why is the system shaped this way?"* (architecture).
- The procedure is involved enough that re-deriving it from scratch wastes time and reintroduces bugs.
- The document will be kept current as the underlying tools evolve.

If the procedure is one or two commands, it doesn't need a workflow document — a slash command or a short note in CLAUDE.md is enough. Workflows earn their place when there are real failure modes, real ordering constraints, and real verification steps to document.

## Voice and structure

Most workflow documents follow this shape:

1. **What this workflow accomplishes** — one paragraph orientation. What you'll have at the end.
2. **Prerequisites** — what must be true before starting (specific build, specific data, specific environment).
3. **Step-by-step procedure** — the actual steps, in order. Each step gets a heading. Each step states what to do, what to expect, and how to verify it worked.
4. **Common failure modes** — when step 3 fails, this is what it looks like and how to fix it. This section is often the most-read part of the document.
5. **What to do with the result** — where the output goes, how it's used downstream, when to discard it.

Cite specific commands, specific file paths, specific expected outputs. *"Run `pixi run -e profiling analyze-profile path/to/profile.bin -- --view top --top 50`"* — not *"analyze the profile."*

The existing `profiling.md` is the prototype. New workflow documents should match its specificity and verification discipline.

## Lifecycle

Workflow documents are **evergreen**. They must stay current as the underlying tools change. When a `pixi.toml` task name changes, when a script's CLI gains a flag, when a build path moves, the workflow that depends on it changes in lockstep.

The synchronization rule (`../style/sync_and_verification.md`) applies. The bundle's `/audit-bundle` command walks several workflow paths against actual source-of-truth files (e.g., pixi task names) to catch drift.

## File naming

Lowercase, underscore-separated, named for the task:

- `profiling.md`
- `sanitizer_iteration.md`
- `release_cutting.md`
- `benchmark_baseline_capture.md`

# Existing workflow documents

- **[CPU profiling](profiling.md)** — gperftools profile capture against a real workload, analysis via the project's `analyze_profile.py` wrapper, diffing against a baseline.

# Cross-references

- **Subsystems:** `../subsystems/` — explanatory documents about *what a subsystem is* live there.
- **Architecture:** `../architecture/` — design overviews and decisions live there.
- **Style:** `../style/cpp_style.md` § "Running locally" — short list of the most-used pixi tasks for build, test, lint, format. The workflow documents here go deeper than that summary.
