# Website Documentation Style

The audience for content in `docs/` is a **Lancet2 user reading the rendered MkDocs site** — bioinformatician, clinician, methods reviewer, or downstream tool developer. They have no source code open. They want answers, not implementation. Every voice and structural rule in this document derives from that.

The companion document for inline code comments is `code_comments.md`. Where the same concept appears in both — the no-jargon rule, the explain-why-not-just-what principle — the rules differ in the bar applied for the specific audience. Do not assume what works for one works for the other.

## Voice and tone

- **Active voice, present tense.** "Lancet2 skips windows with no mutation signal." Not "windows with no mutation signal are skipped by Lancet2."

- **Direct and declarative.** "This flag forces assembly of every window." Not "This flag can be used to force assembly."

- **Technical but accessible.** Define domain terms on first use, then use them freely. Example: *"Cyclomatic Complexity (CC = E − V + 1): the number of independent cycles."* The first occurrence carries the cost of definition; subsequent occurrences read naturally for the technical reader.

- **No hedging.** "The detector terminates immediately." Not "the detector should terminate" or "the detector is expected to terminate." If the behaviour is uncertain or version-dependent, state that explicitly: *"In versions before 2.9.0, this behaviour differed — see CHANGELOG."* Hedging without a stated reason confuses the reader about whether you are describing the system or wishing for it.

- **No unexplained jargon.** Technical terms fall into three categories:
  - **Replace.** Terms that add no precision over plain language. Use "independent" not "orthogonal", "highly correlated" not "collinear", "depth-dependent precision" not "heteroscedastic."
  - **Explain then use.** Terms that are precise and useful but not universally known. Define on first occurrence with a brief parenthetical, then use freely: *"Bessel's correction (divides by n−1 instead of n to avoid underestimating spread from a sample)."*
  - **Use freely.** Terms already explained in context, or self-evident: *"alignment", "coverage", "VCF FORMAT field"*.

  The test for any term: *would a biologist or clinician reading this for the first time have to stop and look it up?* If yes, either explain it or replace it. The bar is stricter here than in code comments because the website audience includes non-developers.

- **Second person for user actions, "Lancet2" for tool behaviour.** *"Set `--min-node-cov` to 5 if you want aggressive pruning."* — *"Lancet2 skips this stage when MD tags are absent."* Avoid "we" and "the system" — the former is presumptuous, the latter is vague.

## Page structure

### Opening paragraph

Every page starts with a single paragraph (2–3 sentences) that tells the reader exactly what this page covers and why it matters. No heading before this paragraph — it sits directly under the `# Title`.

**Good:**

```markdown
# Active Region Detection

Before assembling the de Bruijn graph for a window, Lancet2 runs a fast
pre-scan to check whether the region contains any evidence of variation.
Windows with no mutation signal are skipped entirely, reducing WGS runtime
by ~80%.
```

**Bad:**

```markdown
# Active Region Detection

## Introduction

This page describes the active region detection feature of Lancet2.
Active region detection is an important optimization that helps improve
the performance of the variant calling pipeline.
```

The bad version says nothing the title did not already imply.

### Section hierarchy

Use `##` for major sections, `###` for subsections, `####` for named sub-points (e.g., "Why Match = 0?"). Never skip heading levels — do not jump from `##` to `####`.

### Section separators

Use `---` (horizontal rule) to separate major conceptual boundaries within a page (e.g., between Phase 1 and Phase 2 of the genotyping guide). Do not use `---` between every `##` heading — only where there is a genuine topic shift.

### Closing cross-references

End each page with a `* **Read more:**` or `* **CLI reference:**` line linking to related pages. This creates a navigable web between guides.

```markdown
* **Read more:** [Alignment-Derived Annotations](alignment_annotations.md), [VCF Output Reference](vcf_output.md)
* **CLI reference:** [`--no-active-region`](../reference.md#flags)
```

## Formatting conventions

### Bold

Use bold for:

- **Key terms on first definition:** *"builds a **colored bidirected De Bruijn graph**"*
- **Critical thresholds:** *"if **≥2 reads** show a mismatch at the **same** position"*
- **Emphasis on surprising or important behaviour:** *"the entire graph is **cleared and rebuilt from scratch**"*
- **Algorithmic outcomes:** *"Coverage-invariant: varies **< 3%** above 60×"*

Do not bold entire sentences. Do not use bold for routine emphasis.

### Inline code

Use backticks for:

- CLI flags and parameters: `--min-node-cov`, `-k`
- Code identifiers: `BuildSequence()`, `VariantStore`, `NormalizeVcfParsimony`
- VCF field names: `GRAPH_CX`, `NPBQ`, `SB`
- File formats and extensions: `.vcf.gz`, `.tbi`
- URI schemes: `s3://`, `gs://`
- Literal values: `0x200`, `Q20`

Do not use backticks for general technical terms — use bold instead for first mention of a concept.

### Tables

Use tables for structured comparisons of three or more items. Always left-align columns with `:------` syntax. Three common patterns:

**Parameter tables** (value | default | rationale):

```markdown
| Parameter | Value | Standard Value | Rationale |
|:----------|:------|:---------------|:----------|
| Match | **0** | +2 | Keeps SPOA in the faster int16 SIMD path |
```

**Comparison tables** (property | option A | option B):

```markdown
| Property | Phase 1 (MSA) | Phase 2 (Genotyping) |
|:---------|:--------------|:---------------------|
| **Strategy** | Forgiving | Strict |
```

**Pipeline-stage tables** (stage | operation | purpose):

```markdown
| Stage | Operation | Purpose |
|:------|:----------|:--------|
| 0 | Low-coverage removal | Remove nodes below `--min-node-cov` |
```

### Admonitions

Use MkDocs Material admonitions (`!!!`) for information that must not be missed:

- `!!! warning` — breaking footguns, major performance traps, data loss risks.
- `!!! note` — deterministic behaviour notes, design rationale callouts.
- `!!! tip` — non-obvious optimizations, "why this matters" explanations.

```markdown
!!! warning "BAMs without MD tags cause ~5–10× slower runtime"
    If the input BAM/CRAM files lack the `MD` auxiliary tag, Lancet2
    **automatically disables active region detection** and assembles
    every window.
```

Use admonitions sparingly — at most one or two per page. If everything is a warning, nothing is.

### Numbered lists

Use numbered lists for sequential processes (algorithm steps, pipeline stages). Use bullet lists for unordered sets (use cases, trade-offs, exclusion criteria).

### Code blocks

Use fenced code blocks with language hints for:

- CLI commands: ` ```bash `
- Formulas and pseudocode: ` ``` ` (no language)
- VCF example records: ` ``` ` (no language — they're too wide for syntax highlighting)

### Mathematical notation

Use inline Unicode for mathematical expressions: `O(R × L / k)`, `≥2`, `CC = E − V + 1`. For standalone formulas, use a fenced code block:

````markdown
```
combined = (global_score − local_raw_score − sc_penalty) + (local_pbq_score × local_identity)
```
````

Prefer Unicode symbols (×, −, ≥, ≤, →, π) over ASCII approximations (`*`, `-`, `>=`, `<=`, `->`, `pi`).

## CLI reference page

Each CLI parameter entry follows this template:

```markdown
#### `-x`,`--flag-name`
One-line description of what it does. Default value --> N.
Second line with behavioral detail, trade-offs, or non-obvious implications.
See [Relevant Guide](guides/relevant_guide.md) for the full explanation.
```

Rules:

- First line: what it does + default value.
- Second line: trade-offs ("higher = faster but less sensitive") and/or behavioural detail.
- Third line: cross-link to the guide that explains the underlying algorithm.
- For ranged parameters, show the allowed range: `> [MIN-MAX]. Default value --> N`.

## Cross-linking

### When to link

Link to another page when:

- You mention a concept that has its own dedicated page (e.g., *"See [Graph Complexity](guides/graph_complexity.md)"*).
- A CLI flag has behavioural implications documented elsewhere.
- A VCF FORMAT field has a detailed explanation in the annotations guide.

### Link text

- Use the page's nav title as link text: `[Active Region Detection](active_region.md)` — not `[click here](active_region.md)`.
- For CLI flags, use the backtick-wrapped flag name: `` [`--extract-pairs`](../reference.md#flags) ``.
- For inline "see details" within tables, use `[details](page.md#anchor)`.

### Relative paths

- From `guides/*.md` to another guide: `(other_guide.md)` or `(other_guide.md#section)`.
- From `guides/*.md` to a root page: `(../reference.md#section)`.
- From root pages to guides: `(guides/guide_name.md)` or `(guides/guide_name.md#section)`.

## Images and diagrams

### Dark/light mode support

For images with visible backgrounds, provide both a light and dark variant using Material's URL fragment syntax:

```markdown
![Diagram description](../assets/diagram_name.png#only-light)
![Diagram description](../assets/diagram_name_dark.png#only-dark)
```

Material automatically shows or hides the correct image based on the user's colour scheme.

### File naming

Image assets live in `docs/assets/`. Use descriptive, lowercase, underscore-separated names:

- `pipeline_architecture.png` / `pipeline_architecture_dark.png`
- `01_dbg__chr1_38506673_38507173__low_cov_removal1__k31__comp0.png`

### Alt text

Alt text should be a concise description of the image content: `![Lancet2 Pipeline Architecture]` — not `![image]` or `![diagram]`.

## Experimental features

When documenting functionality that is implemented but not production-ready, use an admonition at the top of the relevant section:

```markdown
!!! warning "Experimental — No ML model support"
    Multi-sample and germline-only modes are functional but **experimental**.
    No pre-trained ML models are currently provided for variant filtering
    in these modes. Variant calls will require custom downstream filtering.
```

State what works, what doesn't, and what the user must do themselves.

## Navigation structure

Pages are organized into four groups in `mkdocs.yml`:

| Group | Purpose | Audience |
|:------|:--------|:---------|
| **Usage Guides** | How to run Lancet2 on your data | New users |
| **How It Works** | Algorithmic internals and design | Developers, methods reviewers |
| **Annotations** | VCF output fields and ML features | Bioinformaticians, ML engineers |
| **CLI Reference** | Exhaustive flag/parameter reference | All users |

When adding a new page:

1. Determine which group it belongs to based on the table above.
2. Within the group, order pages by **pipeline execution order** (How It Works) or **dependency order** (Annotations: VCF Output first, then individual field deep-dives).
3. Add cross-links from related pages and the CLI Reference.
