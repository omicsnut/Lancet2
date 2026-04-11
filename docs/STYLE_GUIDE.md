# Lancet2 Developer Style Guide

This guide defines the conventions for website documentation, inline code comments, C++ code formatting, and static analysis enforcement in Lancet2. All new and edited content must follow these rules. The conventions are enforced by CI (GitHub Actions) and the local `scripts/run_clang_format.py` and `scripts/run_clang_tidy.py` helper scripts.

---

## Core Principles

These principles apply equally to website pages and code comments:

1. **Source code is the single source of truth.** Every claim must be verifiable against the current codebase. Never document behavior from memory or assumption — read the implementation first.

2. **High density, no filler.** Every sentence should teach the reader something they cannot trivially infer. Omit preambles like "In this section, we will discuss..." — start with the content.

3. **Explain the *why*, not just the *what*.** Stating that `Match = 0` is useless without explaining SIMD overflow prevention. Parameters, thresholds, and design decisions always need rationale.

4. **Quantify wherever possible.** Prefer "reduces WGS runtime by ~80%" over "significantly improves performance." Prefer "≥2 reads at the same position" over "multiple reads."

5. **Address the reader's real question.** For every feature, the reader asks: *What does it do? When should I use it? What's the trade-off?* Answer all three.

---

# Part 1 — Website Documentation (MkDocs)

## Voice & Tone

- **Active voice, present tense.** "Lancet2 skips windows with no mutation signal" — not "windows with no mutation signal are skipped by Lancet2."
- **Direct and declarative.** "This flag forces assembly of every window" — not "This flag can be used to force assembly."
- **Technical but accessible.** Define domain terms on first use, then use them freely. Example: "Cyclomatic Complexity (CC = E − V + 1): the number of independent cycles."
- **No hedging.** Write "the detector terminates immediately" — not "the detector should terminate" or "the detector is expected to terminate."
- **Second person for user actions.** "Set `--min-node-cov` to 5 if you want aggressive pruning." Use "Lancet2" (not "we" or "the system") for tool behavior.

---

## Page Structure

### Opening Paragraph

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

### Section Hierarchy

Use `##` for major sections, `###` for subsections, `####` for named sub-points (e.g., "Why Match = 0?"). Never skip heading levels (e.g., don't jump from `##` to `####`).

### Section Separators

Use `---` (horizontal rule) to separate **major conceptual boundaries** within a page (e.g., between Phase 1 and Phase 2 of the genotyping guide). Do not use `---` between every `##` heading — only where there is a genuine topic shift.

### Closing Cross-References

End each page with a `* **Read more:**` or `* **CLI reference:**` line linking to related pages. This creates a navigable web between guides.

```markdown
* **Read more:** [Alignment-Derived Annotations](alignment_annotations.md), [VCF Output Reference](vcf_output.md)
```

```markdown
* **CLI reference:** [`--no-active-region`](../reference.md#flags)
```

---

## Formatting Conventions

### Bold

Use bold for:
- **Key terms on first definition:** "builds a **colored bidirected De Bruijn graph**"
- **Critical thresholds:** "if **≥2 reads** show a mismatch at the **same** position"
- **Emphasis on surprising or important behavior:** "the entire graph is **cleared and rebuilt from scratch**"
- **Algorithmic outcomes:** "Coverage-invariant: varies **< 3%** above 60×"

Do not bold entire sentences or use bold for routine emphasis.

### Inline Code

Use backticks for:
- CLI flags and parameters: `--min-node-cov`, `-k`
- Code identifiers: `BuildSequence()`, `VariantStore`, `NormalizeVcfParsimony`
- VCF field names: `GRAPH_CX`, `NPBQ`, `SB`
- File formats and extensions: `.vcf.gz`, `.tbi`
- URI schemes: `s3://`, `gs://`
- Literal values: `0x200`, `Q20`

Do not use backticks for general technical terms (use bold instead for first mention).

### Tables

Use tables for structured comparisons of 3+ items. Always left-align columns with `:------` syntax. Common table patterns:

**Parameter tables** (value | default | rationale):
```markdown
| Parameter | Value | Standard Value | Rationale |
|:----------|:------|:---------------|:----------|
| Match | **0** | +2 | Prevents 8-bit AVX2 SIMD overflow |
```

**Comparison tables** (property | option A | option B):
```markdown
| Property | Phase 1 (MSA) | Phase 2 (Genotyping) |
|:---------|:--------------|:---------------------|
| **Strategy** | Forgiving | Strict |
```

**Pipeline stage tables** (stage | operation | purpose):
```markdown
| Stage | Operation | Purpose |
|:------|:----------|:--------|
| 0 | Low-coverage removal | Remove nodes below `--min-node-cov` |
```

### Admonitions

Use MkDocs Material admonitions (`!!!`) for information that must not be missed:

- `!!! warning` — Breaking footguns, major performance traps, data loss risks.
- `!!! note` — Deterministic behavior notes, design rationale callouts.
- `!!! tip` — Non-obvious optimizations, "why this matters" explanations.

```markdown
!!! warning "BAMs without MD tags cause ~5–10× slower runtime"
    If the input BAM/CRAM files lack the `MD` auxiliary tag, Lancet2
    **automatically disables active region detection** and assembles
    every window.
```

Use admonitions sparingly — at most 1–2 per page. If everything is a warning, nothing is.

### Numbered Lists

Use numbered lists for **sequential processes** (algorithm steps, pipeline stages). Use bullet lists for **unordered sets** (use cases, trade-offs, exclusion criteria).

### Code Blocks

Use fenced code blocks with language hints for:
- CLI commands: ` ```bash `
- Formulas and pseudocode: ` ``` ` (no language)
- VCF example records: ` ``` ` (no language — they're too wide for syntax highlighting)

### Mathematical Notation

Use inline Unicode for mathematical expressions: `O(R × L / k)`, `≥2`, `CC = E − V + 1`. For standalone formulas, use a fenced code block:

```markdown
```
combined = (global_score − local_raw_score − sc_penalty) + (local_pbq_score × local_identity)
```⁣
```

Prefer Unicode symbols (×, −, ≥, ≤, →, π) over ASCII approximations (*, -, >=, <=, ->, pi).

---

## CLI Reference Page

Each CLI parameter entry follows this template:

```markdown
#### `-x`,`--flag-name`
One-line description of what it does. Default value --> N.
Second line with behavioral detail, trade-offs, or non-obvious implications.
See [Relevant Guide](guides/relevant_guide.md) for the full explanation.
```

Rules:
- First line: what it does + default value.
- Second line: trade-offs ("higher = faster but less sensitive") and/or behavioral detail.
- Third line: cross-link to the guide that explains the underlying algorithm.
- For ranged parameters, show the allowed range: `> [MIN-MAX]. Default value --> N`

---

## Cross-Linking

### When to Link

Link to another page when:
- You mention a concept that has its own dedicated page (e.g., "See [Graph Complexity](graph_complexity.md)")
- A CLI flag has behavioral implications documented elsewhere
- A VCF FORMAT field has a detailed explanation in the annotations guide

### Link Text

- Use the page's nav title as link text: `[Active Region Detection](active_region.md)` — not `[click here](active_region.md)`.
- For CLI flags, use backtick-wrapped flag name: `` [`--extract-pairs`](../reference.md#flags) ``
- For inline "see details" within tables, use `[details](page.md#anchor)`.

### Relative Paths

- From `guides/*.md` to another guide: `(other_guide.md)` or `(other_guide.md#section)`
- From `guides/*.md` to a root page: `(../reference.md#section)`
- From root pages to guides: `(guides/guide_name.md)` or `(guides/guide_name.md#section)`

---

## Images & Diagrams

### Dark/Light Mode Support

For images with visible backgrounds, provide both a light and dark variant using Material's built-in URL fragment syntax:

```markdown
![Diagram description](../assets/diagram_name.png#only-light)
![Diagram description](../assets/diagram_name_dark.png#only-dark)
```

Material automatically shows/hides the correct image based on the user's color scheme.

### File Naming

Image assets live in `docs/assets/`. Use descriptive, lowercase, underscore-separated names:
- `pipeline_architecture.png` / `pipeline_architecture_dark.png`
- `01_dbg__chr1_38506673_38507173__low_cov_removal1__k31__comp0.png`

### Alt Text

Alt text should be a concise description of the image content: `![Lancet2 Pipeline Architecture]` — not `![image]` or `![diagram]`.

---

## Experimental Features

When documenting functionality that is implemented but not production-ready, use an admonition at the top of the relevant section:

```markdown
!!! warning "Experimental — No ML model support"
    Multi-sample and germline-only modes are functional but **experimental**.
    No pre-trained ML models are currently provided for variant filtering
    in these modes. Variant calls will require custom downstream filtering.
```

State what works, what doesn't, and what the user must do themselves.

---

## Navigation Structure

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

---

# Part 2 — Inline Code Comments

Code comments serve a **different audience** than website docs. The reader is a developer with the source file open, trying to understand *what this code does*, *why it was designed this way*, and *what the mathematical/biological intuition is* — in simple, direct terms.

## Code Comment Principles

1. **Explain the intuition, not the syntax.** The code already shows *what* is happening line by line. Comments must explain *why* this approach was chosen and *what would break* if it were done differently.

2. **Distill complexity into simple terms.** A developer reading `mPbqScore += raw * weight` doesn't need "multiply raw by weight." They need: "PBQ-weighted DP score — each base's contribution is scaled by Phred confidence `(1 − 10^(-Q/10))`. Low quality bases contribute less."

3. **Include the math, but make it legible.** Formulas are essential. Annotate each variable so the formula reads like documentation: `Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99`.

4. **State invariants and boundary conditions.** "CRITICAL: tpos coordinates are relative to alignment start, NOT position 0 of the haplotype."

5. **Keep it concise.** Every comment line must earn its place. If a 3-line comment can be condensed to 1 line without losing meaning, condense it. If a comment block exceeds ~15 lines, consider whether the explanation belongs in the website documentation instead, with a cross-reference.

---

## Comment Types

### Block Header Comments

Use `// ====...` boxed headers before major algorithmic functions. State: what the function computes, what its inputs/outputs represent, and any coordinate system or invariant the caller must respect.

```cpp
// ============================================================================
// ComputeLocalScore: evaluate alignment quality in a variant's region.
//
// Given a read→haplotype CIGAR alignment, this function extracts metrics
// for the sub-region of the haplotype that contains the variant:
//
//   1. mPbqScore:  PBQ-weighted DP score. Each position's substitution matrix
//                  contribution is scaled by (1 - ε) where ε = 10^(-PBQ/10).
//
// CRITICAL: tpos coordinates in the CIGAR are relative to the alignment start,
// NOT position 0 of the haplotype.
// ============================================================================
```

### ASCII Art Diagrams

Use ASCII diagrams when spatial relationships are non-trivial — graph topology, coordinate systems, and data structure layouts:

```cpp
//   Haplotype Array :  [ A  B  C  D  E  F  G  H  I ]
//   Variant Region  :           [var_start ... var_end)
//   Alignment       :     [aln_start ... tpos_rel ... ]
```

```cpp
//                          .-->  (T)[3]  --.
//                         /                 \
//   Anchor: (A)[2] ------+                   +-----> Target: (G)[5]
//                         \                 /
//                          `-->  (C)[4]  --'
```

### Pipeline / Architecture Comments

Use `///` Doxygen-style comments with box-drawing characters for multi-stage pipeline annotations before major entry points:

```cpp
/// Pipeline architecture for haplotype assembly:
///
///  ┌─────────────┐
///  │ Outer loop: │  Iterate k from min_k to max_k in steps of mKmerStepLen.
///  │ k-value scan│  If haplotypes are found at any k, stop.
///  └──────┬──────┘
///         │
///         ▼
///  ┌─────────────┐
///  │  BuildGraph │  Build the bidirected de Bruijn graph from reads + reference.
///  └─────────────┘
```

### Inline Formula Comments

Place formula annotations on the line above or next to the computation. Include representative values:

```cpp
/// Convert a Phred quality score to confidence weight: 1 - 10^(-Q/10).
/// Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99, Q=30 → 0.999, Q=40 → 0.9999
inline auto PhredToConfidence(u8 const qual) -> f64 {
```

### Data Structure Visual Tables

Use visual table comments for lookup tables, scoring matrices, and constant arrays:

```cpp
// ┌───────────────────────────────────────────────┐
// │ 5×5 Scoring Matrix for ComputeLocalScore      │
// │ Target (R) × Query (C) | A=0 C=1 G=2 T=3 N=4  │
// ├───────┬───────┬───────┬───────┬───────┬───────┤
// │       │  A(0) │  C(1) │  G(2) │  T(3) │  N(4) │
// │  A(0) │    1  │   -4  │   -4  │   -4  │    0  │
// └───────┴───────┴───────┴───────┴───────┴───────┘
```

### Design Decision Comments

Begin with rationale language and explain what the alternative was:

```cpp
// Skip this k if the reference itself has a repeated k-mer — the de Bruijn
// graph would contain a cycle by construction, making assembly pointless.
```

```cpp
// Quartile CV = (Q3 - Q1) / (Q3 + Q1)
// Robust against outliers unlike standard CV (stddev / mean).
```

### Member Variable Comments

Document VCF field mapping, memory layout, and data structure choice rationale inline:

```cpp
absl::flat_hash_map<std::string, std::vector<usize>> mAltAllelesToHaps;  // Maps ALT allele → haplotype indices
std::string mRefAllele;                                                  // 24B
usize mGenomeStartPos = SIZE_MAX;                                        // 8B
```

---

## Header Files vs Implementation Files

| Location | Comment Style | Content |
|:---------|:-------------|:--------|
| **Header (`.h`)** | `///` Doxygen + `//` blocks above class | Public API contract: what the class/method does, what invariants it maintains, FORMAT field catalog |
| **Implementation (`.cpp`)** | `// ====` boxed headers + inline `//` | Algorithm walkthrough: the *how* and *why*, ASCII diagrams, formula derivations, coordinate system docs |

---

## Constants and Magic Numbers

Every named constant should have a comment explaining:
1. What biological or computational concept it represents
2. Why this specific value was chosen
3. What changes if you modify it

```cpp
static constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1'048'576;  // 2^20 — caps BFS walk-tree expansion
```

If a literal appears directly in code (not via a named constant), it **must** have an inline comment. If it appears more than once, extract it to a named `constexpr`.

---

## What NOT to Comment

- **Don't restate the code.** `++mAligned; // increment aligned count` teaches nothing.
- **Don't comment obvious control flow.** `if (x > 0) // if x is positive` is noise.
- **Don't leave TODO/FIXME without context.** If a TODO stays, it needs: what needs to change, why, and what blocks it.
- **Don't use filler adverbs.** "mathematically guarantees strictly progressively monotonically upward" is 5 adverbs for a sort. Write: "Sort variants by position for ordered VCF output."
- **Don't write essays.** If a comment block exceeds ~15 lines, the explanation likely belongs in the website docs. Add a cross-reference: `// See docs/guides/variant_discovery_genotyping.md for the full scoring derivation.`

---

# Part 3 — Synchronization & Verification

## Synchronization Rule

Every sentence in code comments and every sentence in website documentation **must** be in sync with the actual logic in the codebase. Stale comments are worse than no comments.

When any refactor changes a field name, metric computation, or behavioral semantic, perform an exhaustive codebase-wide audit by grepping for the old names across:
- `src/**/*.h` and `src/**/*.cpp` — code comments and docstrings
- `docs/**/*.md` — website documentation and tables
- VCF header definitions in `pipeline_runner.cpp`
- FORMAT string comments in `variant_call.cpp` and `variant_call.h`

**Common patterns that go stale:**
- Architecture diagram comments referencing VCF FORMAT field lists
- Metric descriptions in tangentially related functions
- `ReadEvidence` struct member comments referencing VCF field names
- Cross-references in one doc page to a metric described in another

---

# Part 4 — C++ Code Style & Linting

Lancet2 enforces code style via **clang-format** (formatting) and **clang-tidy** (static analysis). Both run as CI checks on every push and pull request via `.github/workflows/lint_cpp.yml`. Violations block merging.

## Configuration Files

| Tool | Config |
|:-----|:-------|
| clang-format | `.clang-format` |
| clang-tidy | `.clang-tidy` |

Both configs are extensively commented with rationale for every setting. Read them before modifying.

---

## Formatting Rules (clang-format)

The key formatting decisions that most affect daily coding:

### Line Width & Indentation
- **100-column limit.** Wide enough for descriptive names and templates, narrow enough for side-by-side diff review.
- **2-space indent** for blocks, **4-space continuation** for wrapped args/expressions.
- **Spaces only** — no tabs.

### East Const
Lancet2 enforces **east const** placement: the `const` qualifier appears *after* the type.

```cpp
// ✅ Correct (east const)
int const x = 42;
auto const& ref = container;
std::string_view const seq = read.Sequence();

// ❌ Wrong (west const)
const int x = 42;
const auto& ref = container;
```

The full qualifier order enforced by clang-format is: `inline` → `static` → `constexpr` → *type* → `const` → `volatile` → `restrict`.

### Pointer & Reference Alignment
Left-aligned (Google style): `Type* ptr`, `Type& ref` — not `Type *ptr` or `Type * ptr`.

### Include Ordering
Includes are automatically sorted and grouped into 4 tiers (separated by blank lines):

```cpp
// Tier 0: Main header (auto-detected from filename)
#include "lancet/cbdg/graph.h"

// Tier 1: Project headers
#include "lancet/base/types.h"
#include "lancet/cbdg/kmer.h"

// Tier 2: Third-party headers
#include "absl/container/flat_hash_map.h"
#include "spdlog/fmt/bundled/core.h"

// Tier 3: C++ standard library
#include <algorithm>
#include <vector>

// Tier 4: C standard library wrappers
#include <cmath>
#include <cstdint>
```

Do **not** manually organize includes — let clang-format handle the grouping. If `// clang-format off` is needed to protect include order (e.g., `extern "C"` blocks), document why.

### Quote vs Angle-Bracket Includes
Angle-bracket includes (`<>`) are reserved **exclusively** for C and C++ standard library headers. All other headers — project headers, third-party libraries (Abseil, spdlog, spoa, htslib, etc.) — use quote includes (`""`).

```cpp
// ✅ Correct
#include "absl/container/flat_hash_map.h"  // third-party → quotes
#include "spdlog/fmt/bundled/core.h"       // third-party → quotes
#include <algorithm>                       // C++ stdlib → angle brackets
#include <cmath>                           // C stdlib → angle brackets

// ❌ Wrong
#include <absl/container/flat_hash_map.h>  // third-party must NOT use angle brackets
```

This convention is not enforced by clang-format or clang-tidy — it must be followed manually during code review.

### Brace Style
K&R "Attach" — opening brace on the same line as the statement. Empty bodies compacted to `{}`.

### Integer Literal Separators
Large literals get `'` separators for readability: `1'048'576` (decimal, groups of 3), `0xFF'00'FF` (hex, groups of 2), `0b1111'0000` (binary, groups of 4). Only applied to literals with ≥5 digits.

### Comment Reflowing
clang-format **reflows long comments** to fit within the 100-column limit. This interacts directly with code comment conventions — keep comment lines concise so reflowing doesn't produce awkward breaks. Use `// clang-format off` / `// clang-format on` around ASCII art diagrams and visual tables to prevent reflowing.

---

## Static Analysis (clang-tidy)

### Naming Conventions

The naming convention is enforced by `readability-identifier-naming` checks:

| Element | Convention | Example |
|:--------|:-----------|:--------|
| Classes / Structs | `PascalCase` | `VariantCall`, `AsyncWorker` |
| Methods / Functions | `PascalCase` | `BuildSequence`, `FindEdges` |
| Member Variables | `mPascalCase` | `mKmer`, `mEdges`, `mCurrK` |
| Local Variables | `snake_case` | `new_edge`, `comp_id` |
| Parameters | `snake_case` | `stop_token`, `window_length` |
| Constants (`constexpr` / `static const`) | `UPPER_CASE` | `DEFAULT_MIN_KMER_LEN` |
| Enum Values | `UPPER_CASE` | `PLUS_PLUS`, `FWD`, `UNKNOWN` |
| Type Aliases | `PascalCase` | `NodeID`, `EdgeList` |
| Template Parameters | `PascalCase` | `Args`, `TagResultValue` |
| Namespaces | `snake_case` | `lancet`, `cbdg`, `hts` |
| Macros | `UPPER_CASE` | `LANCET_DEVELOP_MODE` |

**Ignored patterns:** Type aliases `i8`, `u32`, `f64`, `usize` are exempted. Include guard macros (`SRC_LANCET_*_H_`) are exempted. STL iterator methods (`begin`, `end`, `cbegin`, `cend`) are exempted.

### Check Categories

The following check families are enabled (with curated exclusions documented in `.clang-tidy`):

| Category | Coverage | Key Exclusions |
|:---------|:---------|:---------------|
| `bugprone-*` | Bug detection | `-easily-swappable-parameters` (too noisy for graph APIs) |
| `cert-*` | CERT secure coding | `-err58-cpp` (static constexpr never throws) |
| `cppcoreguidelines-*` | Core Guidelines | FFI checks disabled (htslib/minimap2 require raw arrays, varargs) |
| `modernize-*` | C++20 modernization | `-use-trailing-return-type` (breaks lambda deduction in GCC), `-use-ranges` (breaks types without `<=>`) |
| `performance-*` | Performance | All enabled; `string_view` and `Span` exempted from value-param check |
| `readability-*` | Readability & naming | `-magic-numbers` (scientific constants), `-redundant-member-init` (explicit init is convention) |
| `misc-*` | Miscellaneous | `-no-recursion` (graph DFS/BFS is inherently recursive), `-include-cleaner` (breaks include order) |

### Warnings Promoted to Errors

These checks are hard CI errors — they will fail the build:
- All `bugprone-*` (except `exception-escape`)
- All `cert-*`, `performance-*`, `portability-*`
- `modernize-use-override`, `modernize-use-nullptr`, `modernize-use-emplace`
- `readability-container-contains`, `readability-container-size-empty`
- `misc-const-correctness`

### Function Complexity Thresholds

| Metric | Threshold |
|:-------|:----------|
| Cognitive complexity | 35 |
| Statement count | 150 |
| Line count | 200 |
| Parameter count | 6 |
| Nesting depth | 4 |

Functions exceeding these thresholds should be decomposed. If a function legitimately needs to exceed a threshold (e.g., a CIGAR walk), use a scoped `// NOLINTNEXTLINE(readability-function-cognitive-complexity)` with a rationale comment.

### NOLINT Usage Rules

Suppression comments are acceptable **only** when:
1. The code is correct and the check is a false positive (e.g., C FFI interop)
2. The suppression uses `NOLINTNEXTLINE(check-name)` or a `NOLINTBEGIN(check-name)` / `NOLINTEND(check-name)` block — **never** inline `NOLINT` (it pollutes the statement line) and **never** bare/unscoped (omitting the check name)
3. A rationale comment explains *why* the suppression is needed

```cpp
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) -- htslib C FFI requires raw cast
auto const* aptr = reinterpret_cast<unsigned long long const*>(first.data());
```

`NOLINTBEGIN` / `NOLINTEND` blocks are acceptable for groups of related suppressions (e.g., struct members that must be `const&`):

```cpp
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
absl::Span<u8 const> const mQuery;         // 16B
absl::Span<u8 const> const mTarget;        // 16B
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
```

---

## Running Locally

### Formatting
```bash
# Check (dry-run) — reports files that need formatting
python3 scripts/run_clang_format.py

# Fix — applies formatting in-place
python3 scripts/run_clang_format.py --fix

# Check a specific directory or file
python3 scripts/run_clang_format.py src/lancet/cbdg
python3 scripts/run_clang_format.py --fix src/lancet/caller/genotyper.cpp
```

### Static Analysis
```bash
# Requires a build directory with compile_commands.json
python3 scripts/run_clang_tidy.py

# With auto-fix (review changes carefully)
python3 scripts/run_clang_tidy.py --fix

# Custom build directory
python3 scripts/run_clang_tidy.py --build-dir cmake-build-release
```

### CI Behavior
The `lint_cpp.yml` workflow runs both checks on every push and pull request:
1. **clang-format** — dry-run check on all `src/lancet/`, `tests/`, `benchmarks/` files
2. **clang-tidy** — full build + analysis with `compile_commands.json`

Concurrency guards (`lint-${{ github.ref }}`) cancel in-progress runs when new commits are pushed.

---

# Part 5 — Verification Checklist

Before merging any change, verify:

### Website Documentation
- [ ] `pixi run mkdocs build --strict` passes with zero warnings
- [ ] All internal links resolve (no "anchor not found" messages)
- [ ] New pages are added to the `nav:` block in `mkdocs.yml`
- [ ] All numerical values and thresholds are verified against the current source code
- [ ] Cross-links exist between the new page and related existing pages
- [ ] CLI flags mentioned in guides link back to the CLI Reference page
- [ ] Dark mode renders correctly for any new image assets

### Code Comments
- [ ] Every public method has a doc comment in the header
- [ ] Complex algorithms have step-by-step inline comments with ASCII diagrams where helpful
- [ ] All formulas have variable annotations (e.g., `Q=0 → 0.0, Q=20 → 0.99`)
- [ ] No comments that restate the code — every comment explains *why* or *what the intuition is*
- [ ] Constants have rationale comments explaining the chosen value

### Code Style & Linting
- [ ] `python3 scripts/run_clang_format.py` reports zero formatting violations
- [ ] `python3 scripts/run_clang_tidy.py` reports zero warnings
- [ ] All new identifiers follow the naming convention table
- [ ] East const is used throughout (`int const x`, not `const int x`)
- [ ] All NOLINT suppressions are scoped and have rationale comments
- [ ] No bare `NOLINT` — always specify the check name

### Sync
- [ ] VCF header in `pipeline_runner.cpp` matches website documentation
- [ ] FORMAT field comments in `variant_call.cpp` match actual output
- [ ] Codebase-wide grep for any renamed fields — zero stale references remain
