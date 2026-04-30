# Inline Code Comment Style

The audience for inline code comments is **a developer with the source file open** — likely an experienced C++ programmer, possibly a future you returning to a function in six months. They can read the code; what they need is the intuition, the design rationale, and the boundary conditions the code itself cannot communicate.

This is a fundamentally different audience from the website documentation reader. Where website docs explain user-facing behaviour to a clinician or methods reviewer, code comments explain implementation intuition to a developer. Both share the project's Core Principles in `style/README.md`, but the voice and bar differ.

## Code comment principles

1. **Explain the intuition, not the syntax.** The code already shows *what* is happening line by line. Comments must explain *why* this approach was chosen and *what would break* if it were done differently.

2. **Distill complexity into simple terms.** A developer reading `mPbqScore += raw * weight` doesn't need "multiply raw by weight." They need: *"PBQ-weighted DP score — each base's contribution is scaled by Phred confidence `(1 − 10^(−Q/10))`. Low-quality bases contribute less."*

3. **Include the math, but make it legible.** Formulas are essential where the implementation depends on them. Annotate each variable so the formula reads like documentation, and include representative values: `Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99`.

4. **State invariants and boundary conditions.** *"CRITICAL: tpos coordinates are relative to alignment start, NOT position 0 of the haplotype."* The code that maintains the invariant is verifiable; the invariant itself usually is not — it lives in the comment.

5. **Keep it concise.** Every comment line must earn its place. If a 3-line comment can be condensed to 1 line without losing meaning, condense it. If a comment block exceeds ~15 lines, consider whether the explanation belongs in the website documentation or a `subsystems/` deep-dive instead, with a cross-reference from the source.

6. **No unexplained jargon — calibrated for a developer audience.** Three categories:

   - **Replace.** Terms that add no precision over plain language. Use "independent" not "orthogonal", "highly correlated" not "collinear", "depth-dependent precision" not "heteroscedastic".
   - **Explain then use.** Terms that are precise and useful but not universally known. Define on first occurrence with a brief parenthetical, then use freely: *"Bessel's correction (divides by n−1 instead of n to avoid underestimating spread from a sample)."*
   - **Use freely.** Terms a developer audience knows by default — `Euclidean`, `log10`, `branchless`, `SIMD`, `cache line`, `LSB`/`MSB`, `epsilon`, `NaN`. The bar is *would another C++ developer have to look it up?*, not *would a clinician have to look it up?*. The website-docs bar is stricter for the same reason: a different audience.

   Where the same term is used in code comments and website docs, the code-comment version may treat it as known while the website version explains it. That asymmetry is intentional.

## Comment types

Six patterns cover ~90% of comment-writing situations in the project. Use them.

### Block header comments

Use `// ====...` boxed headers before major algorithmic functions. State what the function computes, what its inputs and outputs represent, and any coordinate system or invariant the caller must respect.

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

### ASCII art diagrams

Use ASCII diagrams when spatial relationships are non-trivial — graph topology, coordinate systems, data structure layouts:

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

ASCII diagrams must be wrapped in `// clang-format off` / `// clang-format on` to prevent reflowing. The same applies to data-structure visual tables.

### Pipeline / architecture comments

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

### Inline formula comments

Place formula annotations on the line above or next to the computation. Include representative values:

```cpp
/// Convert a Phred quality score to confidence weight: 1 - 10^(-Q/10).
/// Q=0 → 0.0, Q=10 → 0.9, Q=20 → 0.99, Q=30 → 0.999, Q=40 → 0.9999
inline auto PhredToConfidence(u8 const qual) -> f64 {
```

### Data-structure visual tables

Use visual table comments for lookup tables, scoring matrices, and constant arrays. Wrap in `// clang-format off` / `// clang-format on`:

```cpp
// ┌───────────────────────────────────────────────┐
// │ 5×5 Scoring Matrix for ComputeLocalScore      │
// │ Target (R) × Query (C) | A=0 C=1 G=2 T=3 N=4  │
// ├───────┬───────┬───────┬───────┬───────┬───────┤
// │       │  A(0) │  C(1) │  G(2) │  T(3) │  N(4) │
// │  A(0) │    1  │   -4  │   -4  │   -4  │    0  │
// └───────┴───────┴───────┴───────┴───────┴───────┘
```

### Design decision comments

Begin with rationale language and explain what the alternative was. The comment exists because the choice is non-obvious; if the choice were obvious, the comment would not be needed:

```cpp
// Skip this k if the reference itself has a repeated k-mer — the de Bruijn
// graph would contain a cycle by construction, making assembly pointless.
```

```cpp
// Quartile CV = (Q3 - Q1) / (Q3 + Q1)
// Robust against outliers unlike standard CV (stddev / mean).
```

## Member variable comments

Document VCF field mapping, memory layout, and data structure choice rationale inline on the declaration line:

```cpp
absl::flat_hash_map<std::string, std::vector<usize>> mAltAllelesToHaps;  // Maps ALT allele → haplotype indices
std::string mRefAllele;                                                  // 24B
usize mGenomeStartPos = SIZE_MAX;                                        // 8B
```

The size annotation (`// 8B`, `// 24B`, etc.) serves double duty: it documents the field and makes alignment compliance visible at a glance. The padding-and-alignment rules that govern *which order* members must be declared in are in `cpp_style.md` § "Struct/Class Memory Layout"; this section covers only the comment convention.

## Header files vs implementation files

| Location | Comment style | Content |
|:---------|:--------------|:--------|
| **Header (`.h`)** | `///` Doxygen + `//` blocks above the class | Public API contract: what the class/method does, what invariants it maintains, FORMAT field catalog |
| **Implementation (`.cpp`)** | `// ====` boxed headers + inline `//` | Algorithm walkthrough: the *how* and *why*, ASCII diagrams, formula derivations, coordinate system docs |

The header carries the contract. The implementation carries the walkthrough. A reader who only opens the header should be able to use the class. A reader who opens the implementation should be able to understand or modify it.

## Constants and magic numbers

Every named constant should have a comment explaining:

1. What biological or computational concept it represents.
2. Why this specific value was chosen.
3. What changes if you modify it.

```cpp
static constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1'048'576;  // 2^20 — caps BFS walk-tree expansion
```

If a literal appears directly in code (not via a named constant), it **must** have an inline comment. If it appears more than once, extract it to a named `constexpr`.

## What NOT to comment

- **Don't restate the code.** `++mAligned; // increment aligned count` teaches nothing.
- **Don't comment obvious control flow.** `if (x > 0) // if x is positive` is noise.
- **Don't leave TODO/FIXME without context.** If a TODO stays in the source, it must say what needs to change, why, and what blocks it. A bare `// TODO: fix this` is technical debt that nobody can resolve.
- **Don't use filler adverbs.** "Mathematically guarantees strictly progressively monotonically upward" is five adverbs for a sort. Write: *"Sort variants by position for ordered VCF output."*
- **Don't write essays in source.** If a comment block exceeds ~15 lines, the explanation belongs in `docs_dev/subsystems/` or `docs/guides/` (depending on audience). Cross-reference from the source: `// See docs_dev/subsystems/probe_tracking.md for the full design.`
- **Don't use jargon when a plain word works.** If "independent" conveys the same meaning as "orthogonal", use "independent." If a technical term is essential, define it on first use.

## Cross-references to other style documents

- The synchronization rule for keeping comments aligned with code as the code changes is in `sync_and_verification.md`. Apply it after any rename, metric change, or behavioural-semantics edit.
- The clang-format reflowing behaviour that affects long comment lines is documented in `cpp_style.md` § "Comment Reflowing".
- The comment-in-code-vs-prose-in-website-docs decision (when does an explanation belong in source vs in `docs/`) is in `website_docs.md` § "Page Structure" and in this document's "What NOT to comment" item about essays.
