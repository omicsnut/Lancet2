# Path-scoped rules

Each `.md` file in this directory carries layer-specific design logic
that auto-loads only when Claude is editing files matching the
`paths:` glob in the file's YAML frontmatter. The mechanism is the
path-scoped-rules convention: files here cost zero baseline tokens and
load on-demand only when a tool call touches a matching path.

## What lives here

Six rule files, one per `src/lancet/<layer>/`:

| File          | Auto-loads when editing      | Focus |
|:--------------|:-----------------------------|:------|
| `base.md`     | `src/lancet/base/**`         | Numerical stability, SIMD, statistical tests, type aliases |
| `hts.md`      | `src/lancet/hts/**`          | HTSlib RAII, the `bam1_t` lifetime contract, BGZF streaming |
| `cbdg.md`     | `src/lancet/cbdg/**`         | Bidirected colored de Bruijn graph, k-mer canonicalization, walk enumeration |
| `caller.md`   | `src/lancet/caller/**`       | SPOA scoring, Dirichlet-Multinomial GLs, allele assignment math |
| `core.md`     | `src/lancet/core/**`         | Pipeline orchestration, sharded `VariantStore`, chunked sorted VCF flush |
| `cli.md`      | `src/lancet/cli/**`          | CLI11 surface, VCF header building, BGZF stream lifecycle |

## What does NOT live here

These rules carry layer-specific *design logic*. Code-style rules
(naming, member layout, NOLINT discipline, function-complexity
ceilings) live in `AGENTS.md` and apply uniformly to every layer.
Hook-enforced rules (layer direction, naming, protected paths) are in
`.claude/hooks/`. VCF schema invariants live in
`.claude/agents/vcf-validator.md`. Threading model and per-thread
allocation patterns common to all layers live in `AGENTS.md`. This
directory is for the unique reasoning each layer encodes.

## When to update

Update the relevant rule file when the layer's design changes in a
way that's not captured by code style. New algorithms, new threading
contracts, new invariants between data structures. Do NOT mirror
every refactor — only changes that affect how a future contributor
should reason about the layer.

The `/audit-bundle` command's Pairing 10 ("Path-scoped rules layer
alignment") verifies that each rule file's `paths:` glob still matches
an existing layer directory and that the layer chain in
`validate_layer_direction.py` still includes every layer named here.
Pairing 11 ("Rule source-claim verification") additionally checks that
the named source-code constants and class names each rule cites still
exist in source.
