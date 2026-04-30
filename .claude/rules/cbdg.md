---
description: Lancet2 cbdg/ layer rules — bidirected colored de Bruijn graph, k-mer canonicalization, BCALM2 sign continuity, walk enumeration with arena-allocated walk trees, three-color cycle detection. Load when editing src/lancet/cbdg/.
paths:
  - "src/lancet/cbdg/**"
---

# cbdg/ layer rules

`cbdg/` is the algorithmic core: a bidirected colored de Bruijn graph
(BCDG) for local assembly. The model follows BCALM 2's specification
for bidirected k-mer graphs
(https://github.com/GATB/bcalm/blob/v2.2.3/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md).
This layer is where most subtle bugs hide — get the canonicalization
or sign continuity wrong and downstream variants are silently wrong.

## K-mer canonicalization is lex-min(k-mer, RevComp(k-mer))

`Kmer::Kmer(seq)` chooses the canonical orientation by comparing `seq`
against `RevComp(seq)` and picking the lexicographically smaller. The
implementation in `kmer.cpp::IsCanonicallyPlus` walks inward from
both ends, complementing on the fly — `⌈n/2⌉` iterations, no full
RevComp allocation. **Palindromic k-mers** (where `seq == RevComp(seq)`)
are tagged as `PLUS` by convention.

`Kmer::Identifier()` is `HashStr64(canonical_seq)` — the hash is over
the canonical form, so the same biological k-mer produces the same
identifier regardless of which strand a read came from. Code that
hashes `seq` instead of canonical form (e.g., for a side-table) breaks
this contract and the side-table goes out of sync with the graph.

## Sign and Ordering are NOT the same thing

Two distinct concepts that both use {PLUS, MINUS}:

- **`Kmer::Sign`** (`PLUS`/`MINUS`): does the canonical sequence equal
  the read's original orientation (`PLUS`) or its reverse complement
  (`MINUS`)?
- **`Kmer::Ordering`** (`DEFAULT`/`OPPOSITE`): a query parameter for
  `SignFor`/`SequenceFor` — do you want the canonical orientation or
  its reverse?

`SignFor(DEFAULT)` returns the stored sign; `SignFor(OPPOSITE)`
returns the flipped sign. `SequenceFor(DEFAULT)` returns the canonical
sequence (zero-alloc via `SeqView()`); `SequenceFor(OPPOSITE)`
allocates a fresh `RevComp` string. Use `SeqView()` on hot paths;
allocating `SequenceFor(DEFAULT)` per call is a measurable regression
on chromosome-scale runs.

## EdgeKind encodes both endpoints' signs

`EdgeKind` is one of {`PLUS_PLUS`, `PLUS_MINUS`, `MINUS_PLUS`,
`MINUS_MINUS`} — the source and destination signs together. The
`MergeCords` switch in `kmer.cpp` documents each case with ASCII art
showing how the merged sequence is constructed (overlap region
discarded, non-overlapping suffix or RevComp(suffix) appended). Don't
modify `MergeCords` without working through all four cases on paper —
a bug here produces silently corrupted assemblies.

`MakeFwdEdgeKind`/`SplitIntoSignPair`/`RevEdgeKind` are the canonical
helpers; use them instead of switching on `EdgeKind` directly.
`RevEdgeKind` flips both ends: `PLUS_PLUS ↔ MINUS_MINUS`,
`PLUS_MINUS ↔ PLUS_MINUS` (self-symmetric), `MINUS_PLUS ↔ MINUS_PLUS`
(self-symmetric).

## Three tracking systems coexist on each Node

`Node` carries three independently-updated state buckets:

1. **`mLabel`** (3-bit bitmask, CTRL|CASE|REFERENCE) — which roles
   contributed reads to this k-mer. Used for graph pruning, DOT
   visualization, somatic state classification.
2. **`mCounts`** (`InlinedVector<u32, 2>`) — per-sample read counts
   indexed by sample index. Lazy-grown by
   `IncrementReadSupport(sample_index, ...)` — **the vector grows
   only as far as the highest sample index seen**, so its size is NOT
   the authoritative sample count.
3. **`mRoleCounts`** (2-element `array<u32, 2>`) — per-role
   aggregates (sum across all CTRL samples, sum across all CASE
   samples). `RoleIndex(CTRL) = 0`, `RoleIndex(CASE) = 1`. Fast-path
   for pruning. **`REFERENCE` is never tracked here** — only
   case/control coverage matters for read support.

These three must be updated together. `IncrementReadSupport` is the
only correct entry point; touching one without the others
desynchronizes them. The `mCounts.size()` lazy-growth is the
single most likely place to introduce a subtle bug — see the
Confidence rule below.

## Node::Confidence requires authoritative num_samples

`Node::Confidence(num_samples)` computes
`floor(TotalReadSupport · concordance) + (REFERENCE ? 1 : 0)` where
`concordance = confirming_samples / num_samples`. **The `num_samples`
parameter MUST be the authoritative total from `GraphParams`, NEVER
derived from `mCounts.size()`.** The vector is lazy-grown; nodes that
only saw reads from low-index samples have an undersized `mCounts`,
which would inflate concordance to 1.0 and silently overcount weakly-
supported variants.

The `+1` for reference k-mers is additive (not multiplicative) so it
doesn't distort the coverage scale — REF and ALT path weights remain
directly comparable in the MaxFlow walk scoring.

## MateMer dedup: at most one increment per (qname, kmer, role)

`Graph::BuildGraph` uses `MateMer{qname, kmer_hash, tag_kind}` as the
dedup key. A single read is allowed to contribute at most one
increment per k-mer per sample — without this, paired-end overlaps
double-count k-mers, and shotgun-style amplicon data triple- or
quadruple-counts. The `MateMer` hash combines all three fields; don't
remove any of them when refactoring.

## Walk enumeration: arena-allocated tree, not vector copies

`MaxFlow::NextPath` performs BFS over the bidirected graph using a
flat arena of `WalkTreeNode {edge_ordinal, dst_state, parent_idx,
score}` (16 bytes per arena node). Walks are reconstructed by
parent-pointer traversal only when a walk reaches the sink with
`score > 0`. The old design copied entire walk vectors at each BFS
extension; that produced exponential memory blowup on branchy graphs.

`DEFAULT_GRAPH_TRAVERSAL_LIMIT = 2^20` caps BFS visits per call. A
graph that hits this limit terminates with `mHitTraversalLimit = true`
rather than running forever — the caller can retry with a larger k.
Don't raise the limit casually; pathological graphs (centromeric
repeats) genuinely have >2^20 walks.

The "max-flow" name refers to scoring walks by un-traversed edges:
each `NextPath` call returns a walk that uses at least one new edge,
preferring walks with the most new edges. After enough calls, every
edge is covered and `NextPath` returns `nullopt`. This produces a
deduplicated set of haplotypes ranked by minimum-weight path
confidence.

## Cycle detection: three-color DFS, not backtracking

`cycle_finder.cpp::HasCycle` is three-color DFS (WHITE/GRAY/BLACK).
The old backtracking version was exponential and showed ~51.6s in
profiles; this version is O(V+E) and runs in <1ms. **DFS state is
keyed by `(node_flat_idx, sign)` not just `node_flat_idx`** — in the
bidirected model, a node visited via `+` and via `−` are different
states. The `TraversalIndex` adjacency list is partitioned by `(node,
sign)` precisely so edge iteration respects sign continuity for free.

## Sign continuity: edge.DstSign must equal next edge.SrcSign

This is the BCALM 2 invariant for valid walks in a bidirected graph.
Walks that violate sign continuity are not biologically meaningful
and produce garbage haplotypes. `TraversalIndex`'s flat adjacency
list is structured so iteration naturally enforces continuity — code
that builds a custom adjacency must preserve this property explicitly.

## GraphComplexity gating: AND, not OR

`GraphComplexity::IsComplex()` returns `mCyclomaticComplexity ≥ 50 &&
mNumBranchPoints ≥ 50` (BOTH conditions). The thresholds are derived
from whole-chromosome chr4 profiling (235K component entries):
CC≥50 AND BP≥50 catches 146 pathological windows averaging 5.8s vs.
233K normal windows averaging 414ms — a 14× speedup gate. AND (not
OR) avoids skipping windows where only one metric is elevated; both
have to agree the graph is pathological before triggering k-mer
upsizing.

When `IsComplex()` returns true, the caller retries assembly with a
larger k. Larger k collapses short repeats and merges paralogous
paths, simplifying the graph toward its true biological structure.

## GEI is the ML-feature encoding of complexity

`GraphEntanglementIndex()` returns `log10(1 + (CC × BP × CovCV) /
(UnitigRatio + ε))`. The multiplication is a soft-AND gate: all
three numerator terms must be elevated for the product to spike. The
log10 squash bounds the output to roughly [0, 7]. Don't expose the
raw numerator as a feature directly — the bounded log scale is what
makes it tractable for downstream ML models.
