---
description: Lancet2 hts/ layer rules — HTSlib RAII wrapping, the zero-copy bam1_t lifetime contract, BGZF stream lifecycle, per-thread Extractor invariant. Load when editing src/lancet/hts/.
paths:
  - "src/lancet/hts/**"
---

# hts/ layer rules

`hts/` is the RAII layer over HTSlib's C resources. The library is
not thread-safe at the per-handle level, allocates aggressively in
the C heap, and uses ownership conventions that don't match modern
C++. This layer's job is to translate that into safe, idiomatic C++
without paying allocation cost on the hot path.

## extern "C" wrapping is mandatory at every htslib include

Every `#include "htslib/<x>.h"` is wrapped in `extern "C" { ... }`.
htslib headers are C, and including them without `extern "C"` results
in linker failures on some toolchains (and silent C++ name mangling
elsewhere). Pattern-match the existing wrappers in
`alignment.h`/`extractor.h`/`reference.h`/`bgzf_ostream.h` when adding
a new htslib include.

## RAII via custom deleters — never raw destroy calls

Every htslib C resource has a paired deleter struct in
`extractor.h::detail`: `HtsFileDeleter` (`hts_close`), `SamHdrDeleter`
(`sam_hdr_destroy`), `HtsIdxDeleter` (`hts_idx_destroy`),
`HtsItrDeleter` (`hts_itr_destroy`), `Bam1Deleter` (`bam_destroy1`),
`HtsFilterDeleter` (`hts_filter_free`). Use these via
`std::unique_ptr<T, Deleter>`; never call `*_destroy` or `*_free`
directly. New htslib resources added to this layer follow the same
pattern: deleter struct → `unique_ptr` alias.

## The bam1_t lifetime contract: Alignment is non-owning, invalidated on ++itr

`Alignment` is a zero-copy proxy over the underlying `bam1_t*` that
is structurally managed by `hts::Iterator`. The contract is **the
single most important invariant in this layer**:

1. The pointer becomes invalid after the next `++itr`.
2. `string_view`s returned by `QnameView()`, `CigarData()`, etc. are
   tied to the same lifetime.
3. Data that must outlive the iterator step requires
   `BuildSequence()` / `BuildQualities()`, which deep-copy.
4. **Never** store an `Alignment` object beyond the loop body;
   **never** capture its `string_view` returns past `++itr`.

Violations are silent: the buffer gets reused for the next record, so
the stale `string_view` reads garbage that looks like valid data.
Every existing call site decomposes the `Alignment` into owned types
(typically `cbdg::Read`) inside the loop body before iterating.
Pattern-match this when adding new consumers.

## Per-thread Extractor / Iterator instances are required

HTSlib iterators are not thread-safe — sharing a single
`hts::Extractor` across threads produces silent data corruption (the
read offset is shared mutable state inside htslib). Each worker
thread owns its own `Extractor`/`Iterator`. The per-thread allocation
pattern is encoded in `core/read_collector.h::SampleExtractors`
(per-sample extractor map, one map per thread). New code that touches
BAM/CRAM I/O across multiple threads must construct per-thread
`Extractor` instances; sharing one is a regression even if it
"appears to work."

## Phred quality conversion is constexpr-table-backed

`phred_quality.h::PhredToErrorProb` uses a precomputed constexpr
lookup table sized by `MAX_PHRED_SCORE = 255`. Zero runtime
initialization, thread-safe by construction, branch-free. Do not
replace with a runtime `std::pow(10, -Q/10)` call — every variant
window touches this path. The reverse direction
(`ErrorProbToPhred`) uses `std::log10` since the input is continuous;
that's intentional, not a TODO.

The Phred encoding rule (`Q = −10·log₁₀(P_error)`) is the convention
the entire caller layer assumes. Code reading Phred scores from BAMs
must apply the htslib offset (`qual + 33` is the ASCII-encoded form
in the BAM file; the raw 0..93 range is what `BuildQualities()`
returns).

## BgzfOstream: per-thread, with explicit lifecycle

`BgzfOstream` is **not thread-safe** — concurrent writes to the same
underlying BGZF handle corrupt the block structure. There is exactly
one `BgzfOstream` for the output VCF, owned by `PipelineRunner` on
the main thread; workers send variants to `VariantStore` and only
the main thread flushes through this stream. The lifecycle is:

1. `Open(path, BgzfFormat::VCF)` — opens htslib's BGZF
2. write the header text + `flush()` — header committed before any
   worker starts (this is what makes the output incrementally
   tail-readable)
3. workers add to `VariantStore`; main thread invokes
   `FlushVariantsBeforeWindow` which writes BGZF-blocked records and
   calls `flush()`
4. `Close()` — finalizes the BGZF stream (writes the EOF marker
   block; without this, the file is unreadable by `bcftools`)

A path that returns from `Run()` without calling `Close()` produces a
silently-truncated BGZF file. The pipeline runner uses
`std::exit(EXIT_SUCCESS)` only after `output_vcf.Close()`; preserve
this ordering.

## Cloud URIs validate credentials upfront

`hts::IsCloudUri` + `hts::ValidateCloudAccess` exist because libcurl's
multipart upload threshold is 5MB; small VCF headers won't trigger
the HTTP handshake until `BgzfOstream::Close()` flushes after the
full pipeline run. A 40-hour pipeline that fails authentication at
the end is the worst possible failure mode. The pre-flight
zero-byte PUT validates GCP/S3 credentials before any work starts.
New cloud-output code paths must call `ValidateCloudAccess` before
the long-running work, not after.

## Reference cache size is intentionally large

`reference.cpp` sets `DEFAULT_FAI_CACHE_SIZE = 1 << 28` (256 MiB) so
a full chromosome's blocks load into memory. This matters for
multi-threaded random access — without the cache, threads serialize
on the FAI seek. Don't shrink this without profile data.
