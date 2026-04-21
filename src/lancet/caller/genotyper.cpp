#include "lancet/caller/genotyper.h"

#include "lancet/base/types.h"
#include "lancet/caller/allele_scoring_types.h"
#include "lancet/caller/combined_scorer.h"
#include "lancet/caller/local_scorer.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/variant_set.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/read.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/cigar_utils.h"

// minimap2 C API — POSIX/C headers, not C++ stdlib
extern "C" {
#include "minimap.h"
#include "mmpriv.h"
}

#include "absl/hash/hash.h"
#include "absl/types/span.h"

#include <limits>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

#include <cstdlib>

namespace {

// Free minimap2 alignment results (mm_reg1_t array)
// NOLINTBEGIN(cppcoreguidelines-no-malloc)
inline void FreeMm2Alignment(mm_reg1_t* regs, int const num_regs) {
  if (regs == nullptr) return;
  for (int idx = 0; idx < num_regs; ++idx) std::free(regs[idx].p);
  std::free(regs);
}
// NOLINTEND(cppcoreguidelines-no-malloc)

// Build a CIGAR vector from a minimap2 alignment result, inserting leading and
// trailing soft-clip (S) operations from qs/qe.  Even with end_bonus=INT_MAX,
// ~0.1% of reads retain small clips.  Without explicit S operations the downstream
// ComputeSoftClipPenalty() returns 0, failing to crush chimeric partial alignments.
auto BuildCigar(mm_reg1_t const* hit, int qlen) -> std::vector<lancet::hts::CigarUnit> {
  std::vector<lancet::hts::CigarUnit> cigar;
  if (hit->p == nullptr || hit->p->n_cigar == 0) return cigar;

  auto const nops = static_cast<int>(hit->p->n_cigar);
  // Reserve: core ops + up to 2 soft-clip bookends
  cigar.reserve(static_cast<usize>(nops) + 2);

  // Leading soft-clip: query bases [0, qs) were not aligned
  if (hit->qs > 0) {
    cigar.emplace_back(lancet::hts::CigarOp::SOFT_CLIP, static_cast<u32>(hit->qs));
  }

  // Core CIGAR operations from minimap2
  for (int i = 0; i < nops; ++i) {
    cigar.emplace_back(hit->p->cigar[i]);
  }

  // Trailing soft-clip: query bases [qe, qlen) were not aligned
  if (hit->qe < qlen) {
    cigar.emplace_back(lancet::hts::CigarOp::SOFT_CLIP, static_cast<u32>(qlen - hit->qe));
  }

  return cigar;
}

}  // namespace

namespace lancet::caller {

// ============================================================================
// Constructor: initialize minimap2 for read-to-haplotype allele classification.
//
// Parameter strategy: start from defaults, override 10 fields. We do NOT use
// the 'sr' preset wholesale — its gap/scoring parameters target whole-genome
// mapping where gaps are biology.  Here, biology is baked into the haplotype
// sequences; read-haplotype divergence is sequencer noise and must be penalized
// strictly.  However, we DO set MM_F_SR to activate the SR extension-region
// code path (full-query boundaries + end_bonus reference expansion), which is
// the correct model for 151 bp reads on short synthetic haplotypes.
//
// See scoring_constants.h for the SCORING_* values and
// docs/guides/variant_discovery_genotyping.md for the design rationale.
// ============================================================================
Genotyper::Genotyper() {
  // 0 -> no info, 1 -> error, 2 -> warning, 3 -> debug
  mm_verbose = 1;

  // Start from the default parameter set, then override
  mm_set_opt(nullptr, mIndexingOpts.get(), mMappingOpts.get());

  auto* mopts = mMappingOpts.get();

  // ============================================================================
  // Flags: CIGAR generation + SR extension boundaries
  // ============================================================================
  // MM_F_CIGAR: required for local scoring.
  // MM_F_SR: activates full-query extension (qs0=0, qe0=qlen) and uses
  //   end_bonus to expand the reference extraction region at read edges.
  //   Also disables irrelevant long-read paths (inversion detection,
  //   bw_long re-chain, strand error estimation).  All 10 SR code paths
  //   are verified safe for single-segment haplotype alignment.
  // NOTE: flag |= only sets bit flags — it does not invoke mm_set_opt("sr")
  //   and therefore does not overwrite any parameters set below.
  mopts->flag |= MM_F_CIGAR | MM_F_SR;
  mopts->best_n = 1;  // single best hit per haplotype

  // ============================================================================
  // Scoring: single-affine, strict penalties
  // ============================================================================
  // Gaps/mismatches here are noise, not biology.  Setting q2=q, e2=e
  // disables minimap2's dual-affine (convex) model — no "cheap extension"
  // loophole for long gaps that might bleed reads to the wrong haplotype.
  mopts->a = SCORING_MATCH;
  mopts->b = SCORING_MISMATCH;
  mopts->q = SCORING_GAP_OPEN;
  mopts->e = SCORING_GAP_EXTEND;
  mopts->q2 = SCORING_GAP_OPEN;
  mopts->e2 = SCORING_GAP_EXTEND;

  // ============================================================================
  // Z-Drop: effectively disabled
  // ============================================================================
  // A 300 bp somatic deletion incurs gap penalty O + 300·E = 12 + 900 = 912,
  // exceeding the default zdrop=400.  Setting to 100k prevents DP truncation
  // across any intra-window structural variation up to ~5 kbp.
  mopts->zdrop = 100'000;
  mopts->zdrop_inv = 100'000;

  // ============================================================================
  // Bandwidth: envelope insertions up to ~2 kbp
  // ============================================================================
  // The DP band width constrains how far off-diagonal the SW matrix extends.
  // A 50 bp insertion requires 50 off-diagonal cells.  bw=10k guarantees
  // large germline and somatic insertions never exceed the band.
  mopts->bw = 10'000;

  // ============================================================================
  // Chaining gap: cap search radius to read length
  // ============================================================================
  // max_gap controls the backward search in mg_lchain_dp.  Haplotypes are
  // <1 kbp and reads are 151 bp — no valid chain spans a 200 bp gap.
  // Reducing from the default 5000 eliminates wasted inner-loop iterations.
  mopts->max_gap = 200;

  // ============================================================================
  // End bonus: force extension through edge variants
  // ============================================================================
  // Without this, KSW2's EXTZ_ONLY mode backtracks to the max-score cell
  // and soft-clips the remaining query — hiding variants near read edges.
  // end_bonus only redirects the backtrack endpoint; it does not inflate
  // alignment scores (dp_score always uses ez->max, not mqe + end_bonus)
  // and has zero performance cost.  INT_MAX guarantees the backtrack always
  // reaches the query end, eliminating most of erroneous soft-clips.
  mopts->end_bonus = std::numeric_limits<int>::max();

  // ============================================================================
  // Seeding: dense anchors for short, mutated haplotypes
  // ============================================================================
  // k=11, w=5 places minimizer seeds in densely mutated micro-windows
  // where the default k=15 would fail to produce a continuous exact match.
  // Per-haplotype indexing makes k=11's higher false-positive rate harmless.
  mIndexingOpts->k = 11;
  mIndexingOpts->w = 5;
}

// ============================================================================
// Genotype: main entry point — orchestrates alignment and evidence collection.
//
//  ┌──────────────┐    ┌──────────────┐    ┌─────────────┐
//  │  Haplotypes  │    │    Reads     │    │  VariantSet │
//  │ (REF + ALTs) │    │  (all samps) │    │ (raw vars)  │
//  └───────┬──────┘    └──────┬───────┘    └─────┬───────┘
//          │                  │                  │
//          ▼                  │                  │
//   ResetData()               │                  │
//   (build mm2 indices)       │                  │
//          │                  │                  │
//          │    ┌─────────────┘                  │
//          │    │                                │
//          ▼    ▼                                │
//   AssignReadToAlleles() ◄──────────────────────┘
//   (mm_map per hap,
//    local scoring per var)
//          │
//          ▼
//   AddToTable()
//   (ReadEvidence → VariantSupport)
//          │
//          ▼
//   ┌──────────────┐
//   │    Result    │
//   │  per-variant │
//   │  per-sample  │
//   │  support     │
//   └──────────────┘
// ============================================================================
auto Genotyper::Genotype(Haplotypes hap_seqs, Reads qry_reads, VariantSet const& variant_set)
    -> Result {
  ResetData(hap_seqs);
  Result out_vars_table;

  for (auto const& qry_read : qry_reads) {
    auto allele_assignments = AssignReadToAlleles(qry_read, variant_set);
    AddToTable(out_vars_table, qry_read, allele_assignments);
  }

  return out_vars_table;
}

// ============================================================================
// ResetData: build minimap2 indices for all haplotype sequences.
//
// Each haplotype gets its own index so we can align reads independently
// to REF and each ALT haplotype and compare scores.
// ============================================================================
void Genotyper::ResetData(Haplotypes hap_seqs) {
  mIndices.clear();
  mIndices.reserve(hap_seqs.size());

  auto const* iopts = mIndexingOpts.get();
  for (auto const& hap_seq : hap_seqs) {
    char const* raw_seq = hap_seq.c_str();
    auto* idx_result = mm_idx_str(iopts->w, iopts->k, 0, iopts->bucket_bits, 1, &raw_seq, nullptr);
    mIndices.emplace_back(Minimap2Index(idx_result));
  }

  // Pre-encode haplotype sequences for local scoring.
  // mm_idx stores sequences internally but doesn't expose them via a clean API,
  // so we maintain our own numeric-encoded copies for ComputeLocalScore.
  mEncodedHaplotypes.clear();
  mEncodedHaplotypes.reserve(hap_seqs.size());
  for (auto const& hap_seq : hap_seqs) {
    mEncodedHaplotypes.push_back(EncodeSequence(hap_seq));
  }

  auto* mopts = mMappingOpts.get();
  for (auto const& mm2_idx : mIndices) {
    mm_mapopt_update(mopts, mm2_idx.get());
  }
}

auto Genotyper::AssignReadToAlleles(cbdg::Read const& qry_read, VariantSet const& variant_set)
    -> PerVariantAssignment {
  auto all_alns = AlignToAllHaplotypes(qry_read);
  if (all_alns.empty()) return {};

  auto const qry_quals = qry_read.QualView();
  auto const qry_seq_encoded = EncodeSequence(qry_read.SeqView());
  usize const qry_read_length = qry_read.Length();

  ReadAlnContext const read_ctx{
      .mSeqEncoded = absl::MakeConstSpan(qry_seq_encoded),
      .mBaseQuals = qry_quals,
      .mReadLength = qry_read_length,
  };

  // O(N) PERFORMANCE WIN: Extracted from the variant iterator loop.
  // Calculated exactly once per read.
  u32 const baseline_ref_nm = ComputeHaplotypeEditDistance(
      all_alns, absl::MakeConstSpan(mEncodedHaplotypes[REF_HAP_IDX]),
      absl::MakeConstSpan(qry_seq_encoded), qry_read_length, REF_HAP_IDX);

  PerVariantAssignment allele_assignments;

  // For each haplotype alignment, score all overlapping variants.
  // all_alns is tiny (~2–10 items, one per assembled haplotype).
  for (auto const& aln : all_alns) {
    auto const haplotype = absl::MakeConstSpan(mEncodedHaplotypes[aln.mHapIdx]);

    for (auto const& variant : variant_set) {
      auto const bounds = ExtractHapBounds(variant, aln.mHapIdx);
      if (!bounds || !OverlapsAlignment(aln, *bounds)) continue;

      auto scored = ScoreReadAtVariant(aln, haplotype, read_ctx, *bounds);
      scored.mRefNm = baseline_ref_nm;

      // Reuse the hash probe: find once, then compare, update-in-place
      // or emplace_hint — avoids searching the map twice.
      auto iter = allele_assignments.find(&variant);
      if (iter != allele_assignments.end() &&
          scored.CombinedScore() <= iter->second.CombinedScore()) {
        continue;
      }

      if (iter != allele_assignments.end()) {
        iter->second = scored;
      } else {
        allele_assignments.emplace_hint(iter, &variant, scored);
      }
    }
  }

  return allele_assignments;
}

// ============================================================================
// ExtractHapBounds: resolve a minimap2 alignment's haplotype index to the
// variant's physical coordinates on that specific haplotype.
//
// See genotyper.h for full documentation and ASCII diagram.
// ============================================================================
auto Genotyper::ExtractHapBounds(RawVariant const& variant, usize aln_hap_idx)
    -> std::optional<HapVariantBounds> {
  if (aln_hap_idx == REF_HAP_IDX) {
    return HapVariantBounds{
        .mVarStart = static_cast<i32>(variant.mLocalRefStart0Idx),
        .mVarLen = static_cast<i32>(variant.mRefAllele.size()),
        .mAllele = REF_ALLELE_IDX,
    };
  }
  for (usize alt_pos = 0; alt_pos < variant.mAlts.size(); ++alt_pos) {
    auto const& alt_allele = variant.mAlts[alt_pos];
    auto iter = alt_allele.mLocalHapStart0Idxs.find(aln_hap_idx);
    if (iter != alt_allele.mLocalHapStart0Idxs.end()) {
      return HapVariantBounds{
          .mVarStart = static_cast<i32>(iter->second),
          // Exact length of the ALT allele sequence because
          // multiallelics can have different lengths than the REF allele.
          .mVarLen = static_cast<i32>(alt_allele.mSequence.size()),
          .mAllele = static_cast<AlleleIndex>(alt_pos + 1),
      };
    }
  }
  return std::nullopt;
}

// ============================================================================
// OverlapsAlignment: check if a minimap2 alignment overlaps a variant region.
//
// See genotyper.h for full documentation and ASCII diagram.
// ============================================================================
auto Genotyper::OverlapsAlignment(Mm2AlnResult const& aln, HapVariantBounds const& bounds) -> bool {
  i32 const hap_var_end = bounds.mVarStart + bounds.mVarLen;
  return hap_var_end > aln.mRefStart && bounds.mVarStart < aln.mRefEnd;
}

// ============================================================================
// AlignToAllHaplotypes: map a read against all haplotype sequences.
//
// Uses minimap2's full pipeline (seed → chain → extend) rather than raw ksw2 DP
// because minimizer hashing resolves topological boundaries before falling back
// to dynamic programming.
//
// CRITICAL: No early-exit short-circuit — even if a read perfectly matches one
// haplotype, we must align to all haplotypes to compute correct cross-haplotype
// noise constraints, relative edit distances against the reference baseline,
// and global best-match boundaries.
// ============================================================================
auto Genotyper::AlignToAllHaplotypes(cbdg::Read const& qry_read) -> std::vector<Mm2AlnResult> {
  std::vector<Mm2AlnResult> results;
  results.reserve(mIndices.size());

  int nregs = 0;
  auto* tbuffer = mThreadBuffer.get();
  auto const* map_opts = mMappingOpts.get();
  auto const read_len = static_cast<int>(qry_read.Length());

  for (usize idx = 0; idx < mIndices.size(); ++idx) {
    auto const* hap_mm_idx = mIndices[idx].get();
    auto* regs = mm_map(hap_mm_idx, read_len, qry_read.SeqPtr(), &nregs, tbuffer, map_opts,
                        qry_read.QnamePtr());

    if (regs == nullptr || nregs <= 0) {
      FreeMm2Alignment(regs, nregs);
      continue;
    }

    // Take the top hit only (best_n = 1)
    mm_reg1_t const* top_hit = &regs[0];

    Mm2AlnResult result;
    result.mScore = top_hit->score;
    result.mRefStart = top_hit->rs;  // critical: where alignment starts on haplotype
    result.mRefEnd = top_hit->re;
    result.mIdentity = mm_event_identity(top_hit);
    result.mHapIdx = idx;
    result.mCigar = BuildCigar(top_hit, read_len);

    results.push_back(std::move(result));
    FreeMm2Alignment(regs, nregs);
  }

  return results;
}

// ============================================================================
// AddToTable: record a read's allele assignments into per-variant support.
//
// Constructs ReadEvidence with all available read-level metrics:
//   - base quality at the variant position
//   - original mapping quality
//   - normalized alignment score
//   - strand direction
// All metrics flow through to VCF FORMAT fields via VariantSupport aggregation.
// ============================================================================
void Genotyper::AddToTable(Result& out_vars_table, cbdg::Read const& qry_read,
                           PerVariantAssignment const& allele_assignments) {
  auto const sample_name = qry_read.SampleName();
  auto const rname_hash = absl::HashOf(qry_read.QnameView());
  auto const strand = qry_read.Flag().IsRevStrand() ? Strand::REV : Strand::FWD;

  for (auto const& [var_ptr, assignment] : allele_assignments) {
    // Look up (or create) the per-sample evidence aggregator for this variant.
    // Result is keyed: variant → sample_name → VariantSupport.
    // Default insertion occurs at both tier levels:
    //   rslt[var_ptr]               → inserts an empty SupportArray interface
    //   .FindOrCreate(sample_name)  → creates a unique_ptr<VariantSupport> on first access.
    auto& support = out_vars_table[var_ptr].FindOrCreate(sample_name);

    auto const evidence = VariantSupport::ReadEvidence{
        .mInsertSize = qry_read.InsertSize(),
        .mAlignmentStart = qry_read.StartPos0(),
        .mAlnScore = assignment.CombinedScore(),
        .mFoldedReadPos = assignment.mFoldedReadPos,
        .mRnameHash = static_cast<u32>(rname_hash),
        .mRefNm = assignment.mRefNm,
        .mOwnHapNm = assignment.mOwnHapNm,
        .mAssignedHaplotypeId = assignment.mAssignedHaplotypeId,
        .mAllele = assignment.mAllele,
        .mStrand = strand,
        .mBaseQual = assignment.mBaseQualAtVar,
        .mMapQual = qry_read.MapQual(),
        .mIsSoftClipped = qry_read.IsSoftClipped(),
        .mIsProperPair = qry_read.IsProperPair(),
    };

    support.AddEvidence(evidence);
  }
}

}  // namespace lancet::caller
