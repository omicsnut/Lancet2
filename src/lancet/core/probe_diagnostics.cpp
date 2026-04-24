#include "lancet/core/probe_diagnostics.h"

#include "lancet/base/types.h"
#include "lancet/caller/alt_allele.h"
#include "lancet/caller/per_allele_data.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/caller/support_array.h"
#include "lancet/caller/variant_support.h"
#include "lancet/cbdg/probe_index.h"
#include "lancet/cbdg/probe_tracker.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"

#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lancet::core {
namespace {

// ============================================================================
// VariantMatch: result of searching the MSA VariantSet for a truth probe.
//
// Stores a pointer to the matching RawVariant (null if no match) and the
// 1-based AlleleIndex of the truth ALT allele within that variant.
// ============================================================================
struct VariantMatch {
  // ── 8B Align ────────────────────────────────────────────────────────────
  caller::RawVariant const* mVariant = nullptr;  // 8B

  // ── 4B Align ────────────────────────────────────────────────────────────
  caller::AlleleIndex mTruthAltIdx = 0;  // 4B
};

// ============================================================================
// FindMatchingVariant: search the VariantSet for a probe's exact position
// and allele match. Returns VariantMatch with mVariant==nullptr on miss.
//
// Matching criteria (all three must hold):
//   1. RawVariant.mGenomeChromPos1 == probe.mGenomeStart0 + 1
//   2. RawVariant.mRefAllele == probe.mRef
//   3. Some ALT allele has mSequence == probe.mAlt
// ============================================================================
[[nodiscard]] auto FindMatchingVariant(caller::VariantSet const& variant_set,
                                       cbdg::ProbeVariant const& probe) -> VariantMatch {
  auto const probe_pos1 = probe.mGenomeStart0 + 1;
  for (auto const& raw_var : variant_set) {
    if (raw_var.mGenomeChromPos1 != probe_pos1) continue;
    if (raw_var.mRefAllele != probe.mRef) continue;

    for (usize idx = 0; idx < raw_var.mAlts.size(); ++idx) {
      if (raw_var.mAlts[idx].mSequence != probe.mAlt) continue;

      return {.mVariant = &raw_var, .mTruthAltIdx = static_cast<caller::AlleleIndex>(idx + 1)};
    }
  }

  return {};
}

// ============================================================================
// CountStolenReads: count ALT-carrying reads misassigned to REF or wrong ALT.
//
// For each allele in the genotyper's VariantSupport, check whether any read
// that was assigned to that allele also appears in the probe's read ownership
// set. If that read was assigned to REF → mGenoStolenToRef++. If assigned
// to a different ALT → mGenoStolenToWrongAlt++.
// ============================================================================
void CountStolenReads(caller::VariantSupport const& support,
                      absl::flat_hash_set<u32> const& probe_reads,
                      caller::AlleleIndex truth_alt_idx, cbdg::ProbeKRecord& record) {
  auto const allele_data = support.AlleleData();
  for (usize idx = 0; idx < allele_data.size(); ++idx) {
    for (auto const& [qname_hash, strand] : allele_data[idx].mNameHashes) {
      if (!probe_reads.contains(qname_hash)) continue;

      if (idx == caller::REF_ALLELE_IDX) {
        ++record.mGenoStolenToRef;
      } else if (idx != truth_alt_idx) {
        ++record.mGenoStolenToWrongAlt;
      }
    }
  }
}

// ============================================================================
// CountNonOverlappingReads: count ALT-carrying reads that did not overlap
// the genotyper alignment window at this variant.
// ============================================================================
void CountNonOverlappingReads(caller::VariantSupport const& support,
                              absl::flat_hash_set<u32> const& probe_reads,
                              cbdg::ProbeKRecord& record) {
  absl::flat_hash_set<u32> all_genotyped;
  auto const allele_data = support.AlleleData();

  for (auto const& allele : allele_data) {
    for (auto const& [qname_hash, strand] : allele.mNameHashes) {
      all_genotyped.insert(qname_hash);
    }
  }

  for (auto const qname_hash : probe_reads) {
    if (!all_genotyped.contains(qname_hash)) ++record.mGenoNonOverlapping;
  }
}

// ============================================================================
// FindShiftedMatch: search the VariantSet for a probe with matching alleles
// but different genomic position (indel left-alignment artifact).
// ============================================================================
void FindShiftedMatch(caller::VariantSet const& variant_set, cbdg::ProbeVariant const& probe,
                      cbdg::ProbeKRecord& record) {
  auto const probe_pos1 = probe.mGenomeStart0 + 1;

  for (auto const& raw_var : variant_set) {
    if (raw_var.mRefAllele != probe.mRef) continue;

    bool allele_matched = false;
    for (auto const& alt : raw_var.mAlts) {
      if (alt.mSequence == probe.mAlt) {
        allele_matched = true;
        break;
      }
    }

    if (!allele_matched) continue;

    auto const shift = static_cast<i64>(raw_var.mGenomeChromPos1) - static_cast<i64>(probe_pos1);
    record.mIsMsaShifted = true;
    record.mMsaShiftBp = static_cast<i16>(shift);
    return;
  }
}

}  // namespace

void ProbeDiagnostics::Initialize(std::filesystem::path const& variants_path,
                                  std::filesystem::path const& results_path,
                                  std::shared_ptr<cbdg::ProbeIndex const> probe_index) {
  if (variants_path.empty()) return;
  auto variants = cbdg::ProbeIndex::LoadVariantsFromFile(variants_path);
  if (variants.empty()) return;
  mTracker.LoadVariants(std::move(variants));
  mTracker.SetResultsPath(results_path);
  mTracker.SetProbeIndex(std::move(probe_index));
}

// ============================================================================
// CheckMsaExtraction: classify how each truth probe appeared (or didn't)
// in the MSA-extracted VariantSet.
//
// Three match tiers (in priority order):
//   1. Exact match: same genomic position + same REF/ALT alleles.
//   2. Shifted match: same REF/ALT alleles but different position (indel
//      left-alignment artifact or flanking-sequence normalization).
//   3. Representation match: variant position matches but truth allele was
//      subsumed into a larger MNV representation.
//
// If no match at any tier, the record's MSA flags all remain false, which
// DeriveLostAt classifies as "msa_no_variant."
// ============================================================================
void ProbeDiagnostics::CheckMsaExtraction(caller::VariantSet const& variant_set,
                                          Window const& window, usize /*component_idx*/) {
  auto const variants = mTracker.Variants();
  if (variants.empty()) return;

  auto const chrom_name = window.ChromName();

  for (auto const probe_id : mTracker.ActiveProbeIds()) {
    auto const& probe = variants[probe_id];
    if (probe.mChrom != chrom_name) continue;

    auto& record = mTracker.FindFinalKRecord(probe_id);
    if (record.mHapIndices.empty()) continue;

    // Tier 1: exact match (position + alleles)
    auto const [matching_var, truth_alt_idx] = FindMatchingVariant(variant_set, probe);
    if (matching_var != nullptr) {
      record.mIsMsaExactMatch = true;
      record.mMsaShiftBp = 0;
      continue;
    }

    // Tier 2: shifted match (same alleles, different position)
    FindShiftedMatch(variant_set, probe, record);
    if (record.mIsMsaShifted) continue;

    // Tier 3: representation match (variant subsumed by MNV)
    CheckMsaRepresentationMatch(variant_set, probe, record);
  }
}

// ============================================================================
// CheckMsaRepresentationMatch: detect when a truth probe's allele is subsumed
// within a larger MNV at the same or overlapping position.
//
// Same-length substitutions only (MNV). The truth variant's position maps to
// a specific offset within the MNV's REF/ALT alleles. Both the REF and ALT
// substrings at that offset must match the truth probe's alleles.
// ============================================================================
void ProbeDiagnostics::CheckMsaRepresentationMatch(caller::VariantSet const& variant_set,
                                                   cbdg::ProbeVariant const& probe,
                                                   cbdg::ProbeKRecord& record) {
  auto const probe_pos1 = probe.mGenomeStart0 + 1;
  for (auto const& raw_var : variant_set) {
    for (auto const& allele : raw_var.mAlts) {
      if (allele.mSequence.size() != raw_var.mRefAllele.size()) continue;

      auto const rv_start = raw_var.mGenomeChromPos1;
      auto const rv_end = rv_start + raw_var.mRefAllele.size();

      if (probe_pos1 < rv_start || probe_pos1 >= rv_end) continue;

      auto const offset = probe_pos1 - rv_start;
      if (offset + probe.mRef.size() > raw_var.mRefAllele.size()) continue;
      if (raw_var.mRefAllele.substr(offset, probe.mRef.size()) != probe.mRef) continue;

      if (offset + probe.mAlt.size() > allele.mSequence.size()) continue;
      if (allele.mSequence.substr(offset, probe.mAlt.size()) != probe.mAlt) continue;

      record.mIsMsaRepresentation = true;
      return;
    }
  }
}

// ============================================================================
// CheckGenotyperResult: verify genotyper read assignment for each probe.
//
// For each probe that achieved an exact MSA match (mIsMsaExactMatch==true),
// finds the corresponding RawVariant in the genotyper result map and:
//   1. Records total ALT and REF coverage from the genotyper
//   2. Counts "stolen reads" — ALT-carrying reads misassigned to REF or
//      a different ALT allele
//   3. Counts non-overlapping reads — ALT-carrying reads that didn't
//      appear in the genotyper alignment window at all
// ============================================================================
void ProbeDiagnostics::CheckGenotyperResult(caller::Genotyper::Result const& genotyped,
                                            caller::VariantSet const& variant_set,
                                            usize /*component_idx*/) {
  auto const variants = mTracker.Variants();
  if (variants.empty()) return;

  for (auto const probe_id : mTracker.ActiveProbeIds()) {
    auto& record = mTracker.FindFinalKRecord(probe_id);
    if (!record.mIsMsaExactMatch) continue;

    auto const [matching_var, truth_alt_idx] = FindMatchingVariant(variant_set, variants[probe_id]);
    if (matching_var == nullptr) continue;

    auto const geno_it = genotyped.find(matching_var);
    if (geno_it == genotyped.end()) {
      record.mIsGenoNoResult = true;
      continue;
    }

    auto const& read_index = mTracker.ReadOwnershipIndex();
    auto const reads_it = read_index.find(probe_id);
    auto const* probe_reads = reads_it != read_index.end() ? &reads_it->second : nullptr;

    for (auto const& named_support : geno_it->second) {
      if (!named_support.mData) continue;

      auto const alt_support = named_support.mData->TotalAlleleCov(truth_alt_idx);
      record.mGenoTrueAltReads += static_cast<u16>(alt_support);
      record.mGenoTotalRefReads += static_cast<u16>(named_support.mData->TotalRefCov());

      if (probe_reads == nullptr) continue;
      CountStolenReads(*named_support.mData, *probe_reads, truth_alt_idx, record);
      CountNonOverlappingReads(*named_support.mData, *probe_reads, record);
    }

    record.mIsGenoHasAltSupport = record.mGenoTrueAltReads > 0;
  }
}

}  // namespace lancet::core
