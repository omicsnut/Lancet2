#ifndef SRC_LANCET_HTS_CIGAR_UNIT_H_
#define SRC_LANCET_HTS_CIGAR_UNIT_H_

#include "lancet/base/types.h"

extern "C" {
#include "htslib/sam.h"
}

#include <array>

namespace lancet::hts {

enum class CigarOp : char {
  /// Bases aligned to reference without evidence for indel.
  /// No indication whether the read bases match the reference.
  /// Consumes both query and reference sequences.
  ALIGNMENT_MATCH = 'M',

  /// Bases from the read inserted into the reference.
  /// Consumes only query sequence.
  INSERTION = 'I',

  /// Bases from the reference deleted in the read.
  /// Consumes only reference sequence.
  DELETION = 'D',

  /// Bases from the read have skipped the reference, but have not been deleted.
  /// Consumes only reference sequence.
  REFERENCE_SKIP = 'N',

  /// Bases from the read omitted from alignment, but left in the read.
  /// Consumes only query sequence.
  SOFT_CLIP = 'S',

  /// Bases from the read omitted from alignment and removed from the read.
  /// Consumes neither query nor reference sequence.
  HARD_CLIP = 'H',

  /// Used to represent a padding in both query and reference.
  /// Consumes neither query nor reference sequence.
  ALIGNMENT_PAD = 'P',

  /// Bases aligned and exactly matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MATCH = '=',

  /// Bases aligned but not matching to reference.
  /// Consumes both query and reference sequences.
  SEQUENCE_MISMATCH = 'X',

  /// only present to handle all other cases when alignment mCigar is corrupt
  UNKNOWN_OP = '?'
};

// ── CIGAR consume-query/reference lookup tables ─────────────────────────────
// Branchless O(1) lookups indexed by the underlying CigarOp char value.
// SAM spec §1.4.6: M/I/S/=/X consume query; M/D/N/=/X consume reference.
//
// 128-entry LUT avoids needing to remap char → dense index.
// Only entries at 'M','I','D','N','S','H','P','=','X' positions are set;
// everything else defaults to false (including UNKNOWN_OP '?').
namespace detail {

constexpr auto MakeConsumesRefLut() -> std::array<bool, 128> {
  std::array<bool, 128> lut{};
  lut['M'] = true;  // ALIGNMENT_MATCH
  lut['D'] = true;  // DELETION
  lut['N'] = true;  // REFERENCE_SKIP
  lut['='] = true;  // SEQUENCE_MATCH
  lut['X'] = true;  // SEQUENCE_MISMATCH
  return lut;
}

constexpr auto MakeConsumesQryLut() -> std::array<bool, 128> {
  std::array<bool, 128> lut{};
  lut['M'] = true;  // ALIGNMENT_MATCH
  lut['I'] = true;  // INSERTION
  lut['S'] = true;  // SOFT_CLIP
  lut['='] = true;  // SEQUENCE_MATCH
  lut['X'] = true;  // SEQUENCE_MISMATCH
  return lut;
}

inline constexpr auto CONSUMES_REF = MakeConsumesRefLut();
inline constexpr auto CONSUMES_QRY = MakeConsumesQryLut();

}  // namespace detail

class CigarUnit {
 public:
  // Implicit conversion from BAM u32 cigar encoding enables span-based assign().
  // NOLINTNEXTLINE(google-explicit-constructor)
  CigarUnit(u32 sam_cigop)
      : mLength(bam_cigar_oplen(sam_cigop)),
        mCigOp(static_cast<CigarOp>(bam_cigar_opchr(sam_cigop))) {}

  CigarUnit() = delete;

  [[nodiscard]] auto Operation() const noexcept -> CigarOp { return mCigOp; }
  [[nodiscard]] auto Length() const noexcept -> u32 { return mLength; }

  [[nodiscard]] auto ConsumesReference() const noexcept -> bool {
    auto const idx = static_cast<unsigned char>(mCigOp);
    return idx < detail::CONSUMES_REF.size() && detail::CONSUMES_REF[idx];
  }

  [[nodiscard]] auto ConsumesQuery() const noexcept -> bool {
    auto const idx = static_cast<unsigned char>(mCigOp);
    return idx < detail::CONSUMES_QRY.size() && detail::CONSUMES_QRY[idx];
  }

 private:
  // ── 4B Align ────────────────────────────────────────────────────────────
  u32 mLength;
  // ── 1B Align ────────────────────────────────────────────────────────────
  CigarOp mCigOp;
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_CIGAR_UNIT_H_
