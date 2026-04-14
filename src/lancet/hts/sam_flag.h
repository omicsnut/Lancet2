#ifndef SRC_LANCET_HTS_SAM_FLAG_H_
#define SRC_LANCET_HTS_SAM_FLAG_H_

#include "lancet/base/types.h"

namespace lancet::hts {

enum class Strand : bool { FWD = true, REV = false };

// ============================================================================
// SamFlag: SAM specification bitwise flag interpreter.
//
// Wraps the 16-bit FLAG field from a BAM/CRAM record and provides named
// query methods for each flag bit defined in the SAM specification
// (https://samtools.github.io/hts-specs/SAMv1.pdf, §1.4.2).
//
// Standalone — no dependency on Alignment, bam1_t, or HTSlib. This allows
// downstream types (e.g. cbdg::Read) to query flags without pulling in
// the full Alignment header chain.
// ============================================================================
class SamFlag {
 public:
  explicit SamFlag(u16 flags) : mFlag(flags) {}
  SamFlag() = default;

  [[nodiscard]] auto GetStrand() const noexcept -> Strand;
  [[nodiscard]] auto GetMateStrand() const noexcept -> Strand;
  [[nodiscard]] auto IsFwdStrand() const noexcept -> bool;
  [[nodiscard]] auto IsRevStrand() const noexcept -> bool;
  [[nodiscard]] auto IsMateFwdStrand() const noexcept -> bool;
  [[nodiscard]] auto IsMateRevStrand() const noexcept -> bool;
  [[nodiscard]] auto IsQcFail() const noexcept -> bool;
  [[nodiscard]] auto IsDuplicate() const noexcept -> bool;
  [[nodiscard]] auto IsPrimary() const noexcept -> bool;
  [[nodiscard]] auto IsSecondary() const noexcept -> bool;
  [[nodiscard]] auto IsSupplementary() const noexcept -> bool;
  [[nodiscard]] auto IsMapped() const noexcept -> bool;
  [[nodiscard]] auto IsUnmapped() const noexcept -> bool;
  [[nodiscard]] auto IsMateMapped() const noexcept -> bool;
  [[nodiscard]] auto IsMateUnmapped() const noexcept -> bool;
  [[nodiscard]] auto IsPairedInSequencing() const noexcept -> bool;
  [[nodiscard]] auto IsMappedProperPair() const noexcept -> bool;
  [[nodiscard]] auto IsRead1() const noexcept -> bool;
  [[nodiscard]] auto IsRead2() const noexcept -> bool;
  [[nodiscard]] auto HasFlagsSet(u16 check_flags) const noexcept -> bool;
  [[nodiscard]] auto HasFlagsUnset(u16 check_flags) const noexcept -> bool;

  [[nodiscard]] explicit operator u16() const noexcept { return mFlag; }

  auto operator==(SamFlag const& rhs) const -> bool { return mFlag == rhs.mFlag; }
  auto operator!=(SamFlag const& rhs) const -> bool { return !(rhs == *this); }

 private:
  u16 mFlag = 0;
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_SAM_FLAG_H_
