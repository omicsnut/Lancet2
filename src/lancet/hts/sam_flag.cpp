#include "lancet/hts/sam_flag.h"

extern "C" {
#include "htslib/sam.h"
}

namespace lancet::hts {

auto SamFlag::GetStrand() const noexcept -> Strand {
  return IsFwdStrand() ? Strand::FWD : Strand::REV;
}

auto SamFlag::GetMateStrand() const noexcept -> Strand {
  return IsMateFwdStrand() ? Strand::FWD : Strand::REV;
}

auto SamFlag::IsFwdStrand() const noexcept -> bool {
  return (mFlag & BAM_FREVERSE) == 0;
}
auto SamFlag::IsRevStrand() const noexcept -> bool {
  return (mFlag & BAM_FREVERSE) != 0;
}
auto SamFlag::IsMateFwdStrand() const noexcept -> bool {
  return (mFlag & BAM_FMREVERSE) == 0;
}
auto SamFlag::IsMateRevStrand() const noexcept -> bool {
  return (mFlag & BAM_FMREVERSE) != 0;
}
auto SamFlag::IsQcFail() const noexcept -> bool {
  return (mFlag & BAM_FQCFAIL) != 0;
}
auto SamFlag::IsDuplicate() const noexcept -> bool {
  return (mFlag & BAM_FDUP) != 0;
}
auto SamFlag::IsPrimary() const noexcept -> bool {
  return (mFlag & BAM_FSECONDARY) == 0;
}
auto SamFlag::IsSecondary() const noexcept -> bool {
  return (mFlag & BAM_FSECONDARY) != 0;
}

auto SamFlag::IsSupplementary() const noexcept -> bool {
  return (mFlag & BAM_FSUPPLEMENTARY) != 0;
}

auto SamFlag::IsMapped() const noexcept -> bool {
  return (mFlag & BAM_FUNMAP) == 0;
}
auto SamFlag::IsUnmapped() const noexcept -> bool {
  return (mFlag & BAM_FUNMAP) != 0;
}
auto SamFlag::IsMateMapped() const noexcept -> bool {
  return (mFlag & BAM_FMUNMAP) == 0;
}
auto SamFlag::IsMateUnmapped() const noexcept -> bool {
  return (mFlag & BAM_FMUNMAP) != 0;
}

auto SamFlag::IsPairedInSequencing() const noexcept -> bool {
  return (mFlag & BAM_FPAIRED) != 0;
}

auto SamFlag::IsMappedProperPair() const noexcept -> bool {
  return (mFlag & BAM_FPROPER_PAIR) != 0;
}

auto SamFlag::IsRead1() const noexcept -> bool {
  return (mFlag & BAM_FREAD1) != 0;
}
auto SamFlag::IsRead2() const noexcept -> bool {
  return (mFlag & BAM_FREAD2) != 0;
}

auto SamFlag::HasFlagsSet(u16 check_flags) const noexcept -> bool {
  return (mFlag & check_flags) != 0;
}

auto SamFlag::HasFlagsUnset(u16 check_flags) const noexcept -> bool {
  return (mFlag & check_flags) == 0;
}

}  // namespace lancet::hts
