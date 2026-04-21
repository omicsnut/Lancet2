#ifndef SRC_LANCET_CBDG_READ_H_
#define SRC_LANCET_CBDG_READ_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/sam_flag.h"

#include "absl/types/span.h"

#include <numeric>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

class Read {
 public:
  explicit Read(hts::Alignment const& aln, std::string sample_name, Label::Tag const tag,
                usize const sample_index)
      : mStart0(aln.StartPos0()),
        mInsertSize(aln.InsertSize()),
        mSampleIndex(sample_index),
        mQname(aln.QnameView()),
        mSequence(aln.BuildSequence()),
        mSampleName(std::move(sample_name)),
        mQuality(aln.BuildQualities()),
        mChromIdx(aln.ChromIndex()),
        mSamFlag(aln.FlagRaw()),
        mMapQual(aln.MapQual()),
        mTag(tag) {
    static constexpr u8 DEFAULT_MIN_READ_MAP_QUAL = 20;
    if (aln.MapQual() < DEFAULT_MIN_READ_MAP_QUAL) {
      mPassesAlnFilters = false;
    }

    // Flag read as soft-clipped if total soft-clip bases >= 6% of read length.
    // Uses the original whole-genome alignment CIGAR, not re-alignment CIGAR.
    static constexpr f64 SOFT_CLIP_FRAC_THRESHOLD = 0.06;
    static constexpr auto CLIP_LENGTH = [](auto const& unit) -> u32 {
      return unit.Operation() == hts::CigarOp::SOFT_CLIP ? unit.Length() : 0;
    };
    auto const cigar = aln.CigarData();
    auto const total_clip =
        std::transform_reduce(cigar.cbegin(), cigar.cend(), u32{0}, std::plus<>{}, CLIP_LENGTH);
    auto const clip_frac =
        aln.Length() > 0 ? static_cast<f64>(total_clip) / static_cast<f64>(aln.Length()) : 0.0;
    mIsSoftClipped = clip_frac >= SOFT_CLIP_FRAC_THRESHOLD;
  }

  [[nodiscard]] auto StartPos0() const noexcept -> i64 { return mStart0; }
  [[nodiscard]] auto ChromIndex() const noexcept -> i32 { return mChromIdx; }
  [[nodiscard]] auto Flag() const noexcept -> hts::SamFlag { return hts::SamFlag(mSamFlag); }
  [[nodiscard]] auto MapQual() const noexcept -> u8 { return mMapQual; }

  [[nodiscard]] auto PassesAlnFilters() const noexcept -> bool { return mPassesAlnFilters; }

  [[nodiscard]] auto SrcLabel() const noexcept -> Label { return Label(mTag); }
  [[nodiscard]] auto TagKind() const noexcept -> Label::Tag { return mTag; }
  [[nodiscard]] auto QnamePtr() const noexcept -> char const* { return mQname.c_str(); }
  [[nodiscard]] auto SeqPtr() const noexcept -> char const* { return mSequence.c_str(); }
  [[nodiscard]] auto QnameView() const noexcept -> std::string_view { return mQname; }
  [[nodiscard]] auto SeqView() const noexcept -> std::string_view { return mSequence; }
  [[nodiscard]] auto QualView() const noexcept -> absl::Span<u8 const> { return mQuality; }
  [[nodiscard]] auto Length() const noexcept -> usize { return mSequence.size(); }
  [[nodiscard]] auto SampleName() const noexcept -> std::string_view { return mSampleName; }
  [[nodiscard]] auto SampleIndex() const noexcept -> usize { return mSampleIndex; }
  [[nodiscard]] auto IsSoftClipped() const noexcept -> bool { return mIsSoftClipped; }
  [[nodiscard]] auto InsertSize() const noexcept -> i64 { return mInsertSize; }
  [[nodiscard]] auto IsProperPair() const noexcept -> bool { return (mSamFlag & 0x2) != 0; }

  template <typename HashState>
  friend auto AbslHashValue(HashState hash_state, Read const& read) -> HashState {
    return HashState::combine(std::move(hash_state), read.mSampleName, read.mStart0,
                              read.mSampleIndex, read.mChromIdx, read.mSamFlag, read.mMapQual,
                              static_cast<u8>(read.mTag), read.mQname, read.mSequence,
                              read.mQuality);
  }

  friend auto operator==(Read const& lhs, Read const& rhs) noexcept -> bool {
    return lhs.mSampleName == rhs.mSampleName &&
           lhs.mStart0 == rhs.mStart0 &&
           lhs.mSampleIndex == rhs.mSampleIndex &&
           lhs.mChromIdx == rhs.mChromIdx &&
           lhs.mSamFlag == rhs.mSamFlag &&
           lhs.mMapQual == rhs.mMapQual &&
           lhs.mTag == rhs.mTag &&
           lhs.mQname == rhs.mQname &&
           lhs.mSequence == rhs.mSequence &&
           lhs.mQuality == rhs.mQuality &&
           lhs.mPassesAlnFilters == rhs.mPassesAlnFilters;
  }

  friend auto operator!=(Read const& lhs, Read const& rhs) noexcept -> bool {
    return !(lhs == rhs);
  }

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  i64 mStart0 = -1;          // 8B
  i64 mInsertSize = 0;       // 8B
  usize mSampleIndex = 0;    // 8B
  std::string mQname;        // 32B (8B align)
  std::string mSequence;     // 32B (8B align)
  std::string mSampleName;   // 32B (8B align)
  std::vector<u8> mQuality;  // 24B (8B align)
  // ── 4B Align ────────────────────────────────────────────────────────────
  i32 mChromIdx = -1;  // 4B
  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mSamFlag = 0;  // 2B
  // ── 1B Align ────────────────────────────────────────────────────────────
  u8 mMapQual = 0;                // 1B
  bool mPassesAlnFilters = true;  // 1B
  bool mIsSoftClipped = false;    // 1B
  Label::Tag mTag = Label::CTRL;  // 1B
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_READ_H_
