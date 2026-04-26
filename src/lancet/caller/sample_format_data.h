#ifndef SRC_LANCET_CALLER_SAMPLE_FORMAT_DATA_H_
#define SRC_LANCET_CALLER_SAMPLE_FORMAT_DATA_H_

#include "lancet/base/types.h"

#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"

#include <array>
#include <optional>
#include <string>
#include <utility>

namespace lancet::caller {

// ============================================================================
// SampleFormatData: per-sample VCF FORMAT field payload.
//
// Holds all FORMAT-level annotations for one sample in a VariantCall.
// Members ordered by alignment: 8B → 4B → 2B → 1B.
//
// MEMORY COMPACTION (optional fields):
// The 10 optional FORMAT metrics (FLD, RPCD, BQCD, SDFC, etc.) are stored as a
// compact bitfield + flat array instead of 10 × std::optional<f32>.
//   Before: 10 × 8B (std::optional<f32>) = 80B
//   After:  2B (u16 mask) + 40B (f32[10]) = 42B  (padded to 44B)
// Savings: 36B per sample per variant call.
//
// Cannot use IEEE NaN because -ffast-math (set in cmake/defaults.cmake) implies
// -ffinite-math-only, which lets the compiler assume NaN never exists — so
// std::isnan() is optimized to false and quiet_NaN() may be eliminated.
// ============================================================================
class SampleFormatData {
 public:
  // ── Optional field indices ────────────────────────────────────────────
  // NOLINTNEXTLINE(cppcoreguidelines-use-enum-class)
  enum FieldSlot : u8 {
    FRAG_LEN_DELTA = 0,          // FLD
    READ_POS_COHEN_D = 1,        // RPCD
    BASE_QUAL_COHEN_D = 2,       // BQCD
    MAP_QUAL_COHEN_D = 3,        // MQCD
    ALLELE_MISMATCH_DELTA = 4,   // ASMD
    FRAG_START_ENTROPY = 5,      // FSSE
    ALT_HAP_DISCORD_DELTA = 6,   // AHDD
    HAPLOTYPE_SEG_ENTROPY = 7,   // HSE
    PATH_DEPTH_CV = 8,           // PDCV
    SITE_DEPTH_FOLD_CHANGE = 9,  // SDFC
    NUM_FIELDS = 10,
  };

  // ============================================================================
  // Accessors
  // ============================================================================
  // SetField bridges VariantSupport optional<f64> → compact f32 storage.
  // Narrows f64→f32 and sets the presence bit atomically.
  void SetField(FieldSlot slot, std::optional<f64> const& value) {
    if (value.has_value()) {
      mFieldValues[slot] = static_cast<f32>(value.value());
      mFieldPresence |= static_cast<u16>(1U << slot);
    }
  }

  void SetStrandBias(f32 value) { mStrandBias = value; }
  void SetSoftClipAsym(f32 value) { mSoftClipAsym = value; }
  void SetPolarRadius(f32 value) { mPolarRadius = value; }
  void SetPolarAngle(f32 value) { mPolarAngle = value; }
  void SetTotalDepth(u32 value) { mTotalDepth = value; }
  void SetGenotypeQuality(u32 value) { mGenotypeQuality = value; }
  void SetGenotypeIndices(i16 first, i16 second) { mGenotypeIndices = {first, second}; }
  void SetMissingSupport(bool value) { mIsMissingSupport = value; }

  void SetPhredLikelihoods(absl::InlinedVector<u32, 6> pls) { mPhredLikelihoods = std::move(pls); }
  void SetAlleleDepths(absl::InlinedVector<u16, 4> depths) { mAlleleDepths = std::move(depths); }
  void SetFwdAlleleDepths(absl::InlinedVector<u16, 4> depths) {
    mFwdAlleleDepths = std::move(depths);
  }
  void SetRevAlleleDepths(absl::InlinedVector<u16, 4> depths) {
    mRevAlleleDepths = std::move(depths);
  }
  void SetRmsMappingQualities(absl::InlinedVector<f32, 4> quals) {
    mRmsMappingQualities = std::move(quals);
  }
  void SetNormPosteriorBQs(absl::InlinedVector<f32, 4> quals) {
    mNormPosteriorBQs = std::move(quals);
  }
  void SetContinuousMixtureLods(absl::InlinedVector<f64, 4> lods) {
    mContinuousMixtureLods = std::move(lods);
  }

  // ── Getters ──────────────────────────────────────────────────────────
  [[nodiscard]] auto GetField(FieldSlot slot) const -> std::optional<f32> {
    if ((mFieldPresence & static_cast<u16>(1U << slot)) == 0) return std::nullopt;
    return mFieldValues[slot];
  }

  [[nodiscard]] auto HasField(FieldSlot slot) const -> bool {
    return (mFieldPresence & static_cast<u16>(1U << slot)) != 0;
  }

  [[nodiscard]] auto StrandBias() const -> f32 { return mStrandBias; }
  [[nodiscard]] auto SoftClipAsym() const -> f32 { return mSoftClipAsym; }
  [[nodiscard]] auto PolarRadius() const -> f32 { return mPolarRadius; }
  [[nodiscard]] auto PolarAngle() const -> f32 { return mPolarAngle; }
  [[nodiscard]] auto TotalDepth() const -> u32 { return mTotalDepth; }
  [[nodiscard]] auto GenotypeQuality() const -> u32 { return mGenotypeQuality; }
  [[nodiscard]] auto GenotypeIndices() const -> std::pair<i16, i16> { return mGenotypeIndices; }
  [[nodiscard]] auto IsMissingSupport() const -> bool { return mIsMissingSupport; }

  [[nodiscard]] auto PhredLikelihoods() const -> absl::Span<u32 const> { return mPhredLikelihoods; }
  [[nodiscard]] auto AlleleDepths() const -> absl::Span<u16 const> { return mAlleleDepths; }
  [[nodiscard]] auto FwdAlleleDepths() const -> absl::Span<u16 const> { return mFwdAlleleDepths; }
  [[nodiscard]] auto RevAlleleDepths() const -> absl::Span<u16 const> { return mRevAlleleDepths; }
  [[nodiscard]] auto RmsMappingQualities() const -> absl::Span<f32 const> {
    return mRmsMappingQualities;
  }
  [[nodiscard]] auto NormPosteriorBQs() const -> absl::Span<f32 const> { return mNormPosteriorBQs; }
  [[nodiscard]] auto ContinuousMixtureLods() const -> absl::Span<f64 const> {
    return mContinuousMixtureLods;
  }

  /// VCF string rendering — implemented in sample_format_data.cpp
  [[nodiscard]] auto RenderVcfString() const -> std::string;

 private:
  // ── 8B Align ────────────────────────────────────────────────────────────
  absl::InlinedVector<u32, 6> mPhredLikelihoods;  // biallelic=3, triallelic=6 inline
  absl::InlinedVector<u16, 4> mAlleleDepths;
  absl::InlinedVector<u16, 4> mFwdAlleleDepths;
  absl::InlinedVector<u16, 4> mRevAlleleDepths;
  absl::InlinedVector<f32, 4> mRmsMappingQualities;
  absl::InlinedVector<f32, 4> mNormPosteriorBQs;
  absl::InlinedVector<f64, 4> mContinuousMixtureLods;

  // ── 4B Align ────────────────────────────────────────────────────────────
  f32 mStrandBias{0.0F};
  f32 mSoftClipAsym{0.0F};
  std::array<f32, NUM_FIELDS> mFieldValues{};  // 40B — dense, cache-friendly
  f32 mPolarRadius{0.0F};
  f32 mPolarAngle{0.0F};
  u32 mTotalDepth{0};
  u32 mGenotypeQuality{0};

  // ── 2B Align ────────────────────────────────────────────────────────────
  u16 mFieldPresence{0};  // bitfield — which FieldSlots are valid
  std::pair<i16, i16> mGenotypeIndices = {-1, -1};

  // ── 1B Align ────────────────────────────────────────────────────────────
  bool mIsMissingSupport = false;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_SAMPLE_FORMAT_DATA_H_
