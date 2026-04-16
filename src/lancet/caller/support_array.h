#ifndef SRC_LANCET_CALLER_SUPPORT_ARRAY_H_
#define SRC_LANCET_CALLER_SUPPORT_ARRAY_H_

#include "lancet/caller/variant_support.h"

#include "absl/container/inlined_vector.h"

#include <memory>
#include <string_view>
#include <utility>

namespace lancet::caller {

// ============================================================================
// SupportArray: per-sample variant evidence container.
//
// Maps sample names to VariantSupport objects. Used by VariantCall to access
// per-sample allele evidence when building VCF FORMAT fields and computing
// genotype likelihoods.
//
// InlinedVector<8> avoids heap allocation for the common case of ≤ 8 samples.
// ============================================================================
class SupportArray {
 public:
  struct NamedSupport {
    // ── 8B Alignment ──────────────────────────────────────────────────────
    std::string_view mSampleName;           // 8B (ptr) + 8B (size) = 16B
    std::unique_ptr<VariantSupport> mData;  // 8B
  };

  [[nodiscard]] auto Find(std::string_view sample_name) const -> VariantSupport const*;
  [[nodiscard]] auto FindOrCreate(std::string_view sample_name) -> VariantSupport&;
  [[nodiscard]] auto Extract(std::string_view sample_name) -> std::unique_ptr<VariantSupport>;
  void Insert(std::string_view sample_name, std::unique_ptr<VariantSupport> data) {
    mItems.emplace_back(sample_name, std::move(data));
  }

  [[nodiscard]] auto begin() const { return mItems.begin(); }
  [[nodiscard]] auto end() const { return mItems.end(); }

 private:
  absl::InlinedVector<NamedSupport, 8> mItems;  // 8B aligned
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_SUPPORT_ARRAY_H_
