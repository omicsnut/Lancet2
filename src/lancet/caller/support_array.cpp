#include "lancet/caller/support_array.h"

#include "lancet/caller/variant_support.h"

#include <absl/container/inlined_vector.h>
#include <algorithm>
#include <memory>
#include <ranges>
#include <string_view>

namespace lancet::caller {

auto SupportArray::Find(std::string_view sample_name) const -> VariantSupport const* {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto const* const iter = std::ranges::find_if(mItems, pred);
  return iter != mItems.end() ? iter->mData.get() : nullptr;
}

auto SupportArray::FindOrCreate(std::string_view sample_name) -> VariantSupport& {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto* const iter = std::ranges::find_if(mItems, pred);
  if (iter != mItems.end()) return *iter->mData;

  mItems.push_back(
      NamedSupport{.mSampleName = sample_name, .mData = std::make_unique<VariantSupport>()});

  return *mItems.back().mData;
}

auto SupportArray::Extract(std::string_view sample_name) -> std::unique_ptr<VariantSupport> {
  auto const pred = [&](auto const& item) -> bool { return item.mSampleName == sample_name; };
  auto* const iter = std::ranges::find_if(mItems, pred);
  if (iter != mItems.end()) {
    auto extracted_data = std::move(iter->mData);
    mItems.erase(iter);
    return extracted_data;
  }

  return nullptr;
}

}  // namespace lancet::caller
