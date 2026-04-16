#include "lancet/base/hash.h"

#include "lancet/base/types.h"

#include "absl/hash/internal/city.h"

#include <string_view>

namespace lancet::base {

auto HashStr64(std::string_view str) -> u64 {
  // Abseil doesn't expose a public non-Hashable string hash — CityHash is the only option
  // NOLINTNEXTLINE(abseil-no-internal-dependencies)
  return absl::hash_internal::CityHash64(str.data(), str.length());
}

auto HashStr32(std::string_view str) -> u32 {
  // Abseil doesn't expose a public non-Hashable string hash — CityHash is the only option
  // NOLINTNEXTLINE(abseil-no-internal-dependencies)
  return absl::hash_internal::CityHash32(str.data(), str.length());
}

}  // namespace lancet::base
