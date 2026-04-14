#ifndef SRC_LANCET_BASE_HASH_H_
#define SRC_LANCET_BASE_HASH_H_

#include "lancet/base/types.h"

#include <string_view>

namespace lancet::base {

[[nodiscard]] auto HashStr64(std::string_view str) -> u64;
[[nodiscard]] auto HashStr32(std::string_view str) -> u32;

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_HASH_H_
