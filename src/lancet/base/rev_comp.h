#ifndef SRC_LANCET_BASE_REV_COMP_H_
#define SRC_LANCET_BASE_REV_COMP_H_

#include "lancet/base/types.h"

#include <algorithm>
#include <array>
#include <ranges>
#include <string>
#include <string_view>

namespace lancet::base {

// Constexpr lookup table for DNA complement: A↔T, C↔G, else→N.
constexpr auto MakeDnaComplementTable() -> std::array<char, 256> {
  std::array<char, 256> tbl{};
  for (auto& val : tbl) {
    val = 'N';
  }
  tbl['A'] = 'T';
  tbl['a'] = 't';
  tbl['T'] = 'A';
  tbl['t'] = 'a';
  tbl['C'] = 'G';
  tbl['c'] = 'g';
  tbl['G'] = 'C';
  tbl['g'] = 'c';
  tbl['N'] = 'N';
  tbl['n'] = 'n';
  return tbl;
}

inline constexpr std::array<char, 256> DNA_COMPLEMENT_TABLE = MakeDnaComplementTable();

[[nodiscard]] inline auto RevComp(char const& base) -> char {
  return DNA_COMPLEMENT_TABLE[static_cast<u8>(base)];
}

[[nodiscard]] inline auto RevComp(std::string_view seq) -> std::string {
  std::string result(seq.size(), 'N');
  std::ranges::transform(std::views::reverse(seq), result.begin(), [](char const base) {
    return DNA_COMPLEMENT_TABLE[static_cast<u8>(base)];
  });
  return result;
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_REV_COMP_H_
