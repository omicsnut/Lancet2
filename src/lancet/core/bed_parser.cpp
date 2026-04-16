#include "lancet/core/bed_parser.h"

#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/hts/reference.h"

#include "absl/cleanup/cleanup.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"

extern "C" {
#include "htslib/hts.h"
#include "htslib/kstring.h"
}

#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::core {

auto ParseBedFile(std::filesystem::path const& bed_file, hts::Reference const& ref)
    -> std::vector<hts::Reference::ParseRegionResult> {
  std::vector<hts::Reference::ParseRegionResult> results;

  if (bed_file.empty()) return results;

  // Open the BED file using HTSlib's htsFile API instead of standard C++ streams.
  // Uses HTSlib's built-in libcurl support to stream BED files from
  // s3://, gs://, and https:// URLs.
  htsFile* fptr = hts_open(bed_file.c_str(), "r");
  if (fptr == nullptr) {
    throw std::runtime_error(fmt::format("Could not open bed file: {}", bed_file.string()));
  }

  absl::Cleanup const stream_cleaner = [fptr, &bed_file]() -> void {
    if (hts_close(fptr) < 0) {
      LOG_WARN("Failed to properly close BED file stream: {}", bed_file.string());
    }
  };

  usize line_num = 0;
  i64 region_start = 0;
  i64 region_end = 0;

  std::string curr_chrom;
  std::vector<std::string_view> tokens;
  tokens.reserve(3);

  kstring_t line = KS_INITIALIZE;
  absl::Cleanup const kstring_cleaner = [&line]() -> void { ks_free(&line); };

  // Stream exactly line-by-line avoiding full-file buffering to minimize memory footprint
  while (hts_getline(fptr, '\n', &line) > 0) {
    line_num++;
    std::string_view const line_view(line.s, line.l);

    if (line_view.starts_with('#') || line_view.empty()) continue;

    tokens = absl::StrSplit(line_view, absl::ByChar('\t'));

    if (tokens.size() != 3) {
      auto const msg = fmt::format("Invalid bed line with {} columns at line number {}",
                                   tokens.size(), line_num);
      throw std::runtime_error(msg);
    }

    if (!absl::SimpleAtoi(tokens[1], &region_start) || !absl::SimpleAtoi(tokens[2], &region_end)) {
      auto const msg =
          fmt::format("Could not parse line {} in bed: {}", line_num, bed_file.filename().string());
      throw std::runtime_error(msg);
    }

    curr_chrom = tokens[0];
    if (!ref.FindChromByName(curr_chrom).ok()) {
      auto const msg = fmt::format("Could not find chrom {} from bed file line {} in reference",
                                   tokens[0], line_num);
      throw std::runtime_error(msg);
    }

    results.emplace_back(hts::Reference::ParseRegionResult{
        .mChromName = curr_chrom, .mRegionSpan = {region_start, region_end}});
  }

  return results;
}

}  // namespace lancet::core
