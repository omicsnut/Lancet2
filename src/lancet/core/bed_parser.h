#ifndef SRC_LANCET_CORE_BED_PARSER_H_
#define SRC_LANCET_CORE_BED_PARSER_H_

#include "lancet/hts/reference.h"

#include <filesystem>
#include <vector>

namespace lancet::core {

/// Parse all regions from a BED file using HTSlib's hts_getline streaming.
///
/// Supports local files and cloud URIs (s3://, gs://, https://) via HTSlib's
/// built-in libcurl transport. Validates chromosome names against the reference.
///
/// BED format: tab-separated, 3 columns minimum (chrom, start, end).
/// Lines starting with '#' are treated as comments and skipped.
/// Coordinates are 0-based half-open in BED but stored as 1-based closed
/// via ParseRegionResult's OneBasedClosedOptional convention.
[[nodiscard]] auto ParseBedFile(std::filesystem::path const& bed_file, hts::Reference const& ref)
    -> std::vector<hts::Reference::ParseRegionResult>;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_BED_PARSER_H_
