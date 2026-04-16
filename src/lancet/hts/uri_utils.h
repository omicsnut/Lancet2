#ifndef SRC_LANCET_HTS_URI_UTILS_H_
#define SRC_LANCET_HTS_URI_UTILS_H_

#include <string>
#include <string_view>

namespace lancet::hts {

/// Returns true if the path is a cloud/web URI (gs://, s3://, http(s)://, ftp(s)://).
/// Used by CLI validators and pipeline runner to bypass local filesystem checks.
[[nodiscard]] auto IsCloudUri(std::string_view uri) -> bool;

/// Validates read/write access to a cloud URI by opening and immediately closing via hopen/hclose.
///
/// Returns an empty string on success. On failure, returns a human-readable error message
/// suitable for CLI11 validator output or LOG_CRITICAL downstream.
///
/// This catches authentication failures upfront rather than after a 40-hour pipeline run.
/// AWS/GCP use 5MB+ multipart chunk caching over libcurl — small payloads won't trigger
/// the HTTP handshake until BgzfOstream::Close(), so we force an immediate roundtrip here.
[[nodiscard]] auto ValidateCloudAccess(std::string const& uri, std::string const& mode)
    -> std::string;

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_URI_UTILS_H_
