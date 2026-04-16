#ifndef SRC_LANCET_CORE_INPUT_SPEC_PARSER_H_
#define SRC_LANCET_CORE_INPUT_SPEC_PARSER_H_

#include "lancet/cbdg/label.h"

#include "absl/types/span.h"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::core {

/// Parsed result of a single "--sample <path>:<role>" CLI input spec.
struct ParsedInputSpec {
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::filesystem::path mFilePath;  // 8B+ — resolved BAM/CRAM file path
  // ── 1B Align ────────────────────────────────────────────────────────────
  cbdg::Label::Tag mTag = cbdg::Label::CTRL;  // 1B  — sample role (CTRL or CASE)
};

/// Parse a single "--sample <path>:<role>" spec string.
///
/// Splits on the LAST colon and checks if the suffix is a known role
/// ("control" or "case"). This handles cloud URIs (s3://bucket/file.bam)
/// correctly by only splitting on the rightmost colon.
/// If role is omitted or unrecognized, defaults to CTRL.
[[nodiscard]] auto ParseInputSpec(std::string_view spec) -> ParsedInputSpec;

/// Unify all three CLI input modes into a single vector of ParsedInputSpec.
///
/// Handles --normal paths (CTRL), --tumor paths (CASE), and
/// --sample specs (<path>:<role>). No BAM I/O — pure string parsing.
[[nodiscard]] auto ParseAllInputSpecs(absl::Span<std::filesystem::path const> ctrl_paths,
                                      absl::Span<std::filesystem::path const> case_paths,
                                      absl::Span<std::string const> sample_specs)
    -> std::vector<ParsedInputSpec>;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_INPUT_SPEC_PARSER_H_
