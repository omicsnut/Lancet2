#include "lancet/core/input_spec_parser.h"

#include "lancet/cbdg/label.h"

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

// Returns the graph coloring tag for a known role string,
// or std::nullopt if the string is not a recognized role.
auto TryParseRole(std::string_view const role) -> std::optional<lancet::cbdg::Label::Tag> {
  if (role == "control") return lancet::cbdg::Label::CTRL;
  if (role == "case") return lancet::cbdg::Label::CASE;
  return std::nullopt;
}

}  // namespace

namespace lancet::core {

auto ParseInputSpec(std::string_view spec) -> ParsedInputSpec {
  // Check if the suffix after the LAST colon is a known role string
  // ("control" or "case"). If it is, split there. If not, the entire
  // string is the path and the tag defaults to Label::CTRL.
  // This handles cloud URIs (s3://bucket/file.bam) without false splits.
  auto const colon_pos = spec.rfind(':');
  auto const suffix = (colon_pos != std::string_view::npos && colon_pos < spec.size() - 1)
                          ? spec.substr(colon_pos + 1)
                          : std::string_view{};
  auto const parsed_tag = TryParseRole(suffix);

  auto const fpath = parsed_tag.has_value() ? std::filesystem::path(spec.substr(0, colon_pos))
                                            : std::filesystem::path(spec);
  auto const tag = parsed_tag.value_or(cbdg::Label::CTRL);
  return {.mFilePath = fpath, .mTag = tag};
}

auto ParseAllInputSpecs(absl::Span<std::filesystem::path const> ctrl_paths,
                        absl::Span<std::filesystem::path const> case_paths,
                        absl::Span<std::string const> sample_specs)
    -> std::vector<ParsedInputSpec> {
  std::vector<ParsedInputSpec> results;
  results.reserve(ctrl_paths.size() + case_paths.size() + sample_specs.size());

  // Legacy --normal paths → control samples
  std::ranges::transform(ctrl_paths, std::back_inserter(results),
                         [](auto const& fpath) -> ParsedInputSpec {
                           return {.mFilePath = fpath, .mTag = cbdg::Label::CTRL};
                         });

  // Legacy --tumor paths → case samples
  std::ranges::transform(case_paths, std::back_inserter(results),
                         [](auto const& fpath) -> ParsedInputSpec {
                           return {.mFilePath = fpath, .mTag = cbdg::Label::CASE};
                         });

  // Unified --sample specs (<path>:<role>)
  std::ranges::transform(sample_specs, std::back_inserter(results),
                         [](auto const& spec) { return ParseInputSpec(spec); });

  return results;
}

}  // namespace lancet::core
