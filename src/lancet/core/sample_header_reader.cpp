#include "lancet/core/sample_header_reader.h"

#include "lancet/base/types.h"
#include "lancet/core/input_spec_parser.h"
#include "lancet/core/read_collector.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

#include <absl/types/span.h>
#include <algorithm>
#include <filesystem>
#include <functional>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace lancet::core {

auto MakeSampleList(ReadCollector::Params const& params) -> std::vector<SampleInfo> {
  auto const all_specs =
      ParseAllInputSpecs(params.mCtrlPaths, params.mCasePaths, params.mSampleSpecs);

  std::vector<SampleInfo> results;
  results.reserve(all_specs.size());
  hts::Reference const ref(params.mRefPath);

  // Extract SM read group tag from each BAM/CRAM header → SampleInfo.
  using hts::Alignment::Fields::CORE_QNAME;
  for (auto const& spec : all_specs) {
    hts::Extractor const extractor(spec.mFilePath, ref, CORE_QNAME, {}, true);
    results.emplace_back(extractor.SampleName(), spec.mFilePath, spec.mTag);
  }

  std::ranges::sort(results, std::less<SampleInfo>{});

  // Assign 0-based group indices to identify unique logical samples.
  // Per-lane BAM splits sharing the same SM tag + role get the same index.
  // INVARIANT: group_idx == VCF FORMAT column position for that sample.
  usize group_idx = 0;
  for (usize pos = 0; pos < results.size(); ++pos) {
    if (pos > 0) {
      auto const is_different_name = results[pos].SampleName() != results[pos - 1].SampleName();
      auto const is_different_role = results[pos].TagKind() != results[pos - 1].TagKind();
      if (is_different_name || is_different_role) ++group_idx;
    }
    results[pos].SetSampleIndex(group_idx);
  }

  return results;
}

auto BuildSampleNameList(ReadCollector::Params const& params) -> std::vector<std::string> {
  auto const sinfo_list = MakeSampleList(params);
  std::vector<std::string> results;
  results.reserve(sinfo_list.size());
  std::ranges::transform(
      sinfo_list, std::back_inserter(results),
      [](SampleInfo const& item) -> std::string { return std::string(item.SampleName()); });
  return results;
}

}  // namespace lancet::core
