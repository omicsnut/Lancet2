#ifndef SRC_LANCET_CORE_SAMPLE_HEADER_READER_H_
#define SRC_LANCET_CORE_SAMPLE_HEADER_READER_H_

#include "lancet/core/read_collector.h"
#include "lancet/core/sample_info.h"

#include <string>
#include <vector>

namespace lancet::core {

/// Build sorted, indexed sample list from CLI params.
/// Opens BAM headers to extract SM read group tags, sorts by (TagKind, SampleName),
/// and assigns deterministic sample indices for VCF FORMAT column ordering.
[[nodiscard]] auto MakeSampleList(ReadCollector::Params const& params) -> std::vector<SampleInfo>;

/// Build sample name list for VCF header #CHROM line.
/// Opens BAM headers to extract SM read group tags.
[[nodiscard]] auto BuildSampleNameList(ReadCollector::Params const& params)
    -> std::vector<std::string>;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_SAMPLE_HEADER_READER_H_
