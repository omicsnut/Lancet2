#ifndef SRC_LANCET_CLI_VCF_HEADER_BUILDER_H_
#define SRC_LANCET_CLI_VCF_HEADER_BUILDER_H_

#include "lancet/cli/cli_params.h"

#include <string>

namespace lancet::cli {

/// Builds a complete VCF v4.5 header string from CLI parameters.
/// Opens the reference FASTA to enumerate contigs and conditionally
/// emits SHARED/CTRL/CASE INFO lines when case-control mode is active.
[[nodiscard]] auto BuildVcfHeader(CliParams const& params) -> std::string;

}  // namespace lancet::cli

#endif  // SRC_LANCET_CLI_VCF_HEADER_BUILDER_H_
