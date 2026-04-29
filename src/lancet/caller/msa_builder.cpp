#include "lancet/caller/msa_builder.h"

#include "lancet/base/types.h"

#include "absl/types/span.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

/*
 * NOTE: Both the AlignmentEngine and the spoa::Graph are instantiated once per
 * worker thread within the caller::MsaBuilder object in VariantBuilder. They are
 * passed here by reference. The engine maintains the scoring params and preallocated
 * DP matrix buffer. The graph, which stores the topological biological state,
 * is explicitly Clear()ed per window component. This eliminates the massive OS heap
 * fragmentation penalty from continuous graph initialization/destruction, ensuring
 * strict zero-allocation maximum execution speed across millions of genomic windows.
 */

namespace lancet::caller {

void MsaBuilder::UpdateSpoaState(absl::Span<std::string_view const> sequences,
                                 absl::Span<cbdg::Path::BaseWeights const> weights) {
  mGraph.Clear();
  for (usize idx = 0; idx < sequences.size(); ++idx) {
    // Views point into Path::mSequence strings owned by ComponentResult on the
    // caller's stack — guaranteed to outlive this call. SPOA's const char*
    // overloads take an explicit uint32_t length, so null-termination is irrelevant.
    // NOLINTBEGIN(bugprone-suspicious-stringview-data-usage)
    auto const seq_len = static_cast<std::uint32_t>(sequences[idx].size());
    auto const alignment = mEngine->Align(sequences[idx].data(), seq_len, mGraph);
    mGraph.AddAlignment(alignment, sequences[idx].data(), seq_len, weights[idx]);
    // NOLINTEND(bugprone-suspicious-stringview-data-usage)
  }
}

auto MsaBuilder::BuildGfaString() const -> std::string {
  // https://github.com/rvaser/spoa/pull/36
  // See PR for how to normalize & process the output GFA
  fmt::memory_buffer gfa_buffer;
  fmt::format_to(std::back_inserter(gfa_buffer), "H\tVN:Z:1.0\n");

  // --- 1. Write Segments (Nodes) and Links (Edges) ---
  for (std::unique_ptr<spoa::Graph::Node> const& node : mGraph.nodes()) {
    auto const node_id = node->id + 1;
    auto const node_seq_char = static_cast<char>(mGraph.decoder(node->code));
    fmt::format_to(std::back_inserter(gfa_buffer), "S\t{}\t{}\n", node_id, node_seq_char);

    for (spoa::Graph::Edge const* edge : node->outedges) {
      auto const dest_node_id = edge->head->id + 1;
      fmt::format_to(std::back_inserter(gfa_buffer), "L\t{}\t+\t{}\t+\t0M\n", node_id,
                     dest_node_id);
    }
  }

  // --- 2. Write Walks (Sequence Paths) ---
  for (u32 seq_idx = 0; seq_idx < mGraph.sequences().size(); ++seq_idx) {
    auto const* const hap_prefix = seq_idx == 0 ? "ref" : "hap";
    fmt::format_to(std::back_inserter(gfa_buffer), "P\t{}{}\t", hap_prefix, seq_idx);

    // Trace node successors to construct the sequence path.
    spoa::Graph::Node const* curr_node = mGraph.sequences()[seq_idx];
    bool is_first_node_in_path = true;
    while (curr_node != nullptr) {
      auto const* const node_separator = is_first_node_in_path ? "" : ",";
      fmt::format_to(std::back_inserter(gfa_buffer), "{}{}+", node_separator, curr_node->id + 1);
      is_first_node_in_path = false;
      curr_node = curr_node->Successor(seq_idx);
    }

    fmt::format_to(std::back_inserter(gfa_buffer), "\t*\n");
  }

  return fmt::to_string(gfa_buffer);
}

auto MsaBuilder::BuildFastaString(absl::Span<std::string const> msa_alns) -> std::string {
  fmt::memory_buffer fasta_buffer;
  for (usize idx = 0; idx < msa_alns.size(); ++idx) {
    auto const* const sequence_role_prefix = idx == 0 ? "ref" : "hap";
    fmt::format_to(std::back_inserter(fasta_buffer), ">{}{}\n{}\n", sequence_role_prefix, idx,
                   msa_alns[idx]);
  }
  return fmt::to_string(fasta_buffer);
}

}  // namespace lancet::caller
