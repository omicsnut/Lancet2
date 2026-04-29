#ifndef SRC_LANCET_CALLER_MSA_BUILDER_H_
#define SRC_LANCET_CALLER_MSA_BUILDER_H_

#include "lancet/cbdg/path.h"

#include "absl/types/span.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <memory>
#include <string>
#include <string_view>

namespace lancet::caller {

class MsaBuilder {
 public:
  // ── 8B Align ────────────────────────────────────────────────────────────
  std::unique_ptr<spoa::AlignmentEngine> mEngine;
  spoa::Graph mGraph;

  void UpdateSpoaState(absl::Span<std::string_view const> sequences,
                       absl::Span<cbdg::Path::BaseWeights const> weights);

  /// Render the SPOA alignment graph as a GFA-1.0 document. Caller decides
  /// where the bytes go (typically enqueued for the background flusher).
  [[nodiscard]] auto BuildGfaString() const -> std::string;

  /// Render the multiple sequence alignment as a FASTA document. The
  /// `msa_alns` span comes from `mGraph.GenerateMultipleSequenceAlignment`.
  [[nodiscard]] static auto BuildFastaString(absl::Span<std::string const> msa_alns) -> std::string;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_MSA_BUILDER_H_
