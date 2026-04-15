#ifndef SRC_LANCET_CALLER_VARIANT_SET_H_
#define SRC_LANCET_CALLER_VARIANT_SET_H_

#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

#include "absl/container/btree_set.h"

#include <string_view>

namespace spoa {
class Graph;
}  // namespace spoa

namespace lancet::caller {

// ============================================================================
// VariantSet: Multiallelic Graph Extraction Engine
//
// Extracts multiallelic variants from SPOA directed acyclic graphs by sweeping
// the topology per haplotype. Tracks divergent paths and merges them into
// bundled `RawVariant` outputs with no overlapping biases.
// ============================================================================
class VariantSet {
 public:
  VariantSet(spoa::Graph const& graph, core::Window const& win, usize ref_anchor_start);

  using BTree = absl::btree_set<RawVariant>;

  [[nodiscard]] auto begin() -> BTree::iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto begin() const -> BTree::const_iterator { return mResultVariants.begin(); }
  [[nodiscard]] auto cbegin() const -> BTree::const_iterator { return mResultVariants.cbegin(); }

  [[nodiscard]] auto end() -> BTree::iterator { return mResultVariants.end(); }
  [[nodiscard]] auto end() const -> BTree::const_iterator { return mResultVariants.end(); }
  [[nodiscard]] auto cend() const -> BTree::const_iterator { return mResultVariants.cend(); }

  [[nodiscard]] auto rbegin() -> BTree::reverse_iterator { return mResultVariants.rbegin(); }
  [[nodiscard]] auto rbegin() const -> BTree::const_reverse_iterator {
    return mResultVariants.rbegin();
  }
  [[nodiscard]] auto crbegin() const -> BTree::const_reverse_iterator {
    return mResultVariants.crbegin();
  }

  [[nodiscard]] auto rend() -> BTree::reverse_iterator { return mResultVariants.rend(); }
  [[nodiscard]] auto rend() const -> BTree::const_reverse_iterator {
    return mResultVariants.rend();
  }
  [[nodiscard]] auto crend() const -> BTree::const_reverse_iterator {
    return mResultVariants.crend();
  }

  [[nodiscard]] auto IsEmpty() const -> bool { return mResultVariants.empty(); }
  [[nodiscard]] auto Count() const -> usize { return mResultVariants.size(); }

 private:
  absl::btree_set<RawVariant> mResultVariants;
};

}  // namespace lancet::caller

#endif  // SRC_LANCET_CALLER_VARIANT_SET_H_
