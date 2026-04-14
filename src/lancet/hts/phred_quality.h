#ifndef SRC_LANCET_HTS_PHRED_QUALITY_H_
#define SRC_LANCET_HTS_PHRED_QUALITY_H_

#include "lancet/base/types.h"

namespace lancet::hts {

// Phred quality scores encode sequencing error probability on a logarithmic
// scale: Q = −10 × log₁₀(P_error). A Q30 base has a 1-in-1,000 chance of
// being wrong; Q20 = 1-in-100; Q10 = 1-in-10. Higher Q = more trustworthy.
//
// PhredToErrorProb converts Q back to a probability (e.g., Q30 → 0.001).
// ErrorProbToPhred converts the other direction (e.g., 0.01 → Q20).
//
// Both use the formula P = 10^(−Q/10). PhredToErrorProb uses a precomputed
// constexpr lookup table — thread-safe with zero runtime initialization cost.
static constexpr u8 MAX_PHRED_SCORE = 255;
[[nodiscard]] auto PhredToErrorProb(u32 phred_score) -> f64;
[[nodiscard]] auto ErrorProbToPhred(f64 prob) -> f64;

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_PHRED_QUALITY_H_
