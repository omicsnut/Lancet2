#include "lancet/hts/phred_quality.h"

#include "lancet/base/types.h"

#include <array>

#include <cmath>

namespace lancet::hts {

namespace {

// ============================================================================
// Phred-to-error-probability lookup table — constexpr generated.
//
// Formula: P = 10^(-Q/10)  where Q is the Phred score (0–255).
//
// Decomposition (avoids non-constexpr std::pow):
//   10^(-Q/10) = 10^(-d) × 10^(-r/10)
//   where d = Q ÷ 10  (integer part)
//         r = Q % 10  (remainder)
//
//   10^(-d) is computed as 1 / (10^d), where 10^d uses exact × 10.0 steps.
//   10^(-r/10) uses 10 precomputed fractional constants (r = 0..9).
//
// Precision: 10^N is exact in IEEE 754 for N ≤ 22 (5^22 fits in 53-bit
// mantissa). For Q ≤ 229 (int_exp ≤ 22), each entry is mathematically exact.
// For Q = 230–255 (int_exp = 23–25), the single division introduces at most
// 1–2 ULP deviation — negligible for Phred scoring.
// ============================================================================
constexpr auto MakePhredToErrorProbTable() -> std::array<f64, MAX_PHRED_SCORE + 1> {
  // 10^(-r/10) for r = 0..9. These are the IEEE 754 doubles closest to the
  // exact mathematical values. Verifiable: any calculator with 10^(-r/10).
  //
  //   r=0: 10^( 0.0) = 1.0
  //   r=1: 10^(-0.1) ≈ 0.7943...
  //   r=2: 10^(-0.2) ≈ 0.6309...
  //   r=3: 10^(-0.3) ≈ 0.5011...
  //   r=4: 10^(-0.4) ≈ 0.3981...
  //   r=5: 10^(-0.5) ≈ 0.3162...
  //   r=6: 10^(-0.6) ≈ 0.2511...
  //   r=7: 10^(-0.7) ≈ 0.1995...
  //   r=8: 10^(-0.8) ≈ 0.1584...
  //   r=9: 10^(-0.9) ≈ 0.1258...
  constexpr std::array<f64, 10> FRAC_FACTOR = {
      1.0,
      0.7943282347242815020659182828,  // 10^(-1/10)
      0.6309573444801932494343601366,  // 10^(-2/10)
      0.5011872336272722850015541869,  // 10^(-3/10)
      0.3981071705534972507702523051,  // 10^(-4/10)
      0.3162277660168379331998893544,  // 10^(-5/10)
      0.2511886431509580111085032068,  // 10^(-6/10)
      0.1995262314968879601352455397,  // 10^(-7/10)
      0.1584893192461113485202101373,  // 10^(-8/10)
      0.1258925411794167210423954106,  // 10^(-9/10)
  };

  std::array<f64, MAX_PHRED_SCORE + 1> result{};
  for (u32 phred = 0; phred <= MAX_PHRED_SCORE; ++phred) {
    // 10^(-Q/10) = FRAC_FACTOR[frac_idx] / 10^(int_exp)
    auto const int_exp = phred / 10U;
    auto const frac_idx = phred % 10U;

    // Multiply by 10.0 (exact in IEEE 754) instead of 0.1 (inexact repeating
    // binary fraction). Since 10^N = 2^N × 5^N, the 2^N is absorbed by the
    // exponent field and only 5^N occupies the 53-bit mantissa. 5^22 fits in
    // 53 bits, so 10^0 through 10^22 are computed with zero rounding error.
    // For Q ≤ 229 (int_exp ≤ 22): mathematically exact.
    // For Q = 230–255 (int_exp = 23–25): at most 1–2 ULP from the division.
    f64 pow10 = 1.0;
    for (u32 step = 0; step < int_exp; ++step) pow10 *= 10.0;

    result[phred] = FRAC_FACTOR[frac_idx] / pow10;
  }
  return result;
}

}  // namespace

static constexpr auto LUT_PHRED_TO_ERROR_PROB = MakePhredToErrorProbTable();

auto PhredToErrorProb(u32 phred_score) -> f64 {
  auto const idx =
      phred_score > MAX_PHRED_SCORE ? MAX_PHRED_SCORE : static_cast<usize>(phred_score);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  return LUT_PHRED_TO_ERROR_PROB[idx];  // idx already clamped to [0, MAX_PHRED_SCORE]
}

auto ErrorProbToPhred(f64 prob) -> f64 {
  if (prob == 1.0) {
    return 0;
  }
  if (prob == 0.0) {
    return MAX_PHRED_SCORE;
  }
  static constexpr f64 PHRED_MULTIPLIER = -10.0;
  return PHRED_MULTIPLIER * std::log10(prob);
}

}  // namespace lancet::hts
