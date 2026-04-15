#include "lancet/caller/posterior_base_qual.h"

#include "lancet/base/types.h"
#include "lancet/hts/phred_quality.h"

#include "absl/types/span.h"

#include <algorithm>

#include <cmath>

namespace lancet::caller {

auto ComputeRawPosteriorBaseQual(absl::Span<u8 const> fwd_base_quals,
                                 absl::Span<u8 const> rev_base_quals) -> f64 {
  if (fwd_base_quals.empty() && rev_base_quals.empty()) return 0.0;

  f64 log_err = 0.0;  // Σ log10(ε_i)
  f64 log_ok = 0.0;   // Σ log10(1 - ε_i)

  auto const accumulate_bq = [&log_err, &log_ok](absl::Span<u8 const> quals) -> void {
    for (auto const qual : quals) {
      f64 const eps = hts::PhredToErrorProb(qual);
      log_err += std::log10(std::max(eps, 1e-300));
      log_ok += std::log10(std::max(1.0 - eps, 1e-300));
    }
  };

  accumulate_bq(fwd_base_quals);
  accumulate_bq(rev_base_quals);

  // log-sum-exp: log10(10^a + 10^b) = max(a,b) + log10(1 + 10^(min-max))
  f64 const max_log = std::max(log_err, log_ok);
  f64 const log_sum =
      max_log + std::log10(1.0 + std::pow(10.0, std::min(log_err, log_ok) - max_log));
  f64 const log_posterior_err = log_err - log_sum;
  f64 const posterior_bq = -10.0 * log_posterior_err;

  return posterior_bq;
}

}  // namespace lancet::caller
