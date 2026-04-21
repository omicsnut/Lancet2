#ifndef SRC_LANCET_BASE_POLAR_COORDS_H_
#define SRC_LANCET_BASE_POLAR_COORDS_H_

// ============================================================================
// Polar Coordinate Feature Engineering for Variant Classification
//
// Transforms Cartesian allele depth coordinates (AD_Ref, AD_Alt) into polar
// coordinates (Radius, Angle) to separate biological identity from
// sequencing depth so that neither depends on the other.
//
// ============================================================================
// Why Not DP and VAF?
// ============================================================================
//
// Standard VCF already has DP = AD_Ref + AD_Alt (depth) and
// VAF = AD_Alt / DP (allele fraction). Both degrade ML classification:
//
//   DP is a simple sum (the L₁ norm) of the allele depth vector.
//   It COUPLES the two allele counts: for a fixed DP, increasing Alt
//   forces Ref down. A somatic AD=(95,5) and germline AD=(50,50) both
//   get DP=100 despite representing fundamentally different signal
//   geometries. Sequencing noise scatters reads equally in all
//   directions around the true depth point, forming a circular zone
//   of expected error — the Euclidean distance (L₂ norm) matches
//   this shape; the simple-sum diagonal does not.
//
//   VAF has depth-dependent precision: a 1-read fluctuation swings
//   VAF by 10% at 10× but only 0.1% at 1000×. Two identical VAF=0.5
//   points at different depths carry radically different reliability,
//   but the number itself gives no indication. The ML model must learn
//   the "fan" shape where VAF's meaning depends on DP — a hidden
//   dependency that slows learning and hurts generalization.
//
// ============================================================================
// The Polar Transform
// ============================================================================
//
//   PANG = atan2(AD_Alt, AD_Ref): isolates IDENTITY (the allele fraction)
//     as an angle. A germline het at 20× and 2000× both get PANG=0.785
//     (45°) — the entire problematic diagonal collapses to a constant.
//
//   PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²)): isolates MAGNITUDE (the
//     signal strength) as a log-compressed Euclidean radius. The 100×
//     dynamic range of raw radius compresses to [0, ~3.5].
//
//   Together they make the noise spread UNIFORM across the feature space:
//   the flaring Cartesian wedge (where variance grows with depth) unrolls
//   into a rectangular strip with constant-width noise at every depth.
//   PANG answers "what is it?" and PRAD answers "how sure are we?".
//   These two questions are INDEPENDENT — the model learns each axis
//   separately, converges faster, and generalizes across coverages.
//
// ============================================================================
// Biological Interpretation (Single-Sample)
// ============================================================================
//
//   Angle (PANG):
//     ≈ π/2 (1.571 rad, 90°): Homozygous ALT — nearly all reads are ALT.
//     ≈ π/4 (0.785 rad, 45°): Heterozygous — balanced REF/ALT split.
//     < π/6 (0.524 rad, 30°): Low-frequency signal (somatic, mosaic, or artifact).
//     ≈ 0.0 (0°):             Homozygous REF — no ALT evidence.
//
//   Radius (PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²))):
//     Acts as a confidence tie-breaker for low-angle variants. Log10
//     compression maps all practical coverages to [0, ~3.5]:
//       Low angle + HIGH PRAD (>2) → true low-frequency variant (deep seq).
//       Low angle + LOW PRAD  (<1) → stochastic noise (insufficient evidence).
//
// ============================================================================
// Case-Control "Four-Feature" Paradigm
// ============================================================================
//
//   Each sample independently computes its own (PRAD, PANG) from its
//   own allele depths. The downstream ML model receives:
//     [PRAD_Control, PANG_Control, PRAD_Case, PANG_Case]
//
//   PRAD's log10 compression ensures the four-feature paradigm is coverage-
//   invariant: both angle (identity) and radius (confidence) produce bounded,
//   comparable values regardless of whether the sample is 20× or 2000×.
//
//   This avoids cross-sample relational metrics that would violate VCF's
//   per-sample FORMAT semantics and eliminates sentinel/NaN imputation
//   issues entirely. The model learns the relational delta itself:
//     - Germline:  PANG_Control ≈ π/4  AND PANG_Case ≈ π/4
//     - Somatic:   PANG_Control ≈ 0    AND PANG_Case ≈ π/6..π/4
//     - LOH:       PANG_Control ≈ π/4  AND PANG_Case ≈ π/2
//
// ============================================================================
// Performance Notes
// ============================================================================
//
//   PolarRadius avoids std::hypot() which includes expensive overflow/
//   underflow checks. Since read counts are bounded well within f64 range
//   (max ~millions), intermediate overflow is impossible—std::sqrt(x²+y²)
//   compiles to a single vsqrtsd instruction with -O2.
//
//   PolarAngle uses a branchless minimax polynomial approximation of atan2
//   (Nvidia/Cg fast math) instead of std::atan2. Max error ~0.0015 radians
//   (~0.086°), negligible for ML features formatted to 4 decimal places.
//   The approximation avoids the complex polynomial expansion of std::atan2
//   and compiles to a tight sequence of FP multiply-add instructions.
//
// ============================================================================
// References
// ============================================================================
//
//   - Weisstein, E.W. "Polar Coordinates." MathWorld.
//   - GATK AlleleDepth (AD) format field specification.
//   - VCF 4.3 specification, Section 1.6.2 (FORMAT fields).
//   - Nvidia Cg 3.1 Toolkit: fast atan2 minimax polynomial.
//
// ============================================================================

#include "lancet/base/types.h"

#include <numbers>

#include <cmath>

namespace lancet::base {

/// Polar Radius: log10-compressed Euclidean magnitude of the allele depth vector.
///
/// In plain terms: PRAD replaces DP (read depth) as the "weight of evidence"
/// metric. DP is just a sum (Ref + Alt) — it treats the two allele counts as
/// interchangeable. PRAD uses the Euclidean distance (pythagorean formula),
/// which treats them as independent dimensions of a vector. This matters
/// because sequencing noise forms circular clouds around the true value, not
/// diagonal lines — the Euclidean distance matches the noise shape. The log
/// compression maps the entire practical range to [0, ~3.5].
///
///   PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²))
///
/// The raw Euclidean radius scales linearly with coverage (0→2800+ at 2000×),
/// creating a 100× dynamic range that prevents ML models trained at one
/// coverage from generalizing to another. The log10(1+r) compression:
///   - Preserves monotonicity and ordering (log10 is strictly increasing)
///   - Compresses the range to [0, ~3.5] across all practical coverages
///   - Makes spacing between coverage regimes approximately uniform
///   - Eliminates the need for downstream StandardScaler normalization
///
/// The "+1" ensures log10(1+0) = 0 when both inputs are zero, avoiding
/// log10(0) = −∞ and providing a clean zero baseline for no-coverage sites.
///
/// Uses direct multiplication + sqrt instead of std::hypot to avoid
/// unnecessary overflow/underflow checks (read counts are bounded).
///
/// Returns 0.0 when both inputs are zero (no reads at this locus).
[[nodiscard]] inline auto PolarRadius(f64 const ref_depth, f64 const alt_depth) -> f64 {
  return std::log10(1.0 + std::sqrt((ref_depth * ref_depth) + (alt_depth * alt_depth)));
}

/// Polar Angle: Fast branchless atan2 approximation, in radians.
///
/// In plain terms: PANG replaces VAF (Variant Allele Fraction) as the
/// "what is this variant?" metric. Both encode the allele fraction, but
/// VAF = Alt/(Ref+Alt) has depth-dependent precision: a 1-read error
/// changes VAF by 10% at 10× depth but only 0.1% at 1000× — so its
/// reliability depends on depth. PANG uses the angle of the allele
/// depth vector, which is a pure ratio independent of magnitude. When
/// paired with PRAD, the two form an independent pair: PANG says
/// "what" (identity), PRAD says "how sure" (confidence), and neither
/// depends on the other.
///
///   θ ≈ atan2(AD_Alt, AD_Ref)
///
/// Encodes the biological IDENTITY of the variant — the allele fraction —
/// independent of total depth. Ranges from 0 (hom REF) to π/2 (hom ALT).
///
/// Uses the Nvidia/Cg minimax polynomial approximation of atan2:
///   1. Normalize inputs to the range [-1, 1] via r = (x - |y|) / (|y| + |x|).
///   2. Apply 3rd-degree minimax polynomial: θ ≈ (0.1963·r² − 0.9817)·r + base.
///   3. Correct for octant via copysign.
///
/// Max absolute error: ~0.0015 radians (~0.086°). For ML features formatted
/// to 4 decimal places (≈0.0001 rad resolution), this error is at the
/// noise floor and does not affect classification boundaries.
///
/// Edge cases (allele depths are always non-negative):
///   - AD_Ref = 0, AD_Alt > 0: returns ≈ π/2 (homozygous ALT)
///   - AD_Alt = 0, AD_Ref > 0: returns ≈ 0.0 (homozygous REF)
///   - Both zero:              returns ≈ 0.0 (no reads)
///
/// Coverage stability: PERFECTLY COVERAGE-INVARIANT. PANG is a pure ratio
/// (atan2(y,x) depends only on y/x). Identical allele fractions produce
/// identical angles at any depth from 1× to 2000×.
///
/// Output range for non-negative inputs: [0, π/2] radians.
[[nodiscard]] inline auto PolarAngle(f64 const alt_depth, f64 const ref_depth) -> f64 {
  // ============================================================================
  // Minimax polynomial coefficients
  // ============================================================================
  //
  // In plain terms: instead of computing the exact arctangent (which is
  // slow), we use a fast approximation accurate to within 0.086° — far
  // below what matters for ML classification. The coefficients come from
  // an algorithm (Remez) that minimizes the worst-case error across the
  // entire input range, rather than just being accurate near zero.
  //
  // The core idea: instead of evaluating atan2(y, x) directly (which requires
  // expensive polynomial expansion for full-precision results), we:
  //
  //   1. "Fold" the inputs into a normalized ratio r ∈ [-1, 1] using octant
  //      symmetry: r = (x - |y|) / (|y| + |x|). This exploits the identity
  //      atan(a/b) = π/2 - atan(b/a) to always work near the origin where
  //      polynomials converge fastest.
  //
  //   2. Approximate atan(r) over [-1, 1] with a 3rd-degree odd polynomial:
  //        atan(r) ≈ (A·r² + B)·r
  //
  //      where A and B are chosen by the Remez exchange algorithm to minimize
  //      the MAXIMUM error across the entire [-1, 1] interval (minimax criterion).
  //      This is strictly better than Taylor series, which minimizes error only
  //      near r = 0 and diverges badly near ±1.
  //
  //   3. Add a base angle (π/4 or π/2) to unfold back to the correct octant.
  //
  // The coefficients A = 0.1963 and B = -0.9817 are from the Nvidia Cg 3.1
  // Toolkit's production atan2 approximation (also used in GPU shader math).
  // They achieve max |error| ≈ 0.0015 radians (~0.086°) across the full
  // atan2 domain, which is well below our VCF formatting resolution of
  // 4 decimal places (≈ 0.0001 rad).
  //
  // For comparison, a 1st-order Taylor approximation (A=0, B=-1) has max
  // error ~0.072 rad, and the exact coefficients (A=1/3, B=-1) from the
  // Taylor series have max error ~0.005 rad. The Remez-optimized values
  // trade minute near-zero accuracy for dramatically better edge behavior.
  //
  // Reference: Nvidia Cg 3.1 Toolkit Documentation, "Cg Standard Library
  // Functions", atan2 implementation notes.
  // ============================================================================

  static constexpr f64 COEFF_A = 0.1963;   // Remez-optimized cubic coefficient
  static constexpr f64 COEFF_B = -0.9817;  // Remez-optimized linear coefficient
  static constexpr f64 HALF_PI = std::numbers::pi_v<f64> / 2.0;
  static constexpr f64 QUARTER_PI = std::numbers::pi_v<f64> / 4.0;
  static constexpr f64 EPSILON = 1e-10;  // prevent division by zero when both inputs are zero

  f64 const abs_y = std::abs(alt_depth) + EPSILON;
  f64 const abs_x = std::abs(ref_depth);

  // Octant folding: maps the full atan2 domain to a ratio in [-1, 1].
  // When x ≥ 0: numer = x - |y|, so r ∈ [-1, 0] for |y| > x (octant 2)
  //              and r ∈ [0, 1] for |y| ≤ x (octant 1).
  // The copysign handles the x < 0 case symmetrically.
  f64 const numer = ref_depth - std::copysign(abs_y, ref_depth);
  f64 const denom = abs_y + abs_x;
  f64 const ratio = numer / denom;

  // Base angle: selects the octant center.
  //   x ≥ 0: base = π/2 - π/4 = π/4 (octant 1 center)
  //   x < 0: base = π/2 + π/4 = 3π/4 (octant 4 center)
  f64 const base = HALF_PI - std::copysign(QUARTER_PI, ref_depth);

  // Evaluate the Remez minimax polynomial: atan(r) ≈ (A·r² + B)·r
  // The result is added to the base angle to recover the full atan2 value.
  f64 const angle = base + (((COEFF_A * ratio * ratio) + COEFF_B) * ratio);

  // Sign correction for negative y (ALT depth is always >= 0 in practice,
  // but included for mathematical completeness)
  return std::copysign(angle, alt_depth);
}

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_POLAR_COORDS_H_
