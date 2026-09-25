#include "include/MTSpectra.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace MTSpectra {
namespace {
using Complex = std::complex<double>;
bool finite(Complex value) { return std::isfinite(value.real()) && std::isfinite(value.imag()); }
}

Result estimate(const PackedMatrix &packed, double averages)
{
  if(!std::isfinite(averages) || averages <= 0.)
    throw std::invalid_argument("Spectral AVGT must be finite and positive.");
  const double missing = std::numeric_limits<double>::quiet_NaN();
  Result result;
  result.impedance.fill({missing, missing}); result.tipper.fill({missing, missing});
  result.impedanceError.fill(missing); result.tipperError.fill(missing);

  // A common scale cancels from both the solve and its uncertainty. Use the
  // predictor/reference block so a missing output cannot invalidate good data.
  double scale = 0.;
  for(unsigned i: {0u, 1u, 5u, 6u}) for(unsigned j: {0u, 1u, 5u, 6u}) {
    if(!std::isfinite(packed[i][j])) throw std::invalid_argument("Missing or non-finite spectral reference data.");
    scale = std::max(scale, std::abs(packed[i][j]));
  }
  if(scale == 0.) throw std::invalid_argument("Spectral reference power must be positive.");
  std::array<std::array<Complex, 7>, 7> power;
  for(unsigned i = 0; i < 7; ++i) {
    power[i][i] = {packed[i][i] / scale, 0.};
    for(unsigned j = i + 1; j < 7; ++j) {
      power[i][j] = {packed[j][i] / scale, -packed[i][j] / scale};
      power[j][i] = std::conj(power[i][j]);
    }
  }
  for(unsigned i: {0u, 1u, 5u, 6u})
    if(power[i][i].real() <= 0.) throw std::invalid_argument("Spectral reference power must be positive.");

  // Solve [Zx Zy] S(H,R) = S(output,R) for Ex, Ey and Hz.
  const auto a = power[0][5], b = power[0][6], c = power[1][5], d = power[1][6];
  const auto determinant = a * d - b * c;
  const double size = std::max({std::abs(a), std::abs(b), std::abs(c), std::abs(d)});
  if(size == 0. || std::abs(determinant) <= 64. * std::numeric_limits<double>::epsilon() * size * size)
    throw std::invalid_argument("Singular horizontal/remote spectral matrix.");

  // Diagonal of inv(S(H,R))^H S(R,R) inv(S(H,R)), divided by AVGT.
  const double rx = power[5][5].real(), ry = power[6][6].real();
  const auto cross = power[5][6];
  const double denominator = std::norm(determinant);
  const std::array<double, 2> weights{{
    std::abs(rx * std::norm(d) + ry * std::norm(c) - 2. * std::real(cross * std::conj(d) * c)) / denominator / averages,
    std::abs(ry * std::norm(a) + rx * std::norm(b) - 2. * std::real(cross * std::conj(b) * a)) / denominator / averages
  }};
  for(unsigned output: {3u, 4u, 2u}) {
    if(!finite(power[output][5]) || !finite(power[output][6])) continue;
    const std::array<Complex, 2> transfer{{
      (power[output][5] * d - c * power[output][6]) / determinant,
      (a * power[output][6] - power[output][5] * b) / determinant
    }};
    // Power of output - Zx*Hx - Zy*Hy. Absolute value preserves the established
    // variance convention for spectra with small negative residual roundoff.
    double residual = missing;
    if(finite(power[output][output]) && power[output][output].real() >= 0. &&
       finite(power[0][output]) && finite(power[1][output])) {
      residual = std::abs(power[output][output].real()
        - 2. * std::real(transfer[0] * power[0][output] + transfer[1] * power[1][output])
        + std::norm(transfer[0]) * power[0][0].real() + std::norm(transfer[1]) * power[1][1].real()
        + 2. * std::real(transfer[0] * std::conj(transfer[1]) * power[0][1]));
    }
    for(unsigned component = 0; component < 2; ++component) {
      const unsigned index = output == 2 ? component : (output - 3) * 2 + component;
      const double error = std::sqrt(residual * weights[component]);
      if(output == 2) { result.tipper[index] = transfer[component]; result.tipperError[index] = error; }
      else { result.impedance[index] = transfer[component]; result.impedanceError[index] = error; }
    }
  }
  return result;
}
}
