#ifndef EDITOOLS_MT_SPECTRA_H
#define EDITOOLS_MT_SPECTRA_H

#include <array>
#include <complex>

namespace MTSpectra {
// EDI packed power/cross-power matrix, in Hx, Hy, Hz, Ex, Ey, Rx, Ry order.
// Diagonal: real auto-power. Lower triangle: real cross-power. Upper
// triangle: imaginary cross-power, with S(i,j) = lower(j,i) - i*upper(i,j).
using PackedMatrix = std::array<std::array<double, 7>, 7>;
struct Result {
  // Row-major xx, xy, yx, yy; field units (mV/km)/nT, as in EDI MTSECT.
  std::array<std::complex<double>, 4> impedance;
  std::array<std::complex<double>, 2> tipper;
  // Standard deviations, in the same units as their transfer functions.
  std::array<double, 4> impedanceError;
  std::array<double, 2> tipperError;
};

// Remote-reference estimate and residual-power uncertainties. Averages must
// be positive. Invalid reference channels or singular systems throw; missing
// output channels produce NaN for that output without affecting other outputs.
Result estimate(const PackedMatrix &spectra, double averages);
}
#endif
