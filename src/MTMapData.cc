#include "include/MTMapData.h"
#include <algorithm>
#include <cmath>

namespace MTMapData {
int nearest_period(const dvector &frequencies, double period, double tolerance)
{
  if(!std::isfinite(period) || period <= 0. || !std::isfinite(tolerance) || tolerance < 0.) return -1;
  int best = -1;
  double difference = std::numeric_limits<double>::infinity();
  for(unsigned i = 0; i < frequencies.size(); ++i) {
    if(!std::isfinite(frequencies[i]) || frequencies[i] <= 0.) continue;
    const double relative = std::abs(1. / frequencies[i] - period) / period;
    if(relative <= tolerance && relative < difference) { best = i; difference = relative; }
  }
  return best;
}

bool phase_tensor(const std::array<double, 4> &tensor, PhaseTensor &result)
{
  double scale = 0.;
  for(double value: tensor) {
    if(!std::isfinite(value)) return false;
    scale = std::max(scale, std::abs(value));
  }
  if(scale == 0.) return false;
  const double a = tensor[0] / scale, b = tensor[1] / scale;
  const double c = tensor[2] / scale, d = tensor[3] / scale;
  // Caldwell et al. (2004), Bibby et al. (2005). Four-quadrant angles
  // preserve the orientation of non-symmetric and negative-determinant tensors.
  const double pi1 = .5 * std::hypot(a - d, b + c);
  const double pi2 = .5 * std::hypot(a + d, b - c);
  const double major = pi2 + pi1, minor = pi2 - pi1;
  if(major == 0.) return false;
  const double degrees = 180. / std::acos(-1.);
  const double alpha = .5 * std::atan2(b + c, a - d);
  const double beta = .5 * std::atan2(b - c, a + d);
  result.phiMin = std::atan(minor * scale) * degrees;
  result.phiMax = std::atan(major * scale) * degrees;
  result.skew = beta * degrees;
  result.azimuth = std::fmod((alpha - beta) * degrees + 360., 180.);
  result.axisRatio = std::abs(minor / major);
  return true;
}

bool phase_tensor(const MTStationData &station, unsigned frequency, PhaseTensor &result)
{
  if(frequency >= station.frequencies().size()) return false;
  const auto impedanceMask = station.impedance_mask();
  std::array<double, 4> tensor, real;
  const std::array<RealDataType, 4> types{{PTxx, PTxy, PTyx, PTyy}};
  const std::array<RealDataType, 4> realTypes{{RealZxx, RealZxy, RealZyx, RealZyy}};
  double error;
  bool hasRealImpedance = true;
  for(unsigned c = 0; c < 4; ++c) {
    if(!impedanceMask[c][frequency] || !station.scalar_value(types[c], frequency, tensor[c], error)) return false;
    hasRealImpedance &= station.scalar_value(realTypes[c], frequency, real[c], error);
  }
  if(hasRealImpedance) {
    double scale = 0.;
    for(double value: real) scale = std::max(scale, std::abs(value));
    if(scale > 0.) {
      for(double &value: real) value /= scale;
      if(std::abs(real[0] * real[3] - real[1] * real[2]) < 1e-12) return false;
    }
  }
  return phase_tensor(tensor, result);
}

bool induction_vector(const MTStationData &station, unsigned frequency, bool imaginary,
                      std::array<double, 2> &northEast)
{
  double error;
  return station.scalar_value(imaginary ? ImagTzx : RealTzx, frequency, northEast[0], error) &&
         station.scalar_value(imaginary ? ImagTzy : RealTzy, frequency, northEast[1], error);
}
}
