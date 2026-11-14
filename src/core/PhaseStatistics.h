#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace PhaseStatistics {
inline double separation(double a, double b) { return std::abs(std::remainder(a - b, 360.)); }
// Circular medoid rejects a reversed minority without breaking at +/-180 degrees.
inline bool consensus(const std::vector<double> &phases, double tolerance, double &phase, unsigned &count) {
  if(phases.size() < 3) return false;
  const auto centre = std::min_element(phases.begin(), phases.end(), [&](double a, double b) {
    double da = 0., db = 0.;
    for(double value : phases) { da += separation(a, value); db += separation(b, value); }
    return da < db;
  });
  const double radians = std::acos(-1.) / 180.;
  std::complex<double> sum{}; count = 0;
  for(double value : phases) if(separation(value, *centre) <= tolerance) {
    sum += std::polar(1., value * radians); ++count;
  }
  if(count < 3 || count < .75 * phases.size()) return false;
  phase = std::arg(sum) / radians;
  return true;
}
}
