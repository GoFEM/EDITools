#pragma once
#include <array>
#include <cmath>
#include <geodesic.h>

namespace SurveyGeometry {
inline bool located(const std::array<double, 3> &p) {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::abs(p[0]) <= 90. && std::abs(p[1]) <= 180.;
}
inline double distance(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  if (!located(a) || !located(b)) return INFINITY;
  geod_geodesic earth;
  geod_init(&earth, 6378137., 1. / 298.257223563);
  double metres;
  geod_inverse(&earth, a[0], a[1], b[0], b[1], &metres, nullptr, nullptr);
  return metres / 1000.;
}
}
