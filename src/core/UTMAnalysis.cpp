#include "core/UTMAnalysis.h"
#include <proj.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace UTM {
double longitude(double value) {
  return std::remainder(value, 360.);
}
double Extent::unwrap(double value) const {
  return longitude() + std::remainder(value - longitude(), 360.);
}
int zoneAt(double latitude, double lon) {
  if (!std::isfinite(latitude) || !std::isfinite(lon) || latitude < -80. || latitude > 84.)
    throw std::runtime_error("UTM requires latitudes between 80°S and 84°N.");
  lon = longitude(lon);
  if (latitude >= 56. && latitude < 64. && lon >= 3. && lon < 12.)
    return 32;
  if (latitude >= 72. && lon >= 0. && lon < 42.) {
    if (lon < 9.)
      return 31;
    if (lon < 21.)
      return 33;
    if (lon < 33.)
      return 35;
    return 37;
  }
  return std::max(1, std::min(60, int(std::floor((lon + 180.) / 6.)) + 1));
}
Extent extent(const std::vector<std::array<double, 3>> &points) {
  if (points.empty())
    throw std::runtime_error("Load stations before choosing UTM coordinates.");
  Extent result{0., 0., 90., -90.};
  std::vector<double> longitudes;
  for (const auto &p : points) {
    if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) || p[0] < -80. || p[0] > 84. ||
        p[1] < -180. || p[1] > 180.)
      throw std::runtime_error("UTM needs finite station coordinates between 80°S and 84°N.");
    longitudes.push_back(p[1] < 0. ? p[1] + 360. : p[1]);
    result.south = std::min(result.south, p[0]);
    result.north = std::max(result.north, p[0]);
  }
  std::sort(longitudes.begin(), longitudes.end());
  // Remove the largest empty arc, rather than averaging station density.
  double gap = -1.;
  size_t after = 0;
  for (size_t i = 0; i < longitudes.size(); ++i) {
    const double next = i + 1 < longitudes.size() ? longitudes[i + 1] : longitudes.front() + 360.;
    if (next - longitudes[i] > gap) {
      gap = next - longitudes[i];
      after = (i + 1) % longitudes.size();
    }
  }
  result.west = longitude(longitudes[after]);
  result.east = result.west + 360. - gap;
  return result;
}
std::vector<Region> regions(const Extent &view) {
  std::vector<Region> result;
  const double bands[]{-80., 0., 56., 64., 72., 84.};
  for (int b = 0; b < 5; ++b) {
    const double south = std::max(view.south, bands[b]), north = std::min(view.north, bands[b + 1]);
    if (north <= south)
      continue;
    std::vector<double> cuts;
    for (int z = 0; z <= 60; ++z)
      cuts.push_back(-180. + z * 6.);
    if (b == 2) {
      cuts.erase(std::find(cuts.begin(), cuts.end(), 6.));
      cuts.push_back(3.);
    }
    if (b == 4) {
      for (double v : {6., 12., 18., 24., 30., 36.})
        cuts.erase(std::find(cuts.begin(), cuts.end(), v));
      for (double v : {9., 21., 33.})
        cuts.push_back(v);
    }
    std::sort(cuts.begin(), cuts.end());
    const int first = int(std::floor((view.west + 180.) / 360.)),
              last = int(std::floor((view.east + 180.) / 360.));
    for (int world = first; world <= last; ++world)
      for (size_t i = 0; i + 1 < cuts.size(); ++i) {
        const double west = std::max(view.west, cuts[i] + world * 360.),
                     east = std::min(view.east, cuts[i + 1] + world * 360.);
        if (east > west)
          result.push_back(
              {zoneAt((south + north) / 2., (west + east) / 2.), south >= 0., west, east, south, north});
      }
  }
  return result;
}
Analysis analyze(const std::vector<std::array<double, 3>> &points, const SurveyCoordinates &coordinates) {
  const auto bounds = extent(points);
  if (coordinates.zone < 1 || coordinates.zone > 60)
    throw std::runtime_error("Invalid UTM zone.");
  std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)> context(proj_context_create(),
                                                                       proj_context_destroy);
  proj_log_level(context.get(), PJ_LOG_NONE);
  const auto definition = "+proj=utm +zone=" + std::to_string(coordinates.zone) + " +ellps=WGS84" +
                          (coordinates.north ? "" : " +south");
  std::unique_ptr<PJ, decltype(&proj_destroy)> projection(proj_create(context.get(), definition.c_str()),
                                                          proj_destroy);
  if (!projection)
    throw std::runtime_error("Cannot initialize UTM distortion analysis.");
  Analysis result;
  auto factor = [&](double lat, double lon) {
    proj_errno_reset(projection.get());
    const auto f =
        proj_factors(projection.get(), proj_coord(proj_torad(longitude(lon)), proj_torad(lat), 0., 0.));
    Factor out;
    if (proj_errno(projection.get()) || !std::isfinite(f.meridional_scale) ||
        !std::isfinite(f.parallel_scale) || !std::isfinite(f.meridian_convergence) ||
        f.meridional_scale <= 0. || f.parallel_scale <= 0.) {
      ++result.invalidSamples;
      return out;
    }
    out.scale = (f.meridional_scale + f.parallel_scale) / 2.;
    out.convergence = proj_todeg(f.meridian_convergence);
    result.minimumScale = std::min(result.minimumScale, out.scale);
    result.maximumScale = std::max(result.maximumScale, out.scale);
    result.extentMaxError = std::max(
        result.extentMaxError, std::max(std::abs(f.meridional_scale - 1.), std::abs(f.parallel_scale - 1.)));
    return out;
  };
  for (const auto &p : points) {
    const int z = zoneAt(p[0], p[1]);
    ++result.stationZones[z];
    if (z != coordinates.zone)
      ++result.outsideZone;
    if ((p[0] >= 0.) != coordinates.north)
      ++result.oppositeHemisphere;
    const auto f = factor(p[0], p[1]);
    result.stations.push_back(f);
    if (std::isfinite(f.scale)) {
      result.stationMaxError = std::max(result.stationMaxError, std::abs(f.scale - 1.));
      result.maxConvergence = std::max(result.maxConvergence, std::abs(f.convergence));
    }
  }
  for (int row = 0; row <= 32; ++row)
    for (int col = 0; col <= 32; ++col)
      factor(bounds.south + (bounds.north - bounds.south) * row / 32.,
             bounds.west + (bounds.east - bounds.west) * col / 32.);
  return result;
}
} // namespace UTM
