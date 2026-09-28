#ifndef EDITOOLS_UTM_ANALYSIS_H
#define EDITOOLS_UTM_ANALYSIS_H

#include "SurveyCoordinates.h"
#include <limits>
#include <map>

namespace UTM {
struct Extent {
  // Longitude follows the shortest enclosing arc, possibly extending past 180°.
  double west, east, south, north;
  double longitude() const { return (west + east) / 2.; }
  double latitude() const { return (south + north) / 2.; }
  double unwrap(double longitude) const;
};
struct Region {
  int zone;
  bool north;
  double west, east, south, top;
};
struct Factor {
  double scale = std::numeric_limits<double>::quiet_NaN();
  double convergence = std::numeric_limits<double>::quiet_NaN(); // degrees
};
struct Analysis {
  std::vector<Factor> stations;
  std::map<int, unsigned> stationZones;
  unsigned outsideZone = 0, oppositeHemisphere = 0, invalidSamples = 0;
  double stationMaxError = 0., extentMaxError = 0.; // fractional local scale error
  double minimumScale = std::numeric_limits<double>::infinity(), maximumScale = 0.;
  double maxConvergence = 0.; // absolute degrees at stations
  bool valid() const { return invalidSamples == 0; }
};
double longitude(double value);
int zoneAt(double latitude, double longitude);
Extent extent(const std::vector<std::array<double, 3>> &geographic);
// Actual zone boundaries, including Norway/Svalbard and hemisphere changes.
std::vector<Region> regions(const Extent &view);
// Local ellipsoidal scale and convergence at stations and a 33 × 33 extent grid.
// The sampled extent maximum is an estimate, not a rigorous distortion bound.
Analysis analyze(const std::vector<std::array<double, 3>> &geographic, const SurveyCoordinates &coordinates);
} // namespace UTM
#endif
