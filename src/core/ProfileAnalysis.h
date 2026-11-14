#ifndef EDITOOLS_PROFILE_ANALYSIS_H
#define EDITOOLS_PROFILE_ANALYSIS_H

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace ProfileAnalysis {
using Point = std::array<double, 2>; // East, north, in one common projection (km).
struct Station {
  std::string name;
  Point position;
};
struct Offset {
  std::string name;
  double along, across; // Cumulative distance and signed perpendicular/end-cap distance (km).
  Point foot;
};
double length(const std::vector<Point> &line);
double
medianSpacing(const std::vector<Station> &stations); // Nonzero nearest-neighbour distances; 0 if undefined.
std::vector<Offset> select(const std::vector<Station> &stations, const std::vector<Point> &line,
                           double halfWidth);
// Bounded cells around actual sample positions. Large gaps remain blank.
std::vector<std::array<double, 2>> cells(const std::vector<double> &sorted, double maximumWidth);

struct Direction {
  std::string station;
  double degrees;
};
struct Rose {
  std::vector<double> percent; // Axial bins in [0,180); sum is 100 for nonempty input.
  unsigned samples = 0, stations = 0;
  double mean = 0., concentration = 0.; // Axial circular mean; NaN when undefined.
};
Rose rose(const std::vector<Direction> &directions, unsigned binDegrees, bool equalStations);
struct Strike {
  double degrees = std::numeric_limits<double>::quiet_NaN(), concentration = 0.;
  unsigned samples = 0, stations = 0;
};
// Principal-axis mean modulo 90 degrees, using the fourth circular moment.
// The perpendicular branch remains ambiguous; no data are rotated here.
Strike estimateStrike(const std::vector<Direction> &directions, bool equalStations);
} // namespace ProfileAnalysis
#endif
