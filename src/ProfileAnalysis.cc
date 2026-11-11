#include "include/ProfileAnalysis.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace ProfileAnalysis {
namespace {
double distance(Point a, Point b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }
bool finite(Point p) { return std::isfinite(p[0]) && std::isfinite(p[1]); }
double median(std::vector<double> values) {
  if (values.empty())
    return 0.;
  std::sort(values.begin(), values.end());
  return .5 * (values[(values.size() - 1) / 2] + values[values.size() / 2]);
}
} // namespace
double length(const std::vector<Point> &line) {
  double result = 0.;
  for (unsigned i = 0; i < line.size(); ++i) {
    if (!finite(line[i]))
      throw std::invalid_argument("Profile vertices must be finite.");
    if (i)
      result += distance(line[i - 1], line[i]);
  }
  return result;
}
double medianSpacing(const std::vector<Station> &stations) {
  std::vector<double> nearest;
  for (unsigned i = 0; i < stations.size(); ++i) {
    double best = std::numeric_limits<double>::infinity();
    if (!finite(stations[i].position))
      continue;
    for (unsigned j = 0; j < stations.size(); ++j)
      if (i != j && finite(stations[j].position)) {
        const double d = distance(stations[i].position, stations[j].position);
        if (d > 1e-9)
          best = std::min(best, d);
      }
    if (std::isfinite(best))
      nearest.push_back(best);
  }
  return median(nearest);
}
std::vector<Offset> select(const std::vector<Station> &stations, const std::vector<Point> &line,
                           double halfWidth) {
  if (!std::isfinite(halfWidth) || halfWidth < 0.)
    throw std::invalid_argument("Invalid corridor distance.");
  std::vector<Offset> result;
  if (length(line) <= 1e-9)
    return result;
  for (const auto &site : stations) {
    if (!finite(site.position))
      continue;
    double best = std::numeric_limits<double>::infinity(), travelled = 0.;
    Offset chosen{site.name, 0., 0., {}};
    for (unsigned i = 1; i < line.size(); ++i) {
      const auto &a = line[i - 1], &b = line[i];
      const double dx = b[0] - a[0], dy = b[1] - a[1], segment = distance(a, b);
      if (segment <= 1e-9)
        continue;
      const double t =
          std::max(0., std::min(1., ((site.position[0] - a[0]) * dx + (site.position[1] - a[1]) * dy) /
                                        (segment * segment)));
      const Point foot{{a[0] + t * dx, a[1] + t * dy}};
      const double d = distance(site.position, foot);
      // Equal-distance crossings choose the first segment, so each station occurs once.
      if (d < best) {
        best = d;
        chosen.along = travelled + t * segment;
        chosen.foot = foot;
        chosen.across =
            std::copysign(d, dx * (site.position[1] - foot[1]) - dy * (site.position[0] - foot[0]));
      }
      travelled += segment;
    }
    if (best <= halfWidth + 1e-9)
      result.push_back(chosen);
  }
  std::sort(result.begin(), result.end(), [](const Offset &a, const Offset &b) {
    return a.along == b.along ? a.name < b.name : a.along < b.along;
  });
  return result;
}
std::vector<std::array<double, 2>> cells(const std::vector<double> &sorted, double maximumWidth) {
  if (!std::isfinite(maximumWidth) || maximumWidth <= 0.)
    throw std::invalid_argument("Invalid cell width.");
  std::vector<std::array<double, 2>> result(sorted.size());
  for (unsigned i = 0; i < sorted.size(); ++i)
    if (!std::isfinite(sorted[i]) || (i && sorted[i] < sorted[i - 1]))
      throw std::invalid_argument("Cell positions must be finite and sorted.");
  for (unsigned first = 0; first < sorted.size();) {
    unsigned next = first + 1;
    while (next < sorted.size() && sorted[next] == sorted[first])
      ++next;
    const double left = first ? std::min(maximumWidth, sorted[first] - sorted[first - 1]) : maximumWidth;
    const double right =
        next < sorted.size() ? std::min(maximumWidth, sorted[next] - sorted[first]) : maximumWidth;
    for (unsigned i = first; i < next; ++i)
      result[i] = {{sorted[first] - .5 * left, sorted[first] + .5 * right}};
    first = next;
  }
  return result;
}
Rose rose(const std::vector<Direction> &directions, unsigned binDegrees, bool equalStations) {
  if (!binDegrees || 180 % binDegrees)
    throw std::invalid_argument("Rose bin width must divide 180 degrees.");
  Rose result;
  result.percent.resize(180 / binDegrees);
  result.mean = NAN;
  std::map<std::string, unsigned> counts;
  for (const auto &d : directions)
    if (std::isfinite(d.degrees)) {
      ++counts[d.station];
      ++result.samples;
    }
  result.stations = counts.size();
  double total = 0., x = 0., y = 0.;
  const double radians = std::acos(-1.) / 180.;
  for (const auto &d : directions)
    if (std::isfinite(d.degrees)) {
      const double a = std::fmod(std::fmod(d.degrees, 180.) + 180., 180.);
      const double weight = equalStations ? 1. / counts[d.station] : 1.;
      result.percent[std::min(unsigned(a / binDegrees), unsigned(result.percent.size() - 1))] += weight;
      x += weight * std::cos(2 * a * radians);
      y += weight * std::sin(2 * a * radians);
      total += weight;
    }
  if (total > 0.) {
    for (auto &value : result.percent)
      value *= 100. / total;
    result.concentration = std::hypot(x, y) / total;
    if (result.concentration > 1e-12)
      result.mean = std::fmod(.5 * std::atan2(y, x) / radians + 180., 180.);
  }
  return result;
}
Strike estimateStrike(const std::vector<Direction> &directions, bool equalStations) {
  std::vector<Direction> doubled;
  doubled.reserve(directions.size());
  for (const auto &direction : directions)
    doubled.push_back({direction.station, 2. * direction.degrees});
  const auto axes = rose(doubled, 10, equalStations);
  return {axes.mean * .5, axes.concentration, axes.samples, axes.stations};
}
} // namespace ProfileAnalysis
