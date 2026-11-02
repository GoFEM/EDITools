#include "include/PolarityAnalysis.h"
#include "include/MTSurveyData.h"
#include "include/SurveyGeometry.h"
#include <algorithm>
#include <stdexcept>

namespace PolarityAnalysis {
namespace {
constexpr double pi = 3.14159265358979323846;
struct Sample { double frequency, phase; };
using Curve = std::vector<Sample>;
Curve samples(const MTStationData &station, unsigned component) {
  const RealDataType real[] = {RealZxy, RealZyx}, imag[] = {ImagZxy, ImagZyx}, phase[] = {PhsZxy, PhsZyx};
  Curve result;
  for(unsigned i = 0; i < station.frequencies().size(); ++i) {
    const double f = station.frequencies()[i];
    if(!std::isfinite(f) || f <= 0.) continue;
    double x, y, ex, ey, value, error;
    if(station.scalar_value(real[component], i, x, ex) && station.scalar_value(imag[component], i, y, ey)) {
      error = std::max(ex, ey);
      if(!std::isfinite(error) || error <= 0. || std::hypot(x, y) < 3. * error) continue;
      value = std::atan2(y, x) * 180. / pi;
    } else {
      if(!station.scalar_value(phase[component], i, value, error) || !std::isfinite(error) ||
         error <= 0. || error > 180. / (3. * pi)) continue;
    }
    result.push_back({f, std::remainder(value, 360.)});
  }
  std::sort(result.begin(), result.end(), [](const Sample &a, const Sample &b) { return a.frequency < b.frequency; });
  return result;
}
// Ambiguous near-duplicate samples are excluded instead of counted twice.
const Sample *match(const Curve &curve, double frequency) {
  const auto first = std::lower_bound(curve.begin(), curve.end(), frequency / 1.001,
    [](const Sample &a, double f) { return a.frequency < f; });
  if(first == curve.end() || first->frequency > frequency * 1.001) return nullptr;
  const auto next = std::next(first);
  return next != curve.end() && next->frequency <= frequency * 1.001 ? nullptr : &*first;
}
double separation(double a, double b) { return std::abs(std::remainder(a - b, 360.)); }
bool consensus(const std::vector<double> &phases, double tolerance, double &phase, unsigned &count) {
  if(phases.size() < 3) return false;
  // Circular medoid is robust to one reversed neighbour and the +/-180° seam.
  auto centre = std::min_element(phases.begin(), phases.end(), [&](double a, double b) {
    double da = 0., db = 0.;
    for(double value : phases) { da += separation(a, value); db += separation(b, value); }
    return da < db;
  });
  std::complex<double> sum{};
  count = 0;
  for(double value : phases) if(separation(value, *centre) <= tolerance) {
    sum += std::polar(1., value * pi / 180.); ++count;
  }
  if(count < 3 || count < .75 * phases.size()) return false;
  phase = std::arg(sum) * 180. / pi;
  return true;
}
}
bool Component::supported() const { return points.size() >= 3 && points.size() >= .6 * usable; }
bool Component::candidate() const { return supported() && reversed >= .8 * points.size(); }

std::vector<Station> analyse(const MTSurveyData &survey, const Options &options) {
  if(options.neighbours < 3 || !std::isfinite(options.radiusKm) || options.radiusKm <= 0. ||
     !std::isfinite(options.toleranceDegrees) || options.toleranceDegrees <= 0. || options.toleranceDegrees > 60.)
    throw std::invalid_argument("Invalid polarity comparison settings.");
  const auto names = survey.get_stations_names();
  std::vector<std::array<Curve, 2>> curves;
  for(const auto &name : names) {
    const auto &station = survey.get_station_data(name);
    curves.push_back({{samples(station, 0), samples(station, 1)}});
  }
  std::vector<Station> results;
  for(unsigned i = 0; i < names.size(); ++i) {
    const auto &station = survey.get_station_data(names[i]);
    std::vector<std::pair<double, unsigned>> nearby;
    if(station.active()) for(unsigned j = 0; j < names.size(); ++j) {
      if(i == j || !survey.is_active(names[j])) continue;
      const double km = SurveyGeometry::distance(station.position(), survey.get_station_data(names[j]).position());
      if(km <= options.radiusKm) nearby.emplace_back(km, j);
    }
    std::sort(nearby.begin(), nearby.end());
    if(nearby.size() > options.neighbours) nearby.resize(options.neighbours);
    Station result; result.name = names[i];
    for(unsigned c = 0; c < 2; ++c) {
      auto &component = result.components[c];
      component.usable = curves[i][c].size();
      for(const auto &sample : curves[i][c]) {
        if(match(curves[i][c], sample.frequency) != &sample) continue;
        std::vector<double> phases;
        for(const auto &neighbour : nearby)
          if(const auto *other = match(curves[neighbour.second][c], sample.frequency)) phases.push_back(other->phase);
        double reference; unsigned count;
        if(!consensus(phases, options.toleranceDegrees, reference, count)) continue;
        component.points.push_back({1. / sample.frequency, sample.phase, reference, count});
        if(separation(sample.phase + 180., reference) <= options.toleranceDegrees) ++component.reversed;
      }
    }
    results.push_back(std::move(result));
  }
  return results;
}
}
