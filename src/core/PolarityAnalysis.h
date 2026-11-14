#pragma once
#include <array>
#include <string>
#include <vector>

class MTSurveyData;
namespace PolarityAnalysis {
struct Options {
  unsigned neighbours = 5;
  double radiusKm = 100.;
  double toleranceDegrees = 30.;
};
struct Point {
  double period, phase, referencePhase;
  unsigned neighbours;
};
struct Component {
  unsigned usable = 0, reversed = 0;
  std::vector<Point> points;
  bool supported() const;
  bool candidate() const;
};
struct Station {
  std::string name;
  std::array<Component, 2> components; // XY, YX; raw phase convention, no 180° folding.
};
// Review candidates, never proof of a wiring fault. Uses enabled, unmasked
// samples with SNR >= 3 (or equivalent phase error), period matches <= 0.1%,
// and at least three neighbours with a >=75% circular phase consensus.
std::vector<Station> analyse(const MTSurveyData &survey, const Options &options = {});
}
