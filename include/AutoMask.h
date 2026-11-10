#pragma once
#include "ImpedanceAnalysis.h"
#include <functional>

class MTSurveyData;
namespace AutoMask {
constexpr unsigned criterionCount = 4;
struct Options {
  std::array<bool, criterionCount> criteria{{true, true, true, true}}; // rho+/dispersion, period trend, spatial phase, uncertainty
  std::array<bool, 4> components{{true, true, true, true}}; // XY, YX, Tzx, Tzy
  std::array<double, criterionCount> threshold{{3., 3., 3., 1.}}; // sigma for residuals; limit multiplier for uncertainty
  unsigned minimumFlags = 1, neighbours = 5;
  double maximumGapRatio = 3., radiusKm = 100., periodTolerance = .001;
  double minimumPhaseDifference = 15., phaseAgreement = 30.;
  double maximumRelativeRhoError = .5, maximumPhaseError = 15.;
  double maximumTipperError = .05, minimumTipperDifference = .05, tipperAgreement = .1;
  ImpedanceAnalysis::RhoPlusOptions rhoPlus;
};
struct Point {
  std::string station;
  RealDataType component;
  double frequency;
  std::array<double, criterionCount> scores{{NAN, NAN, NAN, NAN}};
  std::array<std::string, criterionCount> details;
  std::array<bool, criterionCount> flagged{};
  unsigned votes() const;
  bool candidate = false;
};
struct Result {
  std::vector<Point> points;
  unsigned disabled = 0;
  bool cancelled = false;
};
// All enabled survey stations can provide spatial references, even when they
// are not targets. Missing criteria never count as votes. One pass, read-only.
Result analyse(const MTSurveyData &survey, const std::set<std::string> &targets,
               const Options &options = {}, const std::function<bool(unsigned, unsigned)> &progress = {});
const char *componentName(RealDataType component);
using Edits = std::map<std::string, std::vector<std::pair<RealDataType, double>>>;
// Only expands impedance detections to paired diagonals/PT. Tipper masks are independent.
Edits maskPlan(const MTSurveyData &survey, const std::vector<Point> &points,
               const std::vector<bool> &checked, bool pairDiagonals, bool linkTensorMasks);
void apply(MTSurveyData &survey, const Edits &edits);
}
