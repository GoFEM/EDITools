#include "core/AutoMask.h"
#include "core/MTSurveyData.h"
#include "core/PhaseStatistics.h"
#include "core/SurveyGeometry.h"
#include <algorithm>
#include <sstream>

namespace AutoMask {
namespace {
struct Sample {
  double period, frequency;
  std::array<double, 2> value{{NAN, NAN}}, errors{{NAN, NAN}};
};
using Curve = std::vector<Sample>;
bool error(double e) { return std::isfinite(e) && e > 0.; }
bool channelValid(const Sample &s, unsigned channel) { return std::isfinite(s.value[channel]) && error(s.errors[channel]); }
double median(std::vector<double> values) {
  std::sort(values.begin(), values.end()); const auto n = values.size();
  return n % 2 ? values[n/2] : .5*values[n/2-1] + .5*values[n/2];
}
double spread(const std::vector<double> &values) {
  const double centre = median(values); std::vector<double> residuals;
  for(double value : values) residuals.push_back(std::abs(value - centre));
  return 1.4826 * median(residuals);
}
std::string number(double value) { std::ostringstream out; out.precision(4); out << value; return out.str(); }
Curve samples(const MTStationData &station, unsigned c) {
  const RealDataType first[]{RhoZxy, RhoZyx, RealTzx, RealTzy}, second[]{PhsZxy, PhsZyx, ImagTzx, ImagTzy}; Curve result;
  for(unsigned i = 0; i < station.frequencies().size(); ++i) {
    const double f = station.frequencies()[i], p = 1. / f;
    if(!std::isfinite(f) || f <= 0. || !std::isfinite(p) || p <= 0.) continue;
    Sample s; s.frequency = f; s.period = p;
    double value, e;
    if(station.scalar_value(first[c], i, value, e) && (c >= 2 || value > 0.)) { s.value[0] = value; s.errors[0] = e; }
    if(station.scalar_value(second[c], i, value, e)) { s.value[1] = value; s.errors[1] = e; }
    result.push_back(s);
  }
  std::sort(result.begin(), result.end(), [](const Sample &a, const Sample &b) { return a.period < b.period; });
  return result;
}
const Sample *match(const Curve &curve, double period, double tolerance) {
  const auto first = std::lower_bound(curve.begin(), curve.end(), period/(1.+tolerance),
                                    [](const Sample &s, double p) { return s.period < p; });
  if(first == curve.end() || first->period > period*(1.+tolerance)) return nullptr;
  const auto next = std::next(first);
  return next != curve.end() && next->period <= period*(1.+tolerance) ? nullptr : &*first;
}
// Leave the tested point out of a robust local linear trend in log period.
// Four reference samples, two on each side, avoid endpoint extrapolation.
double trendScore(const Curve &curve, unsigned i, bool phase, bool tipper) {
  std::vector<double> x, y, errors;
  const auto &target = curve[i];
  const auto valid = [phase, tipper](const Sample &s) {
    if(tipper) return channelValid(s, phase);
    return phase ? channelValid(s, 1) : std::isfinite(s.value[0]) && s.value[0] > 0. && error(s.errors[0]) && error(s.errors[0]/s.value[0]);
  };
  if(!valid(target)) return NAN;
  for(unsigned j : {i-2, i-1, i+1, i+2}) {
    const auto &s = curve[j]; if(!valid(s)) return NAN;
    x.push_back(std::log(s.period/target.period));
    y.push_back(tipper ? s.value[phase] : phase ? s.value[1] : std::log(s.value[0]));
    errors.push_back(tipper ? s.errors[phase] : phase ? s.errors[1] : s.errors[0]/s.value[0]);
  }
  if(phase && !tipper) for(unsigned j = 1; j < y.size(); ++j) y[j] = y[j-1] + std::remainder(y[j]-y[j-1], 360.);
  std::vector<double> slopes, intercepts;
  for(unsigned a = 0; a < x.size(); ++a) for(unsigned b = a+1; b < x.size(); ++b)
    slopes.push_back((y[b]-y[a])/(x[b]-x[a]));
  const double slope = median(slopes);
  for(unsigned j = 0; j < x.size(); ++j) intercepts.push_back(y[j]-slope*x[j]);
  const double predicted = median(intercepts);
  const double residual = tipper ? std::abs(target.value[phase]-predicted) : phase ? PhaseStatistics::separation(target.value[1], predicted) : std::abs(std::log(target.value[0])-predicted);
  const double uncertainty = std::hypot(tipper ? target.errors[phase] : phase ? target.errors[1] : target.errors[0]/target.value[0],
                                      std::hypot(median(errors), spread(intercepts)));
  return residual / uncertainty;
}
void periodTest(Point &point, const Curve &curve, unsigned i, const Options &options, bool tipper) {
  point.details[1] = "Needs two usable periods on each side; no masked gaps or extrapolation.";
  if(i < 2 || i+2 >= curve.size()) return;
  for(unsigned j = i-2; j < i+2; ++j) {
    const double ratio = curve[j+1].period/curve[j].period;
    if(ratio <= 1.+1e-12 || ratio > options.maximumGapRatio) {
      point.details[1] = "Ambiguous periods or a bounding gap exceeds the limit."; return;
    }
  }
  const double rho = trendScore(curve, i, false, tipper), phase = trendScore(curve, i, true, tipper);
  if(std::isfinite(rho)) point.scores[1] = rho;
  if(std::isfinite(phase)) point.scores[1] = std::isfinite(point.scores[1]) ? std::max(point.scores[1], phase) : phase;
  if(std::isfinite(point.scores[1])) point.details[1] = std::string("Leave-one-out trend: ") + (tipper ? "real " : "log-rho ") + number(rho) + (tipper ? " sigma; imaginary " : " sigma; phase ") + number(phase) + " sigma. Four reference periods.";
}
// Poor precision is distinct from a large normalized residual: a noisy datum
// can agree with almost any model when its propagated uncertainty is large.
void uncertaintyTest(Point &point, const Sample &sample, const Options &options, bool tipper) {
  if(tipper) {
    for(unsigned k = 0; k < 2; ++k) if(std::isfinite(sample.value[k]) && error(sample.errors[k]))
      point.scores[3] = std::isfinite(point.scores[3]) ? std::max(point.scores[3], sample.errors[k]/options.maximumTipperError) : sample.errors[k]/options.maximumTipperError;
    point.details[3] = "Tipper errors: real " + number(sample.errors[0]) + "; imaginary " + number(sample.errors[1]) +
      ". Absolute limit " + number(options.maximumTipperError*options.threshold[3]) + ". No division by amplitude: zero tippers are valid.";
    return;
  }
  const double relative = std::isfinite(sample.value[0]) && sample.value[0] > 0. && error(sample.errors[0]) ? sample.errors[0]/sample.value[0] : NAN;
  const double phaseError = channelValid(sample, 1) ? sample.errors[1] : NAN;
  if(std::isfinite(relative)) point.scores[3] = relative/options.maximumRelativeRhoError;
  if(std::isfinite(phaseError)) point.scores[3] = std::isfinite(point.scores[3]) ?
    std::max(point.scores[3], phaseError/options.maximumPhaseError) : phaseError/options.maximumPhaseError;
  point.details[3] = "Main-data uncertainty: rho error " + (std::isfinite(relative) ? number(100.*relative) + "%" : "unavailable") +
    "; phase error " + (std::isfinite(phaseError) ? number(phaseError) + " degrees" : "unavailable") +
    ". Limits: " + number(100.*options.maximumRelativeRhoError*options.threshold[3]) + "% or " +
    number(options.maximumPhaseError*options.threshold[3]) + " degrees. Score is error / base limit, not sigma. Errors are not modified.";
}
}
unsigned Point::votes() const { return std::count(flagged.begin(), flagged.end(), true); }

Result analyse(const MTSurveyData &survey, const std::set<std::string> &targets,
               const Options &options, const std::function<bool(unsigned, unsigned)> &progress) {
  for(double threshold : options.threshold) if(!error(threshold)) throw std::invalid_argument("Thresholds must be finite and positive.");
  if(!std::isfinite(options.maximumGapRatio) || options.maximumGapRatio <= 1. || !error(options.radiusKm) ||
     !std::isfinite(options.periodTolerance) || options.periodTolerance < 0. || options.periodTolerance > .1 ||
     !error(options.minimumPhaseDifference) || options.minimumPhaseDifference > 180. ||
     !error(options.phaseAgreement) || options.phaseAgreement > 60. || options.neighbours < 3 ||
     !error(options.maximumTipperError) || !error(options.minimumTipperDifference) || !error(options.tipperAgreement) ||
     !error(options.maximumRelativeRhoError) || !error(options.maximumPhaseError) ||
     options.minimumFlags < 1 || options.minimumFlags > criterionCount)
    throw std::invalid_argument("Invalid auto-masking settings.");
  for(const auto &name : targets) if(!survey.is_station_present(name)) throw std::invalid_argument("Selected station is no longer present.");
  std::map<std::string, std::array<Curve, 4>> curves;
  for(const auto &name : survey.get_stations_names()) if(survey.is_active(name))
    for(unsigned c = 0; c < options.components.size(); ++c) if(options.components[c]) curves[name][c] = samples(survey.get_station_data(name), c);
  Result result; unsigned done = 0;
  for(const auto &name : targets) {
    if(progress && !progress(done++, targets.size())) { result.points.clear(); result.cancelled = true; return result; }
    const auto &station = survey.get_station_data(name);
    if(!station.active()) { ++result.disabled; continue; }
    std::vector<std::pair<double, std::string>> neighbours;
    if(options.criteria[2]) for(const auto &other : curves) if(other.first != name) {
      const double km = SurveyGeometry::distance(station.position(), survey.get_station_data(other.first).position());
      if(km <= options.radiusKm) neighbours.emplace_back(km, other.first);
    }
    std::sort(neighbours.begin(), neighbours.end());
    if(neighbours.size() > options.neighbours) neighbours.resize(options.neighbours);
    for(unsigned c = 0; c < options.components.size(); ++c) if(options.components[c]) {
      const auto &curve = curves.at(name)[c]; const bool tipper = c >= 2;
      std::map<double, double> rhoScores; std::string rhoStatus = "Disabled.";
      if(options.criteria[0] && tipper) rhoStatus = "Rho+ does not apply to tippers.";
      if(options.criteria[0] && !tipper) try {
        const auto fit = ImpedanceAnalysis::rhoPlus(ImpedanceAnalysis::curve(station, c == 0 ? ImpedanceAnalysis::Quantity::XY : ImpedanceAnalysis::Quantity::YX), options.rhoPlus);
        if(!fit.converged || !fit.gridStable) rhoStatus = "Provisional rho+ fit; excluded from automatic masking.";
        else {
          rhoStatus = "Joint log-rho/phase residual from a stable rho+ fit; shared dispersion criterion, not a separate Kramers-Kronig vote.";
          for(unsigned i = 0; i < fit.observed.size(); ++i)
            rhoScores[fit.observed[i].period] = std::max(std::abs(fit.rhoResidual[i]), std::abs(fit.phaseResidual[i]));
        }
      } catch(const std::exception &e) { rhoStatus = e.what(); }
      for(unsigned i = 0; i < curve.size(); ++i) {
        const auto &s = curve[i]; if(!std::isfinite(s.value[0]) && !std::isfinite(s.value[1])) continue;
        Point point; point.station = name; point.component = std::array<RealDataType, 4>{{RealZxy, RealZyx, RealTzx, RealTzy}}[c]; point.frequency = s.frequency;
        point.details.fill("Disabled."); point.details[0] = rhoStatus;
        if(rhoScores.count(s.period)) point.scores[0] = rhoScores.at(s.period);
        else if(options.criteria[0] && !rhoScores.empty()) point.details[0] = "Sample lacks usable complex impedance or positive propagated errors.";
        if(options.criteria[1]) periodTest(point, curve, i, options, tipper);
        bool spatialDeparture = false;
        if(options.criteria[2] && !tipper) {
          point.details[2] = "Needs three nearby, enabled stations with unique matching periods, positive phase errors and 75% phase agreement.";
          if(channelValid(s, 1) && match(curve, s.period, options.periodTolerance) == &s) {
            std::vector<double> phases, errors;
            for(const auto &neighbour : neighbours)
              if(const auto *other = match(curves.at(neighbour.second)[c], s.period, options.periodTolerance))
                if(channelValid(*other, 1)) { phases.push_back(other->value[1]); errors.push_back(other->errors[1]); }
            double reference; unsigned agreeing;
            if(PhaseStatistics::consensus(phases, options.phaseAgreement, reference, agreeing)) {
              std::vector<double> residuals;
              for(double phase : phases) residuals.push_back(std::remainder(phase-reference, 360.));
              const double departure = PhaseStatistics::separation(s.value[1], reference);
              point.scores[2] = departure / std::hypot(s.errors[1], std::hypot(median(errors), spread(residuals)));
              spatialDeparture = departure >= options.minimumPhaseDifference;
              point.details[2] = number(departure) + " degrees from " + std::to_string(agreeing) + "/" + std::to_string(phases.size()) +
                " agreeing neighbours; minimum departure " + number(options.minimumPhaseDifference) + " degrees. Phase only; no amplitude comparison.";
            }
          }
        }
        if(options.criteria[2] && tipper) {
          point.details[2] = "Needs three matched neighbours and 75% agreement of complex tipper values.";
          if(match(curve, s.period, options.periodTolerance) == &s &&
             channelValid(s, 0) && channelValid(s, 1)) {
            std::vector<const Sample *> refs;
            std::array<std::vector<double>, 2> values, errors;
            for(const auto &neighbour : neighbours)
              if(const auto *other = match(curves.at(neighbour.second)[c], s.period, options.periodTolerance))
                if(channelValid(*other, 0) && channelValid(*other, 1)) {
                  refs.push_back(other);
                  for(unsigned k = 0; k < 2; ++k) { values[k].push_back(other->value[k]); errors[k].push_back(other->errors[k]); }
                }
            if(refs.size() >= 3) {
              const std::array<double, 2> centre{{median(values[0]), median(values[1])}};
              const auto agreeing = std::count_if(refs.begin(), refs.end(), [&](const Sample *r) {
                return std::hypot(r->value[0]-centre[0], r->value[1]-centre[1]) <= options.tipperAgreement;
              });
              if(agreeing >= 3 && agreeing >= .75*refs.size()) {
                double score = 0.;
                for(unsigned k = 0; k < 2; ++k) score = std::max(score, std::abs(s.value[k]-centre[k]) /
                  std::hypot(s.errors[k], std::hypot(median(errors[k]), spread(values[k]))));
                const double departure = std::hypot(s.value[0]-centre[0], s.value[1]-centre[1]);
                point.scores[2] = score; spatialDeparture = departure >= options.minimumTipperDifference;
                point.details[2] = "Complex tipper departure " + number(departure) + " from " + std::to_string(agreeing) +
                  "/" + std::to_string(refs.size()) + " agreeing neighbours; minimum " + number(options.minimumTipperDifference) + ".";
              }
            }
          }
        }
        if(options.criteria[3]) uncertaintyTest(point, s, options, tipper);
        for(unsigned k = 0; k < criterionCount; ++k) point.flagged[k] = options.criteria[k] && std::isfinite(point.scores[k]) &&
          point.scores[k] > options.threshold[k] && (k != 2 || spatialDeparture);
        point.candidate = point.votes() >= options.minimumFlags;
        result.points.push_back(std::move(point));
      }
    }
  }
  if(progress && !progress(done, targets.size())) { result.points.clear(); result.cancelled = true; }
  return result;
}

const char *componentName(RealDataType component) {
  switch(component) { case RealZxy: return "XY"; case RealZyx: return "YX"; case RealTzx: return "Tzx"; case RealTzy: return "Tzy"; default: return "?"; }
}
Edits maskPlan(const MTSurveyData &survey, const std::vector<Point> &points,
               const std::vector<bool> &checked, bool pairDiagonals, bool linkTensorMasks) {
  if(checked.size() != points.size()) throw std::invalid_argument("Mask selection size mismatch.");
  std::map<std::string, std::set<std::pair<RealDataType, double>>> unique;
  for(unsigned i = 0; i < points.size(); ++i) if(checked[i]) {
    const auto &p = points[i]; auto &changes = unique[p.station];
    for(auto type: MTStationData::mask_edit_types(p.component, linkTensorMasks, pairDiagonals))
      changes.emplace(type, p.frequency);
  }
  Edits result;
  for(const auto &entry : unique) {
    const auto &station = survey.get_station_data(entry.first);
    const auto z = station.impedance_mask(true), pt = station.phase_tensor_mask(true), t = station.tipper_mask(true);
    for(const auto &change : entry.second) {
      const auto &f = station.frequencies(); const auto found = std::find(f.begin(), f.end(), change.second);
      const auto &mask = change.first == RealTzx || change.first == RealTzy ? t : change.first >= PTxx && change.first <= PTyy ? pt : z;
      if(found != f.end() && mask[type_to_column_table.at(change.first)][found-f.begin()]) result[entry.first].push_back(change);
    }
  }
  return result;
}
void apply(MTSurveyData &survey, const Edits &edits) {
  for(const auto &entry : edits) for(const auto &edit : entry.second)
    survey.get_station_data(entry.first).set_data_mask(edit.first, edit.second, false);
}
}
