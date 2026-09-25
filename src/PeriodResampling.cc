#include "include/PeriodResampling.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace {
using namespace PeriodResampling;
const double missing = std::numeric_limits<double>::quiet_NaN();
bool same_period(double a, double b) { return std::abs(a - b) <= 32. * std::numeric_limits<double>::epsilon() * std::max(a, b); }
bool finite(std::complex<double> z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
bool valid_error(double e) { return std::isfinite(e) && e >= 0.; }
struct Point {
  double period, error;
  std::complex<double> value;
  bool enabled;
  unsigned index;
  bool hasValue;
};
struct Value {
  std::complex<double> value{missing, missing};
  double error = missing;
  Record record;
};
Value interpolate(const std::vector<Point> &points, double target, const Options &options)
{
  Value result;
  result.record.period = target;
  for(const auto &p: points) if(same_period(p.period, target) && (!p.enabled || p.hasValue)) {
    result.value = p.value; result.error = p.error;
    result.record.lower = result.record.upper = p.period;
    result.record.status = !p.enabled ? Status::Masked : (finite(p.value) ? Status::Existing : Status::Missing);
    return result;
  }
  const Point *left = nullptr, *right = nullptr;
  bool hasData = false;
  for(const auto &p: points) if(p.enabled && finite(p.value)) {
    hasData = true;
    if(p.period < target) left = &p;
    else { right = &p; break; }
  }
  if(!left || !right) {
    result.record.status = hasData ? Status::Outside : Status::Missing;
    return result;
  }
  result.record.lower = left->period; result.record.upper = right->period;
  if(options.maskedBarriers)
    for(const auto &p: points) if(!p.enabled && p.period > left->period && p.period < right->period) {
      result.record.status = Status::Masked; return result;
    }
  const double span = std::log(right->period) - std::log(left->period);
  if(span > std::log(options.maximumRatio) + 32. * std::numeric_limits<double>::epsilon()) {
    result.record.status = Status::Gap; return result;
  }
  if(!valid_error(left->error) || !valid_error(right->error)) {
    result.record.status = Status::InvalidError; return result;
  }
  const double weight = (std::log(target) - std::log(left->period)) / span;
  result.value = (1. - weight) * left->value + weight * right->value;
  result.error = (1. - weight) * left->error + weight * right->error;
  result.record.status = finite(result.value) && valid_error(result.error) ? Status::Interpolated : Status::Missing;
  return result;
}
}

// Keep access to the station's raw errors and masks here so exact samples can be copied losslessly.
struct PeriodResamplingAccess {
  static std::vector<Point> points(const MTStationData &s, unsigned c) {
    std::vector<Point> result;
    for(unsigned i = 0; i < s.freqs.size(); ++i) {
      if(!std::isfinite(s.freqs[i]) || s.freqs[i] <= 0.) continue;
      const double period = 1. / s.freqs[i];
      if(!std::isfinite(period) || period <= 0.) continue;
      if(c < 4) result.push_back({period, s.Z_err_floor[c][i], s.Z[c][i], s.Z_mask[c][i], i,
                                 std::isfinite(s.Z[c][i].real()) || std::isfinite(s.Z[c][i].imag())});
      else if(c < 6) result.push_back({period, s.T_err_floor[c-4][i], s.T[c-4][i], s.T_mask[c-4][i], i,
                                      std::isfinite(s.T[c-4][i].real()) || std::isfinite(s.T[c-4][i].imag())});
      else result.push_back({period, s.PT_err[c-6][i], {s.PT[c-6][i], 0.}, s.PT_mask[c-6][i], i, std::isfinite(s.PT[c-6][i])});
    }
    std::sort(result.begin(), result.end(), [](const Point &a, const Point &b) { return a.period < b.period; });
    return result;
  }
  static void copy_sample(MTStationData &out, unsigned to, const MTStationData &in, unsigned from,
                          std::array<Value, 10> &values) {
    for(unsigned c = 0; c < 10; ++c) {
      bool preserve = false;
      if(c < 4) {
        preserve = !in.Z_mask[c][from] || std::isfinite(in.Z[c][from].real()) || std::isfinite(in.Z[c][from].imag());
        if(preserve) {
          out.Z[c][to] = in.Z[c][from]; out.Z_err[c][to] = in.Z_err[c][from];
          out.Z_err_floor[c][to] = in.Z_err_floor[c][from]; out.Z_mask[c][to] = in.Z_mask[c][from];
        }
        if(!in.Z_mask[c][from] || std::isfinite(in.Rho[c][from])) {
          out.Rho[c][to] = in.Rho[c][from]; out.Rho_err[c][to] = in.Rho_err[c][from]; out.Z_mask[c][to] = in.Z_mask[c][from];
        }
        if(!in.Z_mask[c][from] || std::isfinite(in.Phs[c][from])) {
          out.Phs[c][to] = in.Phs[c][from]; out.Phs_err[c][to] = in.Phs_err[c][from]; out.Z_mask[c][to] = in.Z_mask[c][from];
        }
      } else if(c < 6) {
        const auto k = c - 4;
        preserve = !in.T_mask[k][from] || std::isfinite(in.T[k][from].real()) || std::isfinite(in.T[k][from].imag());
        if(preserve) {
          out.T[k][to] = in.T[k][from]; out.T_err[k][to] = in.T_err[k][from];
          out.T_err_floor[k][to] = in.T_err_floor[k][from]; out.T_mask[k][to] = in.T_mask[k][from];
        }
      } else {
        const auto k = c - 6;
        preserve = !in.PT_mask[k][from] || std::isfinite(in.PT[k][from]);
        if(preserve) { out.PT[k][to] = in.PT[k][from]; out.PT_err[k][to] = in.PT_err[k][from]; out.PT_mask[k][to] = in.PT_mask[k][from]; }
      }
      if(preserve) {
        auto &record = values[c].record;
        record.lower = record.upper = 1. / in.freqs[from];
        const bool enabled = c < 4 ? in.Z_mask[c][from] : (c < 6 ? in.T_mask[c-4][from] : in.PT_mask[c-6][from]);
        const bool complete = c < 4 ? finite(in.Z[c][from]) : (c < 6 ? finite(in.T[c-4][from]) : std::isfinite(in.PT[c-6][from]));
        record.status = !enabled ? Status::Masked : (complete ? Status::Existing : Status::Missing);
      }
    }
  }
  static std::shared_ptr<MTSurveyData> run(const MTSurveyData &source, const std::vector<double> &targets, const Options &options) {
    auto result = std::make_shared<MTSurveyData>(source);
    result->m_survey_name = source.get_survey_name() + " (resampled)";
    result->m_response_observations.clear();
    result->m_resampling_source = std::make_shared<MTSurveyData>(source);
    result->m_resampling_info = {};
    auto &info = result->m_resampling_info;
    info.method = "Linear real/imaginary transfer functions versus log(period); weighted endpoint standard errors";
    info.maximumRatio = options.maximumRatio; info.maskedBarriers = options.maskedBarriers; info.periods = targets;
    std::set<double> frequencies;
    for(double p: targets) frequencies.insert(1. / p);
    for(auto &entry: result->m_stations_data) {
      const auto &input = source.get_station_data(entry.first);
      std::array<std::vector<Point>, 10> series;
      bool hasImpedance = false;
      for(unsigned c = 0; c < 10; ++c) {
        series[c] = points(input, c);
        if(c < 4) for(const auto &p: series[c])
          hasImpedance |= std::isfinite(p.value.real()) || std::isfinite(p.value.imag());
      }
      MTStationData output;
      output.set_size(frequencies.size(), true); output.set_frequencies(frequencies);
      output.station_name = input.station_name; output.file_name = input.file_name;
      output.location = input.location; output.is_active = input.is_active; output.error_floor = input.error_floor;
      std::vector<std::array<Value, 10>> values(output.freqs.size());
      for(unsigned f = 0; f < output.freqs.size(); ++f) {
        const double target = 1. / output.freqs[f];
        for(unsigned c = 0; c < 10; ++c) {
          auto &v = values[f][c]; v = interpolate(series[c], target, options);
          v.record.station = entry.first; v.record.component = c;
          const bool enabled = available(v.record.status);
          if(c < 4) {
            output.Z[c][f] = enabled ? v.value : std::complex<double>(missing, missing);
            output.Z_err[c][f] = output.Z_err_floor[c][f] = v.error; output.Z_mask[c][f] = enabled;
          } else if(c < 6) {
            output.T[c-4][f] = enabled ? v.value : std::complex<double>(missing, missing);
            output.T_err[c-4][f] = output.T_err_floor[c-4][f] = v.error; output.T_mask[c-4][f] = enabled;
          }
        }
      }
      output.calculate_apparent_resistivity(); output.calculate_phase(); output.calculate_phase_tensor();
      output.propagate_rho_phase_error(); output.propagate_phase_tensor_error();
      for(unsigned f = 0; f < output.freqs.size(); ++f) {
        const double target = 1. / output.freqs[f];
        auto &v = values[f];
        for(unsigned c = 6; c < 10; ++c) {
          if(hasImpedance) {
            Record record = v[c].record;
            record.lower = target; record.upper = target; record.status = Status::Derived;
            for(unsigned z = 0; z < 4; ++z) {
              if(!available(v[z].record.status)) record.status = v[z].record.status;
              if(v[z].record.lower > 0.) record.lower = std::min(record.lower, v[z].record.lower);
              record.upper = std::max(record.upper, v[z].record.upper);
            }
            // Tensor masks are independent of impedance masks and must survive derivation.
            for(const auto &p: series[c]) if(!p.enabled &&
                (same_period(p.period, target) || (options.maskedBarriers && p.period >= record.lower && p.period <= record.upper)))
              record.status = Status::Masked;
            const double a = output.Z[0][f].real(), b = output.Z[1][f].real();
            const double cc = output.Z[2][f].real(), d = output.Z[3][f].real();
            const double scale = std::max({std::abs(a), std::abs(b), std::abs(cc), std::abs(d)});
            if(available(record.status) && (!(scale > 0.) || !std::isfinite(output.PT[c-6][f]) ||
               std::abs((a / scale) * (d / scale) - (b / scale) * (cc / scale)) < 1e-12)) record.status = Status::Missing;
            output.PT_mask[c-6][f] = available(record.status);
            if(!available(record.status)) output.PT[c-6][f] = missing;
            v[c].record = record;
          } else {
            output.PT[c-6][f] = available(v[c].record.status) ? v[c].value.real() : missing;
            output.PT_err[c-6][f] = v[c].error;
            output.PT_mask[c-6][f] = available(v[c].record.status);
            // A directly supplied tensor needs no impedance, whose masks must remain enabled for maps.
            output.Z_mask[c-6][f] = true;
          }
        }
        for(const auto &p: series[0]) if(same_period(p.period, target)) {
          copy_sample(output, f, input, p.index, v);
          break;
        }
        for(const auto &item: v) info.records.push_back(item.record);
      }
      entry.second = std::move(output);
    }
    return result;
  }
};

namespace PeriodResampling {
const char *component_name(unsigned c) {
  static const char *names[] = {"Zxx", "Zxy", "Zyx", "Zyy", "Tzx", "Tzy", "PTxx", "PTxy", "PTyx", "PTyy"};
  if(c >= 10) throw std::invalid_argument("Invalid component.");
  return names[c];
}
const char *status_name(Status s) {
  switch(s) {
    case Status::Existing: return "Existing";
    case Status::Interpolated: return "Interpolated";
    case Status::Derived: return "Derived from impedance";
    case Status::Masked: return "Masked / mask barrier";
    case Status::Outside: return "Outside component range";
    case Status::Gap: return "Gap too large";
    case Status::InvalidError: return "Invalid endpoint error";
    default: return "Missing / incomplete";
  }
}
bool available(Status s) { return s == Status::Existing || s == Status::Interpolated || s == Status::Derived; }
std::vector<double> validate_periods(std::vector<double> periods) {
  if(periods.empty() || periods.size() > 2000) throw std::invalid_argument("Specify between 1 and 2000 target periods.");
  for(double p: periods) if(!std::isfinite(p) || p <= 0. || !std::isfinite(1. / p) || 1. / p <= 0.)
    throw std::invalid_argument("Periods must be finite, positive seconds with representable frequencies.");
  std::sort(periods.begin(), periods.end());
  periods.erase(std::unique(periods.begin(), periods.end(), same_period), periods.end());
  return periods;
}
std::vector<double> logarithmic_grid(double minimum, double maximum, unsigned pointsPerDecade) {
  validate_periods({minimum, maximum});
  if(maximum < minimum || pointsPerDecade < 1) throw std::invalid_argument("Set ordered period limits and at least one point per decade.");
  const double steps = (std::log10(maximum) - std::log10(minimum)) * pointsPerDecade;
  if(steps > 1999.) throw std::invalid_argument("Grid exceeds 2000 periods; reduce density or range.");
  std::vector<double> result;
  for(unsigned i = 0; i <= static_cast<unsigned>(steps); ++i)
    result.push_back(std::exp(std::log(minimum) + std::log(10.) * double(i) / pointsPerDecade));
  if(!same_period(result.back(), maximum)) result.push_back(maximum);
  else result.back() = maximum;
  return validate_periods(result);
}
std::vector<Sample> layout(const MTStationData &station, unsigned c) {
  component_name(c);
  std::vector<Sample> result;
  for(const auto &p: PeriodResamplingAccess::points(station, c))
    result.push_back({p.period, !station.active() || !p.enabled ? Status::Masked : (finite(p.value) ? Status::Existing : Status::Missing)});
  return result;
}
std::shared_ptr<MTSurveyData> resample(const MTSurveyData &source, const std::vector<double> &periods, const Options &options) {
  const auto targets = validate_periods(periods);
  if(!std::isfinite(options.maximumRatio) || options.maximumRatio <= 1.) throw std::invalid_argument("Maximum bounding-period ratio must exceed 1.");
  if(source.n_stations() == 0) throw std::invalid_argument("Load a survey before resampling.");
  if(static_cast<double>(source.n_stations()) * targets.size() > 200000.) throw std::invalid_argument("Preview exceeds 200000 station-period cells; reduce the target grid.");
  return PeriodResamplingAccess::run(source, targets, options);
}
}
