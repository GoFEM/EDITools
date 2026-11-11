#include "include/ImpedanceAnalysis.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace ImpedanceAnalysis {
namespace {
const double pi = std::acos(-1.), mu = 4e-7 * pi;
using Complex = std::complex<double>;
bool finite(Complex z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
Complex root(Complex z) {
  auto value = std::sqrt(z);
  return value.imag() < 0. ? -value : value; // Explicit branch: phase in [0, 180] degrees.
}
}
const char *name(Quantity q) {
  switch(q) {
  case Quantity::XX: return "Zxx";
  case Quantity::XY: return "Zxy";
  case Quantity::YX: return "−Zyx";
  case Quantity::YY: return "Zyy";
  case Quantity::SSQ: return "SSQ";
  case Quantity::Berdichevsky: return "Berdichevsky (Zxy − Zyx)/2";
  case Quantity::Determinant: return "Determinant √det(Z)";
  }
  throw std::invalid_argument("Unknown impedance quantity.");
}
bool Sample::usable() const { return std::isfinite(period) && period > 0. && finite(impedance) && std::abs(impedance) > 0.; }
double Sample::resistivity() const { return usable() ? std::norm(impedance) * period / (2. * pi * mu) : NAN; }
double Sample::phase() const { return usable() ? std::arg(impedance) * 180. / pi : NAN; }
Sample sample(const MTStationData &station, Quantity quantity, unsigned f, double clockwiseDegrees) {
  Sample sample;
  if(f >= station.frequencies().size() || !std::isfinite(clockwiseDegrees)) return sample;
  const double frequency = station.frequencies()[f];
  if(!std::isfinite(frequency) || frequency <= 0.) return sample;
  sample.period = 1. / frequency;
  const std::array<RealDataType, 4> re{{RealZxx, RealZxy, RealZyx, RealZyy}}, im{{ImagZxx, ImagZxy, ImagZyx, ImagZyy}};
  const std::array<RealDataType, 4> rho{{RhoZxx, RhoZxy, RhoZyx, RhoZyy}}, phase{{PhsZxx, PhsZxy, PhsZyx, PhsZyy}};
  std::array<Complex, 4> z; std::array<double, 4> errors; std::array<bool, 4> valid;
  for(unsigned c = 0; c < 4; ++c) {
    double r = NAN, i = NAN, e = NAN, ei = NAN;
    valid[c] = station.scalar_value(re[c], f, r, e) && station.scalar_value(im[c], f, i, ei);
    z[c] = {r, i}; errors[c] = e;
  }
  const unsigned component = unsigned(quantity);
  const double angle = std::remainder(clockwiseDegrees, 180.) * pi / 180.;
  if(component < 4 && angle != 0.) {
    double c = std::cos(angle), s = std::sin(angle);
    if(std::abs(c) < 1e-14) c = 0.;
    if(std::abs(s) < 1e-14) s = 0.;
    const double rotation[2][2]{{c, s}, {-s, c}};
    Complex value{}; double variance = 0.;
    for(unsigned k = 0; k < 4; ++k) {
      const double weight = rotation[component / 2][k / 2] * rotation[component % 2][k % 2];
      if(weight == 0.) continue;
      if(!valid[k]) return sample;
      value += weight * z[k];
      variance += std::isfinite(errors[k]) && errors[k] >= 0. ? weight * weight * errors[k] * errors[k] : NAN;
    }
    sample.impedance = (quantity == Quantity::YX ? -1. : 1.) * value;
    sample.error = std::sqrt(variance);
    if(sample.usable()) {
      const double relative = sample.error / std::abs(sample.impedance);
      sample.rhoError = 2. * sample.resistivity() * relative;
      sample.phaseError = relative * 180. / pi;
    }
  } else if(component < 4) {
    if(valid[component]) { sample.impedance = (quantity == Quantity::YX ? -1. : 1.) * z[component]; sample.error = errors[component];
      double value;
      station.scalar_value(rho[component], f, value, sample.rhoError);
      station.scalar_value(phase[component], f, value, sample.phaseError); }
  } else {
    const bool all = std::all_of(valid.begin(), valid.end(), [](bool v) { return v; });
    std::array<Complex, 4> derivative{};
    if(quantity == Quantity::Berdichevsky && valid[1] && valid[2]) {
      sample.impedance = (z[1] - z[2]) * .5; derivative[1] = .5; derivative[2] = -.5;
    } else if(all && quantity == Quantity::SSQ) {
      sample.impedance = root((z[0]*z[0] + z[1]*z[1] + z[2]*z[2] + z[3]*z[3]) * .5);
      if(sample.usable()) for(unsigned c = 0; c < 4; ++c) derivative[c] = z[c] / (2. * sample.impedance);
    } else if(all && quantity == Quantity::Determinant) {
      sample.impedance = root(z[0]*z[3] - z[1]*z[2]);
      if(sample.usable()) derivative = {{z[3], -z[2], -z[1], z[0]}};
      if(sample.usable()) for(auto &d: derivative) d /= 2. * sample.impedance;
    }
    if(sample.usable()) {
      double variance = 0.;
      for(unsigned c = 0; c < 4; ++c) if(std::abs(derivative[c]) > 0.) {
        if(!std::isfinite(errors[c]) || errors[c] < 0.) { variance = NAN; break; }
        variance += std::norm(derivative[c]) * errors[c] * errors[c];
      }
      sample.error = std::sqrt(variance);
      const double relative=sample.error/std::abs(sample.impedance);
      sample.rhoError=2.*sample.resistivity()*relative; sample.phaseError=relative*180./pi;
    }
  }
  return sample;
}
std::vector<Sample> curve(const MTStationData &station, Quantity quantity) {
  std::vector<Sample> result;
  for(unsigned f = 0; f < station.frequencies().size(); ++f) {
    auto value = sample(station, quantity, f);
    if(value.period > 0.) result.push_back(value);
  }
  std::sort(result.begin(), result.end(), [](const Sample &a, const Sample &b) { return a.period < b.period; });
  return result;
}
std::vector<RegionalPoint> regionalMean(const std::vector<std::vector<Sample>> &curves, double tolerance, unsigned minimum) {
  if(!std::isfinite(tolerance) || tolerance < 0. || tolerance > .25 || minimum == 0)
    throw std::invalid_argument("Invalid regional averaging settings.");
  struct Entry { Sample value; unsigned station; };
  std::vector<Entry> entries;
  for(unsigned s = 0; s < curves.size(); ++s) for(const auto &v: curves[s])
    if(std::isfinite(v.period) && v.period > 0.) entries.push_back({v, s});
  std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) { return a.value.period < b.value.period; });
  std::vector<RegionalPoint> result;
  for(size_t begin = 0; begin < entries.size();) {
    size_t end = begin + 1;
    while(end < entries.size() && (entries[end].value.period == entries[begin].value.period ||
          (entries[end].value.period / entries[begin].value.period - 1.) <= tolerance)) ++end;
    std::map<unsigned, unsigned> counts;
    for(size_t j = begin; j < end; ++j) ++counts[entries[j].station];
    std::vector<double> logRho, phase, periods;
    Complex direction{};
    for(size_t j = begin; j < end; ++j) if(counts[entries[j].station] == 1) {
      const auto &v = entries[j].value;
      if(!v.usable() || !std::isfinite(v.resistivity()) || v.resistivity() <= 0.) continue;
      logRho.push_back(std::log(v.resistivity())); phase.push_back(std::arg(v.impedance));
      periods.push_back(v.period); direction += std::polar(1., phase.back());
    }
    // Keep a row even below the count threshold, making holes explicit in plots.
    RegionalPoint point; point.period = entries[begin + (end - begin - 1) / 2].value.period;
    point.stations = logRho.size();
    if(!periods.empty()) point.period = periods[(periods.size()-1)/2];
    if(point.stations >= minimum) {
      double mean = 0.; for(double v: logRho) mean += v / point.stations;
      point.resistivity = std::exp(mean);
      if(point.stations > 1) {
        double variance = 0.; for(double v: logRho) variance += (v - mean) * (v - mean);
        point.logRhoSpread = std::sqrt(variance / (point.stations - 1));
      }
      // Choose a common phase branch before taking the complex geometric mean.
      if(std::abs(direction) > 1e-8 * point.stations) {
        const double center = std::arg(direction);
        double average = 0.; for(double v: phase) average += (center + std::remainder(v - center, 2*pi)) / point.stations;
        point.phase = std::remainder(average, 2*pi) * 180. / pi;
      }
    }
    result.push_back(point); begin = end;
  }
  return result;
}
namespace {
// Integral of 1/(lambda + i*omega) over a cell, with its real part referenced
// to omega=1. This removes the divergent constant in the infinite upper tail.
Complex kernel(double a, double b, double omega) {
  auto real = [omega](double edge) { return std::log(std::hypot(edge, omega)) - std::log(std::hypot(edge, 1.)); };
  const double hiReal = std::isinf(b) ? 0. : real(b);
  const double hiImag = std::isinf(b) ? 0. : std::atan2(omega, b);
  return {hiReal - real(a), hiImag - std::atan2(omega, a)};
}
struct Solution { Eigen::VectorXd x; bool converged = false; };
// Active-set bounded least squares, implemented here; no external solver code.
Solution boundedLeastSquares(const Eigen::MatrixXd &a, const Eigen::VectorXd &b) {
  const int n = a.cols();
  Solution result; result.x = Eigen::VectorXd::Zero(n);
  std::vector<bool> free(n, false);
  const Eigen::VectorXd norms = a.colwise().norm();
  const double tolerance = 1e-10 * std::max(1., b.norm());
  for(unsigned iteration = 0; iteration < 5000; ++iteration) {
    const Eigen::VectorXd gradient = a.transpose() * (b - a * result.x);
    int enter = -1; double largest = tolerance;
    for(int j = 0; j < n; ++j) if(!free[j] && norms[j] > 1e-14) {
      const double violation = (result.x[j] == 0. ? gradient[j] : -gradient[j]) / norms[j];
      if(violation > largest) { largest = violation; enter = j; }
    }
    if(enter < 0) {
      result.converged = true;
      for(int j = 0; j < n; ++j) if(free[j] && norms[j] > 1e-14 && std::abs(gradient[j])/norms[j] > 10*tolerance)
        result.converged = false;
      break;
    }
    free[enter] = true;
    bool solved = false;
    for(int inner = 0; inner < 2*n + 1; ++inner) {
      std::vector<int> indices; Eigen::VectorXd rhs = b;
      for(int j = 0; j < n; ++j) { if(free[j]) indices.push_back(j); else if(result.x[j] != 0.) rhs -= a.col(j); }
      if(indices.empty()) { solved = true; break; }
      Eigen::MatrixXd active(a.rows(), indices.size());
      for(unsigned j = 0; j < indices.size(); ++j) active.col(j) = a.col(indices[j]);
      auto qr = active.colPivHouseholderQr(); qr.setThreshold(1e-12);
      const Eigen::VectorXd candidate = qr.solve(rhs);
      if(!candidate.allFinite()) return result;
      double step = 1.; int hit = -1; double bound = 0.;
      for(unsigned j = 0; j < indices.size(); ++j) {
        const int k = indices[j]; const double delta = candidate[j] - result.x[k];
        if(candidate[j] < 0. || candidate[j] > 1.) {
          const double limit = candidate[j] < 0. ? 0. : 1.;
          const double fraction = (limit - result.x[k]) / delta;
          if(fraction < step) { step = std::max(0., fraction); hit = k; bound = limit; }
        }
      }
      for(unsigned j = 0; j < indices.size(); ++j) {
        const int k = indices[j]; result.x[k] += step * (candidate[j] - result.x[k]);
        result.x[k] = std::max(0., std::min(1., result.x[k]));
      }
      if(hit < 0) { solved = true; break; }
      result.x[hit] = bound; free[hit] = false;
    }
    if(!solved) return result;
  }
  return result;
}
struct Fit {
  std::vector<double> edges;
  Eigen::VectorXd density;
  double offset = 0., frequencyScale = 1., chi = 0.; bool converged = false;
  Complex logC(double period) const {
    const double w = 2*pi / period / frequencyScale;
    Complex value(offset, 0.);
    for(int j = 0; j < density.size(); ++j) value += density[j] * kernel(edges[j], edges[j+1], w);
    return value;
  }
  Sample evaluate(double period) const {
    Sample value; value.period = period;
    value.impedance = Complex(0., mu * 2*pi / period) * std::exp(logC(period));
    return value;
  }
};
Fit fit(const std::vector<Sample> &samples, const std::vector<std::array<double,2>> &errors, unsigned density) {
  Fit result;
  const double low = 2*pi / samples.back().period, high = 2*pi / samples.front().period;
  result.frequencyScale = std::sqrt(low) * std::sqrt(high);
  const double a = std::log10(low / result.frequencyScale) - 4., b = std::log10(high / result.frequencyScale) + 4.;
  const unsigned cells = std::max(1u, unsigned(std::ceil((b-a) * density)));
  result.edges.push_back(0.);
  for(unsigned j = 0; j <= cells; ++j) result.edges.push_back(std::pow(10., a + (b-a)*j/cells));
  result.edges.push_back(std::numeric_limits<double>::infinity());
  const unsigned n = samples.size(), m = result.edges.size() - 1;
  Eigen::MatrixXd aMatrix(2*n, m); Eigen::VectorXd data(2*n), constant = Eigen::VectorXd::Zero(2*n);
  for(unsigned i = 0; i < n; ++i) {
    const double omega = 2*pi / samples[i].period;
    data[2*i] = std::log(std::abs(samples[i].impedance) / (mu * omega)) / errors[i][0];
    data[2*i+1] = (std::arg(samples[i].impedance) - pi/2.) / errors[i][1];
    constant[2*i] = 1. / errors[i][0];
    for(unsigned j = 0; j < m; ++j) {
      const auto k = kernel(result.edges[j], result.edges[j+1], omega / result.frequencyScale);
      aMatrix(2*i,j) = k.real()/errors[i][0]; aMatrix(2*i+1,j) = k.imag()/errors[i][1];
    }
  }
  const double norm = constant.squaredNorm();
  const Eigen::RowVectorXd means = constant.transpose() * aMatrix / norm;
  const double mean = constant.dot(data) / norm;
  Eigen::MatrixXd centered = aMatrix - constant * means;
  const auto solution = boundedLeastSquares(centered, data - constant * mean);
  result.density = solution.x; result.converged = solution.converged;
  result.offset = mean - (means * result.density)[0];
  result.chi = (aMatrix * result.density + constant * result.offset - data).squaredNorm();
  return result;
}
}
std::vector<unsigned> rhoPlusOutliers(const RhoPlusResult &result, double threshold) {
  if(!std::isfinite(threshold) || threshold <= 0.)
    throw std::invalid_argument("Residual threshold must be positive and finite.");
  std::vector<unsigned> indices;
  if(!result.converged || !result.gridStable) return indices;
  if(result.rhoResidual.size()!=result.observed.size() || result.phaseResidual.size()!=result.observed.size())
    throw std::invalid_argument("Incomplete rho+ residuals.");
  for(unsigned i=0;i<result.observed.size();++i)
    if(std::isfinite(result.rhoResidual[i]) && std::isfinite(result.phaseResidual[i]) &&
       (std::abs(result.rhoResidual[i])>threshold || std::abs(result.phaseResidual[i])>threshold)) indices.push_back(i);
  return indices;
}
RhoPlusResult rhoPlus(const std::vector<Sample> &samples, const RhoPlusOptions &options) {
  if(options.intervalsPerDecade < 4 || options.intervalsPerDecade > 48)
    throw std::invalid_argument("Invalid rho+ fit settings.");
  RhoPlusResult result;
  for(const auto &s: samples) {
    if(!s.usable() || !std::isfinite(s.resistivity()) || s.resistivity()<=0.) { ++result.excluded; continue; }
    const double amplitudeError=s.rhoError/(2.*s.resistivity()), phaseError=s.phaseError*pi/180.;
    if(!std::isfinite(amplitudeError) || amplitudeError<=0. || !std::isfinite(phaseError) || phaseError<=0.) {
      ++result.invalidErrors; continue;
    }
    result.observed.push_back(s);
  }
  std::sort(result.observed.begin(), result.observed.end(), [](const Sample &a, const Sample &b) { return a.period < b.period; });
  if(result.observed.size() < 3 || result.observed.front().period == result.observed.back().period)
    throw std::invalid_argument("Rho+ needs at least three usable complex-impedance samples with positive propagated resistivity/phase errors spanning distinct periods.");
  std::vector<std::array<double,2>> errors;
  for(const auto &s:result.observed) errors.push_back({{s.rhoError/(2.*s.resistivity()),s.phaseError*pi/180.}});
  const auto coarse = fit(result.observed, errors, options.intervalsPerDecade);
  const auto fine = fit(result.observed, errors, 2 * options.intervalsPerDecade);
  result.converged = coarse.converged && fine.converged;
  result.refinementChange = std::abs(coarse.chi - fine.chi);
  result.gridStable = result.refinementChange <= .01 * std::max(1., fine.chi);
  result.chiSquare = fine.chi; result.normalizedRms = std::sqrt(fine.chi / (2 * result.observed.size()));
  for(unsigned i = 0; i < result.observed.size(); ++i) {
    const auto &s = result.observed[i]; const auto prediction = fine.evaluate(s.period);
    result.fitted.push_back(prediction);
    result.rhoResidual.push_back(std::log(prediction.resistivity() / s.resistivity()) / (2. * errors[i][0]));
    result.phaseResidual.push_back((prediction.phase() - s.phase()) * pi/180. / errors[i][1]);
  }
  const double lo = std::log(result.observed.front().period), hi = std::log(result.observed.back().period);
  for(unsigned i = 0; i <= 200; ++i) result.response.push_back(fine.evaluate(std::exp(lo + (hi-lo)*i/200.)));
  return result;
}
}
