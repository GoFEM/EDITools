#ifndef FIT_STATISTICS_H
#define FIT_STATISTICS_H

#include "MTSurveyData.h"
#include <cmath>
#include <limits>

namespace FitStatistics {
enum class ErrorSource { Response, Observed };
struct Bin {
  std::size_t count = 0;
  double squares = 0.;
  void add(double residual) { ++count; squares += residual * residual; }
  double rms() const { return count ? std::sqrt(squares / count) : std::numeric_limits<double>::quiet_NaN(); }
};
struct Result {
  Bin total;
  std::vector<double> residuals;
  std::map<double, Bin> periods;
  std::map<std::string, Bin> components, stations;
  std::map<std::string, std::map<double, Bin>> station_periods;
  std::size_t available = 0;
};
std::string component_name(RealDataType type);
Result compare(const MTSurveyData &observed, const MTSurveyData &response, ErrorSource errors);
}
#endif
