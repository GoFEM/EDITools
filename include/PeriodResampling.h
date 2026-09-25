#ifndef PERIOD_RESAMPLING_H
#define PERIOD_RESAMPLING_H

#include "MTSurveyData.h"

namespace PeriodResampling {
struct Options {
  double maximumRatio = 2.;
  bool maskedBarriers = true;
};
struct Sample {
  double period;
  Status status;
};
const char *component_name(unsigned component);
const char *status_name(Status status);
bool available(Status status);
std::vector<double> validate_periods(std::vector<double> periods);
std::vector<double> logarithmic_grid(double minimum, double maximum, unsigned pointsPerDecade);
std::vector<Sample> layout(const MTStationData &station, unsigned component);
// Returns an independent survey; source data and an interpolation audit are embedded in it.
std::shared_ptr<MTSurveyData> resample(const MTSurveyData &source,
                                     const std::vector<double> &periods, const Options &options);
}
#endif
