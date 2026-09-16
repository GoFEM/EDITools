#include "include/FitStatistics.h"
#include <algorithm>

namespace FitStatistics {
std::string component_name(RealDataType type)
{
  const auto name = Datum::convert_type_to_string(type);
  if(name.compare(0, 4, "Real") == 0 || name.compare(0, 4, "Imag") == 0) return name.substr(4);
  if(name.compare(0, 3, "Rho") == 0) return "Rho " + name.substr(4);
  if(name.compare(0, 3, "Phs") == 0) return "Phase " + name.substr(4);
  return name;
}

Result compare(const MTSurveyData &observed, const MTSurveyData &response, ErrorSource errors)
{
  Result result;
  auto add = [&](const MTResponseData::Scalar &row) {
    ++result.available;
    if(!observed.is_station_present(row.receiver)) return;
    const auto &station = observed.get_station_data(row.receiver);
    const auto &frequencies = station.frequencies();
    auto match = std::find(frequencies.begin(), frequencies.end(), row.frequency);
    if(match == frequencies.end()) {
      double closest = 1e-3;
      for(auto it = frequencies.begin(); it != frequencies.end(); ++it) {
        const double distance = std::abs(*it - row.frequency) / row.frequency;
        if(distance < closest) { closest = distance; match = it; }
      }
    }
    if(match == frequencies.end()) return;
    double value, observed_error;
    if(!station.scalar_value(row.type, std::distance(frequencies.begin(), match), value, observed_error)) return;
    const double error = errors == ErrorSource::Response ? row.error : observed_error;
    if(!std::isfinite(error) || error <= 0 || !std::isfinite(row.value)) return;
    double difference = row.value - value;
    if(row.type == PhsZxx || row.type == PhsZxy || row.type == PhsZyx || row.type == PhsZyy)
      difference = std::remainder(difference, 360.);
    const double residual = difference / error;
    if(!std::isfinite(residual)) return;
    // Observed frequency keys align small differences in output precision.
    const double period = 1. / *match;
    result.total.add(residual);
    result.residuals.push_back(residual);
    result.periods[period].add(residual);
    result.components[component_name(row.type)].add(residual);
    result.stations[row.receiver].add(residual);
    result.station_periods[row.receiver][period].add(residual);
  };

  if(!response.response_observations().empty()) {
    for(const auto &row: response.response_observations()) add(row);
  } else {
    // Older projects store component arrays only. Positive component errors
    // provide the best available indication of explicitly supplied scalars.
    for(const auto &name: response.get_stations_names()) {
      const auto &station = response.get_station_data(name);
      for(unsigned f = 0; f < station.frequencies().size(); ++f)
        for(const auto &component: type_to_column_table) {
          double value, error;
          if(station.scalar_value(component.first, f, value, error) && std::isfinite(error) && error > 0)
            add({station.frequencies()[f], name, component.first, value, error});
        }
    }
  }
  return result;
}
}
