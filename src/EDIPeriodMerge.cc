#include "include/EDIPeriodMerge.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

struct EDIPeriodMergeAccess {
  static void relabel(MTStationData &station, unsigned index, double frequency) {
    const double ratio = station.freqs.at(index) / frequency;
    station.freqs[index] = frequency;
    for(unsigned c = 0; c < 4; ++c) {
      station.Rho[c][index] *= ratio;
      station.Rho_err[c][index] *= ratio;
    }
  }
};

namespace EDIPeriodMerge {
Plan analyze(const MTSurveyData &survey, const std::set<std::string> &incoming, double tolerance)
{
  if(!std::isfinite(tolerance) || tolerance < 0. || tolerance > .05)
    throw std::invalid_argument("Period tolerance must be between 0 and 5%.");
  struct Point { double period, frequency; std::string station; unsigned index; bool incoming; };
  std::vector<Point> points;
  std::set<double> before, after;
  Plan plan;
  for(const auto &name: survey.get_stations_names()) {
    const auto &freqs = survey.get_station_data(name).frequencies();
    for(unsigned i = 0; i < freqs.size(); ++i) {
      const double period = 1. / freqs[i];
      if(!std::isfinite(period) || period <= 0.) { ++plan.invalid; continue; }
      points.push_back({period, freqs[i], name, i, incoming.count(name) != 0});
      before.insert(period); ++plan.samples;
      plan.minimum = std::min(plan.minimum, period); plan.maximum = std::max(plan.maximum, period);
    }
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) {
    if(a.period != b.period) return a.period < b.period;
    if(a.station != b.station) return a.station < b.station;
    return a.index < b.index;
  });
  for(unsigned begin = 0; begin < points.size();) {
    unsigned end = begin + 1;
    while(end < points.size() && (points[end].period == points[begin].period ||
          (points[end].period - points[begin].period) / points[begin].period < tolerance)) ++end;
    std::map<std::string, unsigned> stations;
    std::map<double, double> periods, anchors;
    bool hasIncoming = false;
    for(unsigned i = begin; i < end; ++i) {
      ++stations[points[i].station]; periods.emplace(points[i].period, points[i].frequency);
      if(!points[i].incoming) anchors.emplace(points[i].period, points[i].frequency);
      hasIncoming |= points[i].incoming;
    }
    Group group; group.minimum = points[begin].period; group.maximum = points[end - 1].period;
    group.periods = periods.size(); group.stations = stations.size();
    group.collision = std::any_of(stations.begin(), stations.end(), [](const auto &s) { return s.second > 1; });
    group.existingConflict = anchors.size() > 1;
    auto representative = periods.begin(); std::advance(representative, (periods.size() - 1) / 2);
    const double frequency = anchors.empty() ? representative->second : anchors.begin()->second;
    group.target = 1. / frequency;
    const bool merge = hasIncoming && periods.size() > 1 && !group.collision && !group.existingConflict;
    if(hasIncoming && (periods.size() > 1 || group.collision)) {
      plan.groups.push_back(group);
      if(merge) ++plan.mergeable; else ++plan.blocked;
    }
    for(unsigned i = begin; i < end; ++i) {
      const auto &point = points[i];
      if(merge && point.incoming && point.frequency != frequency) {
        plan.changes.push_back({point.station, point.index, point.frequency, frequency});
        after.insert(1. / frequency);
      } else after.insert(point.period);
    }
    begin = end;
  }
  plan.before = before.size(); plan.after = after.size();
  return plan;
}
void apply(MTSurveyData &survey, const Plan &plan)
{
  for(const auto &change: plan.changes)
    if(survey.get_station_data(change.station).frequencies().at(change.index) != change.before)
      throw std::invalid_argument("The survey changed after the merge preview.");
  for(const auto &change: plan.changes)
    EDIPeriodMergeAccess::relabel(survey.get_station_data(change.station), change.index, change.after);
}
}
