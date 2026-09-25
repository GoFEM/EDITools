#ifndef EDI_PERIOD_MERGE_H
#define EDI_PERIOD_MERGE_H

#include "MTSurveyData.h"

namespace EDIPeriodMerge {
struct Change { std::string station; unsigned index; double before, after; };
struct Group {
  double minimum = 0., maximum = 0., target = 0.;
  unsigned periods = 0, stations = 0;
  bool collision = false, existingConflict = false;
};
struct Plan {
  unsigned before = 0, after = 0, samples = 0, invalid = 0, mergeable = 0, blocked = 0;
  double minimum = std::numeric_limits<double>::infinity(), maximum = 0.;
  std::vector<Group> groups;
  std::vector<Change> changes;
};
// Only incoming stations may change. Groups have a bounded total relative span,
// not a chain of pairwise matches. Existing survey periods are fixed anchors.
Plan analyze(const MTSurveyData &survey, const std::set<std::string> &incoming, double tolerance);
// Relabel frequencies without averaging/interpolation. Keep Z/T, errors and
// masks; update apparent resistivity and its error for the changed frequency.
void apply(MTSurveyData &survey, const Plan &plan);
}
#endif
