#ifndef EDITOOLS_SURVEY_REPORT_H
#define EDITOOLS_SURVEY_REPORT_H

#include "MTSurveyData.h"
#include <QString>
#include <functional>

class QWidget;
namespace SurveyReport {
enum class PeriodState { Full, Partial, Masked, Missing };
struct PeriodSummary {
  // Rows: impedance, tipper, phase tensor. Columns follow PeriodState.
  std::array<std::array<unsigned, 4>, 3> counts{};
  std::map<double, std::array<PeriodState, 3>> periods;
  double minimum = std::numeric_limits<double>::infinity(), maximum = 0.;
};
// Each positive, stored period counts once per data family. Full requires every
// scalar; Partial requires any usable scalar. Error magnitudes are not checked.
PeriodSummary summarize_periods(const MTStationData &station);
struct Options {
  struct Axis {
    bool autoscale = true;
    double lower = 0., upper = 1.;
  };
  QString title = "Survey data report", source, responseName;
  bool includeDisabled = true, errorBars = true, phaseWrap = true;
  std::array<Axis, 4> axes;
  std::array<std::array<bool, 4>, 4> components{{{{true, true, true, true}}, {{true, true, true, true}},
                                             {{true, true, true, true}}, {{true, true, true, true}}}};
};
// One overview page, then one page per included station. Returns false on cancellation.
// Publication is atomic: failure/cancellation preserves any existing destination.
bool write_pdf(const QString &path, const MTSurveyData &survey, const MTSurveyData *response,
               const Options &options, const std::function<bool(unsigned, unsigned)> &progress = {});
void show_dialog(QWidget *parent, const MTSurveyData &survey, const std::map<std::string, MTSurveyData> &responses,
                 const QString &selectedResponse, const Options &defaults, QString &lastDirectory);
}
#endif
