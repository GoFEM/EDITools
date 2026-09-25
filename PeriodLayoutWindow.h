#ifndef PERIOD_LAYOUT_WINDOW_H
#define PERIOD_LAYOUT_WINDOW_H

#include "include/PeriodResampling.h"
#include "qcustomplot.h"
#include <QDialog>
#include <functional>

class QComboBox;
class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;
class QCheckBox;
class QPushButton;
class QLabel;
class QTableWidget;

class PeriodLayoutWindow : public QDialog {
public:
  explicit PeriodLayoutWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> createSurvey);
  void setData(const std::shared_ptr<MTSurveyData> &survey);
private:
  std::vector<double> targetPeriods() const;
  void invalidate();
  void preview();
  void updatePlots();
  void drawLayout(QCustomPlot *plot, bool target);
  void exportReport(const PeriodResampling::Info &info);
  std::shared_ptr<MTSurveyData> survey, result;
  std::function<void(std::shared_ptr<MTSurveyData>)> createSurvey;
  QComboBox *gridMode, *reference, *group, *component;
  QLineEdit *minimum, *maximum, *custom;
  QSpinBox *density;
  QDoubleSpinBox *gap;
  QCheckBox *barriers;
  QPushButton *apply, *report, *openSource, *sourceReport;
  QLabel *summary, *sourceSummary, *history;
  QCustomPlot *sourcePlot, *targetPlot, *coverage;
  QTableWidget *counts;
};
#endif
