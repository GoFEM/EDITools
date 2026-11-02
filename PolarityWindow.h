#pragma once
#include "include/PolarityAnalysis.h"
#include <QDialog>
#include <functional>
#include <memory>

class MTSurveyData;
class SurveyHistory;
class QTableWidget;
class QSpinBox;
class QDoubleSpinBox;
class QCheckBox;
class QPushButton;
class QLabel;
class QCustomPlot;

class PolarityWindow : public QDialog {
public:
  PolarityWindow(QWidget *parent, SurveyHistory &history, std::function<void()> changed);
  void setData(std::shared_ptr<MTSurveyData> data, const QString &selected = {});
private:
  void analyse();
  void select();
  void preview();
  void apply();
  std::array<bool, 5> selection() const;
  std::shared_ptr<MTSurveyData> survey;
  SurveyHistory &history;
  std::function<void()> changed;
  std::vector<PolarityAnalysis::Station> results;
  QTableWidget *stations;
  QSpinBox *neighbours;
  QDoubleSpinBox *radius, *tolerance;
  std::array<QCheckBox *, 5> channels;
  std::array<QCustomPlot *, 2> plots;
  QPushButton *applyButton, *suggestion;
  QLabel *status;
};
