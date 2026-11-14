#pragma once
#include "core/Datum.h"
#include <QWidget>
#include <QVector>
#include <array>

class MTSurveyData;
class QCustomPlot;
class QCPGraph;

// The same edit plan drives the table, plotted markers and accepted changes.
using GroupMaskEdits = std::map<std::string, std::vector<std::pair<RealDataType, double>>>;

class GroupMaskPlots : public QWidget {
public:
  GroupMaskPlots(QWidget *parent, const MTSurveyData &survey, const std::set<std::string> &stations);
  void setEdits(const GroupMaskEdits &edits, bool unmask);
private:
  struct Curve {
    std::string station;
    RealDataType type;
    unsigned pane, component;
    QCPGraph *observed, *masked, *pending;
    QVector<double> periods, values;
    std::vector<double> frequencies;
  };
  void fitAxes();
  std::array<QCustomPlot *, 4> plots;
  std::vector<Curve> curves;
};
