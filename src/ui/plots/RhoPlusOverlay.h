#ifndef EDITOOLS_RHO_PLUS_OVERLAY_H
#define EDITOOLS_RHO_PLUS_OVERLAY_H

#include "core/ImpedanceAnalysis.h"
#include "qcustomplot.h"
#include <QTimer>
#include <functional>

// Owns only the fit cache and refresh timer; plots own the non-selectable curves.
class RhoPlusOverlay {
public:
  RhoPlusOverlay(QCustomPlot *rho, QCustomPlot *phase, std::function<void(const QString &)> status);
  void hide();
  void set_colors(const std::array<std::string, 8> &colors);
  void update(const MTStationData *station, bool enabled, bool phaseWrap,
              const std::array<bool,4> &rhoVisible, const std::array<bool,4> &phaseVisible,
              const ImpedanceAnalysis::RhoPlusOptions &options);
private:
  void calculate();
  void draw();
  struct Component {
    std::vector<ImpedanceAnalysis::Sample> input, response;
    bool pending = false;
    bool provisional = false;
    QString status;
  };
  QCustomPlot *rho, *phase;
  std::array<QCPCurve *,4> rhoLines, phaseLines;
  std::array<Component,4> components;
  std::array<bool,4> rhoVisible{}, phaseVisible{};
  ImpedanceAnalysis::RhoPlusOptions options;
  bool enabled = false, phaseWrap = true;
  QTimer timer;
  std::function<void(const QString &)> report;
};
#endif
