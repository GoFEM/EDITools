#ifndef PHASE_PLOT_H
#define PHASE_PLOT_H

#include "MTDataPlot.h"

class PhasePlot: public MTDataPlot
{
public:
  PhasePlot(QCustomPlot* plot);

  void set_observed_data(MTStationData &data, bool rescaleAxes = true);
  void set_predicted_data(const MTStationData &data, bool rescaleAxes = false);

private:
  void set_layout();
};

#endif // PHASE_PLOT_H
