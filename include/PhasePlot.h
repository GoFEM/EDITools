#ifndef PHASE_PLOT_H
#define PHASE_PLOT_H

#include "MTDataPlot.h"

class PhasePlot: public MTDataPlot
{
public:
  PhasePlot(QCustomPlot* plot);

  void set_observed_data(MTStationData &data, bool rescaleAxes = true);
  void set_predicted_data(const MTStationData &data, bool rescaleAxes = false);
  void set_phase_wrap(bool on);
  bool phase_wrap() const;

private:
  void set_layout();
  void wrap_yx_yy_to_first_quadrant(std::vector<dvector> &phase) const;
  void set_phase_graph_data(const std::vector<std::vector<bool>> &mask,
                            const std::vector<dvector> &phase,
                            const std::vector<dvector> &phase_err,
                            const std::vector<dvector> &appRes,
                            const dvector &frequencies);
  void set_phase_graph_responses(const std::vector<dvector> &phase,
                                 const std::vector<dvector> &appRes,
                                 const dvector &frequencies);

  bool m_phaseWrap;
};

#endif // PHASE_PLOT_H
