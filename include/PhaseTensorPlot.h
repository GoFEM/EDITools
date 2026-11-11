/*
 * The EDI Tools application.
 *
 * Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef PHASE_TENSOR_PLOT_H
#define PHASE_TENSOR_PLOT_H

#include "MTDataPlot.h"

class PhaseTensorPlot: public MTDataPlot
{
public:
  PhaseTensorPlot(QCustomPlot* plot);

  void set_observed_data(MTStationData &data, bool rescaleAxes = true);
  void set_predicted_data(const MTStationData &data, bool rescaleAxes = false);
  void set_ellipse_mode(bool on);
  void set_view_style(const TFViewStyle &style) override;
  void clear_observed_data() override;
  void clear_predicted_data() override;
  void showPointToolTip(QMouseEvent *event) override;

private:
  void set_layout();
  void apply_component_visibility() override;
  void set_legend_component_visible(unsigned component, bool visible) override;
  std::vector<RealDataType> get_graph_data_types(const QCPGraph *graph) const override;
  void set_ellipses(const MTStationData &data, bool predicted);
  class Ellipses;
  Ellipses *m_ellipses = nullptr;
  QCPColorScale *m_colorScale;
  bool m_ellipseMode = false;
};

#endif // PHASE_TENSOR_PLOT_H
