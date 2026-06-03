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

#ifndef TIPPER_PLOT_H
#define TIPPER_PLOT_H

#include "MTDataPlot.h"

class TipperPlot: public MTDataPlot
{
public:
  TipperPlot(QCustomPlot* plot);

  void set_observed_data(MTStationData &data, bool rescaleAxes = true);
  void set_predicted_data(const MTStationData &data, bool rescaleAxes = false);
  void set_arrow_mode(bool on);
  bool arrow_mode() const;

private:
  std::vector<RealDataType> get_graph_data_types(const QCPGraph *graph) const override;
  void set_layout();
  void clear_arrow_items();
  void clear_graph_data();
  void update_legend();
  void update_axis_style();
  void set_arrow_selection_data(const std::vector<std::vector<bool>> &mask,
                                const dvector &frequencies);
  void draw_arrow(double period, double xComponent, double yComponent, const QPen &pen);
  void draw_reference_arrow();
  void draw_tipper_arrows(const std::vector<dvector> &tipper,
                          const std::vector<std::vector<bool>> *mask,
                          const dvector &frequencies,
                          bool predicted);

  bool m_arrowMode;
  std::vector<QCPItemLine*> m_arrowItems;
  std::vector<QCPAbstractItem*> m_referenceItems;
};

#endif // PHASE_PLOT_H
