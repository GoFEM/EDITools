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

#ifndef MAP_PLOT_H
#define MAP_PLOT_H

#include "qcustomplot.h"

#include "MTSurveyData.h"

class MapPlot
{
public:
  MapPlot(QCustomPlot* plot, QWidget *parent);

  void set_data(const std::vector<std::array<double, 3>> &locations);
  void set_selected_points(const std::vector<std::array<double, 3>> &locations);
  void get_point_value(const unsigned idx, double &key, double &value) const;

private:
  void set_layout();

private:
  QCustomPlot *m_plot;
  QWidget *m_parent;
};

#endif // MAP_PLOT_H
