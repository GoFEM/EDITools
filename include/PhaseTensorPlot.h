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

private:
  void set_layout();
};

#endif // PHASE_TENSOR_PLOT_H
