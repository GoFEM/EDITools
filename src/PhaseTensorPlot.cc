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

#include "include/PhaseTensorPlot.h"

PhaseTensorPlot::PhaseTensorPlot(QCustomPlot *plot):
  MTDataPlot(plot)
{
  set_layout();
}

void PhaseTensorPlot::set_observed_data(MTStationData &data, bool rescaleAxes)
{
  m_data = &data;

  std::vector<dvector> phase, phase_err;
  m_data->get_phase_tensor(phase, phase_err);
  const auto mask = m_data->phase_tensor_mask();

  set_graph_data(mask, phase, phase_err, data.frequencies());

  if(rescaleAxes)
  {
    m_plot->rescaleAxes();
    QCPRange xrange = m_plot->xAxis->range();
    m_plot->xAxis->setRange(xrange.lower / 2., xrange.upper * 2.);
  }

  m_plot->replot();
}

void PhaseTensorPlot::set_predicted_data(const MTStationData &data, bool rescaleAxes)
{
  std::vector<dvector> phase, phase_err;
  data.get_phase_tensor(phase, phase_err);

  set_graph_responses(phase, data.frequencies());

  if(rescaleAxes)
  {
    m_plot->rescaleAxes();
    QCPRange xrange = m_plot->xAxis->range();
    m_plot->xAxis->setRange(xrange.lower / 2., xrange.upper * 2.);
  }

  m_plot->replot();
}

void PhaseTensorPlot::set_layout()
{
  QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
  m_plot->xAxis->setTicker(logTicker);
  m_plot->xAxis->setLabel("T [s]");
  m_plot->yAxis->setLabel("Phase Tensor");
  m_plot->xAxis->setScaleType(QCPAxis::stLogarithmic);

  set_layout_generic({"XX", "XY", "YX", "YY"},
                     {"XX", "XY", "YX", "YY"});

  m_name2type = {{"XX", PTxx}, {"XY", PTxy}, {"YX", PTyx}, {"YY", PTyy}};
}
