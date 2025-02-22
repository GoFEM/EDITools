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

#include "TipperPlot.h"

TipperPlot::TipperPlot(QCustomPlot *plot):
  MTDataPlot(plot)
{
  set_layout();
}

void TipperPlot::set_observed_data(MTStationData &data, bool rescaleAxes)
{
  m_data = &data;
  std::vector<dvector> tipper, tipper_err;
  m_data->get_tipper(tipper, tipper_err);
  const auto mask = m_data->tipper_mask();

  set_graph_data(mask, tipper, tipper_err, data.frequencies());

  if(rescaleAxes)
  {
    m_plot->rescaleAxes();
    QCPRange xrange = m_plot->xAxis->range();
    m_plot->xAxis->setRange(xrange.lower / 2., xrange.upper * 2.);
  }

  m_plot->replot();
}

void TipperPlot::set_predicted_data(const MTStationData &data, bool rescaleAxes)
{
  std::vector<dvector> tipper, tipper_err;
  data.get_tipper(tipper, tipper_err);

  set_graph_responses(tipper, data.frequencies());

  if(rescaleAxes)
  {
    m_plot->rescaleAxes();
    QCPRange xrange = m_plot->xAxis->range();
    m_plot->xAxis->setRange(xrange.lower / 2., xrange.upper * 2.);
  }

  m_plot->replot();
}

void TipperPlot::set_layout()
{
  QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
  m_plot->xAxis->setTicker(logTicker);
  m_plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
  m_plot->xAxis->setLabel("T [s]");
  m_plot->yAxis->setLabel("Tipper");

  set_layout_generic({"Re Tzx", "Re Tzy", "Im Tzx", "Im Tzy"},
                     {"Re Tzx", "Re Tzy", "Im Tzx", "Im Tzy"});

  m_name2type = {{"Re Tzx", RealTzx}, {"Re Tzy", RealTzy},
                 {"Im Tzx", ImagTzx}, {"Im Tzy", ImagTzy}};
}
