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

#include "include/PhasePlot.h"

namespace
{
const double minPhaseApparentResistivity = 1e-9;
}

PhasePlot::PhasePlot(QCustomPlot *plot):
  MTDataPlot(plot), m_phaseWrap(true)
{
  set_layout();
}

void PhasePlot::set_observed_data(MTStationData &data, bool rescaleAxes)
{
  m_data = &data;

  std::vector<dvector> phase, phase_err;
  m_data->get_phase(phase, phase_err);
  std::vector<dvector> appRes, appRes_err;
  m_data->get_apparent_resistivity(appRes, appRes_err);
  wrap_yx_yy_to_first_quadrant(phase);
  const auto mask = m_data->impedance_mask();

  set_phase_graph_data(mask, phase, phase_err, appRes, data.frequencies());

  apply_axis_ranges(rescaleAxes, true, QCPRange(0, 90));
  m_plot->replot();
}

void PhasePlot::set_predicted_data(const MTStationData &data, bool rescaleAxes)
{
  std::vector<dvector> phase, phase_err;
  data.get_phase(phase, phase_err);
  std::vector<dvector> appRes, appRes_err;
  data.get_apparent_resistivity(appRes, appRes_err);
  wrap_yx_yy_to_first_quadrant(phase);

  set_phase_graph_responses(phase, appRes, data.frequencies());

  apply_axis_ranges(rescaleAxes, true, QCPRange(0, 90));
  m_plot->replot();
}

void PhasePlot::set_layout()
{
  QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
  m_plot->xAxis->setTicker(logTicker);
  m_plot->xAxis->setLabel("T [s]");
  m_plot->yAxis->setLabel(QString::fromWCharArray(L"Phase [\u00b0]"));
  m_plot->xAxis->setScaleType(QCPAxis::stLogarithmic);

  set_layout_generic({"XX", "XY", "YX", "YY"},
                     {"XX", "XY", "YX", "YY"});

  m_name2type = {{"XX", PhsZxx}, {"XY", PhsZxy}, {"YX", PhsZyx}, {"YY", PhsZyy}};
}

void PhasePlot::set_phase_wrap(bool on)
{
  m_phaseWrap = on;
}

bool PhasePlot::phase_wrap() const
{
  return m_phaseWrap;
}

void PhasePlot::wrap_yx_yy_to_first_quadrant(std::vector<dvector> &phase) const
{
  if(!m_phaseWrap || phase.size() < 4)
    return;

  for(unsigned component: {2u, 3u})
  {
    for(auto &value: phase[component])
    {
      if(value >= -180.0 && value <= -90.0)
        value += 180.0;
      else if(value >= 180.0 && value <= 270.0)
        value -= 180.0;
    }
  }
}

void PhasePlot::set_phase_graph_data(const std::vector<std::vector<bool>> &mask,
                                     const std::vector<dvector> &phase,
                                     const std::vector<dvector> &phase_err,
                                     const std::vector<dvector> &appRes,
                                     const dvector &frequencies)
{
  for(unsigned i = 0; i < phase.size(); ++i)
  {
    QVector<double> x, y, yerr;
    for(unsigned j = 0; j < frequencies.size(); ++j)
    {
      if(!mask[i][j] || appRes[i][j] <= minPhaseApparentResistivity)
        continue;

      x.push_back(1.0 / frequencies[j]);
      y.push_back(phase[i][j]);
      yerr.push_back(phase_err[i][j]);
    }

    m_plot->graph(i)->setData(x, y);
    m_errorBars[i]->setData(yerr);
  }

  for(unsigned i = 0; i < phase.size(); ++i)
  {
    QVector<double> x, y, yerr;
    for(unsigned j = 0; j < frequencies.size(); ++j)
    {
      if(mask[i][j] || appRes[i][j] <= minPhaseApparentResistivity)
        continue;

      x.push_back(1.0 / frequencies[j]);
      y.push_back(phase[i][j]);
      yerr.push_back(phase_err[i][j]);
    }

    m_plot->graph(i + phase.size())->setData(x, y);
    m_errorBars[i + phase.size()]->setData(yerr);
  }
}

void PhasePlot::set_phase_graph_responses(const std::vector<dvector> &phase,
                                          const std::vector<dvector> &appRes,
                                          const dvector &frequencies)
{
  for(unsigned i = 0; i < phase.size(); ++i)
  {
    QVector<double> x, y;
    for(unsigned j = 0; j < frequencies.size(); ++j)
    {
      if(appRes[i][j] <= minPhaseApparentResistivity)
        continue;

      x.push_back(1.0 / frequencies[j]);
      y.push_back(phase[i][j]);
    }

    m_plot->graph(i + phase.size()*2)->setData(x, y);
  }
}
