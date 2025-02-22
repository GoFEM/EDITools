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

#include "include/MapPlot.h"

MapPlot::MapPlot(QCustomPlot *plot, QWidget* parent):
  m_plot(plot), m_parent(parent)
{
  set_layout();
}

void MapPlot::set_data(const std::vector<std::array<double, 3>> &locations)
{
  QVector<double> x(locations.size()),
                  y(locations.size());

  for(unsigned i = 0; i < locations.size(); ++i)
  {
    x[i] = locations[i][1];
    y[i] = locations[i][0];
  }

  m_plot->graph(0)->setData(x, y);
  m_plot->rescaleAxes();
  m_plot->replot();
}

void MapPlot::set_selected_points(const std::vector<std::array<double, 3>> &locations)
{
  QVector<double> x(locations.size()),
                  y(locations.size());

  for(unsigned i = 0; i < locations.size(); ++i)
  {
    x[i] = locations[i][1];
    y[i] = locations[i][0];
  }

  m_plot->graph(1)->setData(x, y);
  m_plot->rescaleAxes();
  m_plot->replot();
}

void MapPlot::get_point_value(const unsigned idx, double &key, double &value) const
{
  key = m_plot->graph(0)->data()->at(idx)->key;
  value = m_plot->graph(0)->data()->at(idx)->value;
}

void MapPlot::set_layout()
{
  m_plot->xAxis->setLabel("Longitude");
  m_plot->yAxis->setLabel("Latitude");

  QCPGraph* graph = m_plot->addGraph();
  graph->setName("Stations");
  graph->setLineStyle(QCPGraph::lsNone);
  graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::black, 5));
  graph->setSelectable(QCP::stSingleData);

  m_parent->connect(graph, SIGNAL(selectionChanged(const QCPDataSelection &)),
                    m_parent, SLOT(stationSelected(const QCPDataSelection &)));

  graph = m_plot->addGraph();
  graph->setName("Selected");
  graph->setLineStyle(QCPGraph::lsNone);
  graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::red, 7));

  m_plot->setNoAntialiasingOnDrag(true);
  m_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

  m_plot->legend->setVisible(true);
  m_plot->legend->setBrush(QBrush(QColor(255,255,255,100)));
  m_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
}
