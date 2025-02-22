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

#include "MTDataPlot.h"

MTDataPlot::MTDataPlot(QCustomPlot *plot):
  m_plot(plot), m_associated_plot(nullptr)
{
  m_contextMenu = new QMenu;
  m_contextMenu->addAction(tr("Mask"), this, &MTDataPlot::maskSelectedData);
  m_contextMenu->addAction(tr("Inverted Mask"), this, &MTDataPlot::invMaskSelectedData);
  m_contextMenu->addAction(tr("Unmask"), this, &MTDataPlot::unmaskSelectedData);

  m_plot->setFocusPolicy(Qt::StrongFocus);

  connect(m_plot, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*)));
}

void MTDataPlot::set_associated_plot(MTDataPlot &plot)
{
  m_associated_plot = &plot;
}

void MTDataPlot::set_error_bars_visible(bool on)
{
  for(auto &eb: m_errorBars)
    eb->setVisible(on);
}

void MTDataPlot::set_masking_mode(bool on)
{
  if(on)
    m_plot->setSelectionRectMode(QCP::srmSelect);
  else
    m_plot->setSelectionRectMode(QCP::srmNone);
}

void MTDataPlot::dataSelected(bool selected)
{
  if(!selected)
    return;

  m_plot->replot();

//  QCPDataSelection selection = graph->selection();
//  std::cout << graph->name().toStdString() << "\t" << selection.dataPointCount() << std::endl;
}

void MTDataPlot::plotContextRequest(QPoint pos)
{
  m_contextMenu->popup(m_plot->mapToGlobal(pos));
}

void MTDataPlot::showPointToolTip(QMouseEvent *event)
{
  double x = m_plot->xAxis->pixelToCoord(event->pos().x());
  double y = m_plot->yAxis->pixelToCoord(event->pos().y());

  m_plot->setToolTip(QString("X: %1\nY: %2").arg(x).arg(y));
}

void MTDataPlot::set_graph_data(const std::vector<std::vector<bool> > &mask,
                                const std::vector<std::vector<double> > &data,
                                const std::vector<std::vector<double> > &data_err,
                                const std::vector<double> &frequencies)
{
  // Active data
  {
    for(unsigned i = 0; i < data.size(); ++i)
    {
      QVector<double> x, y, yerr;
      for (unsigned j = 0; j < frequencies.size(); ++j)
      {
        if(!mask[i][j])
          continue;

        x.push_back(1.0 / frequencies[j]);
        y.push_back(data[i][j]);
        yerr.push_back(data_err[i][j]);
      }

      m_plot->graph(i)->setData(x, y);
      m_errorBars[i]->setData(yerr);
    }
  }

  // Masked data
  {
    for(unsigned i = 0; i < data.size(); ++i)
    {
      QVector<double> x, y, yerr;
      for (unsigned j = 0; j < frequencies.size(); ++j)
      {
        if(mask[i][j])
          continue;

        x.push_back(1.0 / frequencies[j]);
        y.push_back(data[i][j]);
        yerr.push_back(data_err[i][j]);
      }

      m_plot->graph(i + data.size())->setData(x, y);
      m_errorBars[i + data.size()]->setData(yerr);
    }
  }
}

void MTDataPlot::set_graph_responses(const std::vector<std::vector<double> > &data,
                                     const std::vector<double> &frequencies)
{
  for(unsigned i = 0; i < data.size(); ++i)
  {
    QVector<double> x, y;
    for (unsigned j = 0; j < frequencies.size(); ++j)
    {
      x.push_back(1.0 / frequencies[j]);
      y.push_back(data[i][j]);
    }

    m_plot->graph(i + data.size()*2)->setData(x, y);
  }
}

void MTDataPlot::set_layout_generic(const std::vector<QString> &data_graph_names,
                                    const std::vector<QString> &masked_graph_names)
{
  std::vector<QColor> colors = {Qt::cyan, Qt::red, Qt::blue, Qt::green};

  for(unsigned i = 0; i < data_graph_names.size(); ++i)
  {
    QCPGraph* graph = m_plot->addGraph();
    graph->setName(data_graph_names[i]);
    graph->setLineStyle(QCPGraph::lsNone);
    graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, colors[i], colors[i], 5));
    graph->setSelectable(QCP::stMultipleDataRanges);

    QCPScatterStyle style = graph->selectionDecorator()->scatterStyle();
    style.setSize(9);
    graph->selectionDecorator()->setScatterStyle(style, QCPScatterStyle::spSize);
    graph->selectionDecorator()->setUsedScatterProperties(QCPScatterStyle::spSize);

    connect(graph, SIGNAL(selectionChanged(bool)),
            this, SLOT(dataSelected(bool)));

    m_errorBars.push_back(new QCPErrorBars(m_plot->xAxis, m_plot->yAxis));
    m_errorBars.back()->removeFromLegend();
    m_errorBars.back()->setPen(QPen(colors[i]));
    m_errorBars.back()->setSelectable(QCP::stNone);
    m_errorBars.back()->setDataPlottable(graph);
  }

  for(unsigned i = 0; i < masked_graph_names.size(); ++i)
  {
    QCPGraph* graph = m_plot->addGraph();
    graph->setName(masked_graph_names[i]);
    graph->setLineStyle(QCPGraph::lsNone);
    graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::gray, 5));
    graph->setSelectable(QCP::stMultipleDataRanges);
    graph->selectionDecorator()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::gray, 8));
    graph->removeFromLegend();

    QCPScatterStyle style = graph->selectionDecorator()->scatterStyle();
    style.setSize(9);
    graph->selectionDecorator()->setScatterStyle(style, QCPScatterStyle::spSize);
    graph->selectionDecorator()->setUsedScatterProperties(QCPScatterStyle::spSize);

    connect(graph, SIGNAL(selectionChanged(bool)), this, SLOT(dataSelected(bool)));

    m_errorBars.push_back(new QCPErrorBars(m_plot->xAxis, m_plot->yAxis));
    m_errorBars.back()->removeFromLegend();
    m_errorBars.back()->setPen(QPen(Qt::gray));
    m_errorBars.back()->setSelectable(QCP::stNone);
    m_errorBars.back()->setDataPlottable(graph);
  }

  for(unsigned i = 0; i < data_graph_names.size(); ++i)
  {
    QCPGraph* graph = m_plot->addGraph();
    graph->setName(data_graph_names[i] + " Predicted");
    graph->setLineStyle(QCPGraph::lsLine);
    graph->setPen(colors[i]);
    graph->setSelectable(QCP::stNone);
    graph->removeFromLegend();
  }

  m_plot->setNoAntialiasingOnDrag(true);
  m_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iMultiSelect);
  m_plot->setSelectionRectMode(QCP::srmSelect);

  m_plot->legend->setVisible(true);
  m_plot->legend->setBrush(QBrush(QColor(255,255,255,100)));
  m_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop); //

  m_plot->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(m_plot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(plotContextRequest(QPoint)));
}

void MTDataPlot::mask_selected_data(bool on)
{
  QList<QCPGraph*> selectedGraphs = m_plot->selectedGraphs();
  for(auto &graph: selectedGraphs)
  {
    QCPDataSelection selection = graph->selection();

    auto itype = m_name2type.find(graph->name().toStdString());
    const RealDataType type = itype->second;

    for (QCPDataRange dataRange: selection.dataRanges())
    {
      auto begin = graph->data()->at(dataRange.begin());
      auto end = graph->data()->at(dataRange.end());
      for (auto it = begin; it != end; ++it)
      {
        m_data->set_data_mask(type, 1. / it->key, on);
      }
    }
  }

  set_observed_data(*m_data, false);
  if(m_associated_plot)
    m_associated_plot->set_observed_data(*m_associated_plot->m_data, false);
}

void MTDataPlot::maskSelectedData()
{
  mask_selected_data(false);
}

void MTDataPlot::invMaskSelectedData()
{
  auto frequencies = m_data->frequencies();

  QList<QCPGraph*> selectedGraphs = m_plot->selectedGraphs();
  for(auto &graph: selectedGraphs)
  {
    QCPDataSelection selection = graph->selection();

    auto itype = m_name2type.find(graph->name().toStdString());
    const RealDataType type = itype->second;

    dvector selected_freqs;

    for (QCPDataRange dataRange: selection.dataRanges())
    {
      auto begin = graph->data()->at(dataRange.begin());
      auto end = graph->data()->at(dataRange.end());
      for (auto it = begin; it != end; ++it)
        selected_freqs.push_back(1. / it->key);
    }

    for(unsigned i = 0; i < frequencies.size(); ++i)
    {
      bool found = false;
      for(unsigned j = 0; j < selected_freqs.size(); ++j)
        if(fabs(selected_freqs[j] - frequencies[i]) / frequencies[i] < 1e-3)
        {
          found = true;
          break;
        }

      m_data->set_data_mask(type, frequencies[i], found);
    }
  }

  set_observed_data(*m_data, false);

  if(m_associated_plot)
    m_associated_plot->set_observed_data(*m_associated_plot->m_data, false);
}

void MTDataPlot::unmaskSelectedData()
{
  mask_selected_data(true);
}
