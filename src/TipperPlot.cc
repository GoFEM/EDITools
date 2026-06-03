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

#include "include/TipperPlot.h"

namespace
{
constexpr double tipperArrowPixelsPerUnit = 90.0;
}

TipperPlot::TipperPlot(QCustomPlot *plot):
  MTDataPlot(plot), m_arrowMode(false)
{
  set_layout();
  update_legend();
  update_axis_style();
}

void TipperPlot::set_observed_data(MTStationData &data, bool rescaleAxes)
{
  m_data = &data;
  std::vector<dvector> tipper, tipper_err;
  m_data->get_tipper(tipper, tipper_err);
  const auto mask = m_data->tipper_mask();

  clear_arrow_items();
  if(m_arrowMode)
  {
    clear_graph_data();
    set_arrow_selection_data(mask, data.frequencies());
    draw_tipper_arrows(tipper, &mask, data.frequencies(), false);
    draw_reference_arrow();
  }
  else
  {
    set_graph_data(mask, tipper, tipper_err, data.frequencies());
  }

  apply_axis_ranges(rescaleAxes, m_arrowMode, QCPRange(-1, 1));

  m_plot->replot();
}

void TipperPlot::set_predicted_data(const MTStationData &data, bool rescaleAxes)
{
  std::vector<dvector> tipper, tipper_err;
  data.get_tipper(tipper, tipper_err);

  if(m_arrowMode)
    draw_tipper_arrows(tipper, nullptr, data.frequencies(), true);
  else
    set_graph_responses(tipper, data.frequencies());

  apply_axis_ranges(rescaleAxes, m_arrowMode, QCPRange(-1, 1));

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

void TipperPlot::set_arrow_mode(bool on)
{
  m_arrowMode = on;
  update_legend();
  update_axis_style();
}

bool TipperPlot::arrow_mode() const
{
  return m_arrowMode;
}

std::vector<RealDataType> TipperPlot::get_graph_data_types(const QCPGraph *graph) const
{
  if(m_arrowMode)
  {
    const int idx = graph->property("tipperArrowRole").toInt();
    if(idx == 1)
      return {RealTzx, RealTzy};
    if(idx == 2)
      return {ImagTzx, ImagTzy};
  }

  return MTDataPlot::get_graph_data_types(graph);
}

void TipperPlot::clear_arrow_items()
{
  for(auto *item: m_arrowItems)
    m_plot->removeItem(item);
  m_arrowItems.clear();

  for(auto *item: m_referenceItems)
    m_plot->removeItem(item);
  m_referenceItems.clear();
}

void TipperPlot::clear_graph_data()
{
  for(int i = 0; i < m_plot->graphCount(); ++i)
    m_plot->graph(i)->data()->clear();

  for(auto *errorBar: m_errorBars)
    errorBar->setData(QVector<double>());
}

void TipperPlot::update_legend()
{
  const std::vector<QString> pointNames = {"Re Tzx", "Re Tzy", "Im Tzx", "Im Tzy"};
  const auto &pointColors = PlotColors::tipperScalarColors();

  if(m_arrowMode)
  {
    for(int i = 0; i < m_plot->graphCount(); ++i)
      m_plot->graph(i)->removeFromLegend();

    if(m_plot->graphCount() >= 3)
    {
      m_plot->graph(0)->setName("Real");
      m_plot->graph(0)->setProperty("tipperArrowRole", 1);
      m_plot->graph(0)->setPen(QPen(PlotColors::tipperReal()));
      m_plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                        PlotColors::tipperReal(),
                                                        PlotColors::tipperReal(), 5));
      m_plot->graph(0)->addToLegend();

      m_plot->graph(2)->setName("Imaginary");
      m_plot->graph(2)->setProperty("tipperArrowRole", 2);
      m_plot->graph(2)->setPen(QPen(PlotColors::tipperImaginary()));
      m_plot->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                        PlotColors::tipperImaginary(),
                                                        PlotColors::tipperImaginary(), 5));
      m_plot->graph(2)->addToLegend();

      if(m_errorBars.size() >= 3)
      {
        m_errorBars[0]->setPen(QPen(PlotColors::tipperReal()));
        m_errorBars[2]->setPen(QPen(PlotColors::tipperImaginary()));
      }

      if(m_plot->graphCount() >= 7)
      {
        m_plot->graph(4)->setName("Real masked");
        m_plot->graph(4)->setProperty("tipperArrowRole", 1);
        m_plot->graph(4)->setPen(QPen(PlotColors::tipperMasked()));
        m_plot->graph(4)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                          PlotColors::tipperMasked(),
                                                          PlotColors::tipperMasked(), 5));

        m_plot->graph(6)->setName("Imaginary masked");
        m_plot->graph(6)->setProperty("tipperArrowRole", 2);
        m_plot->graph(6)->setPen(QPen(PlotColors::tipperMasked()));
        m_plot->graph(6)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                          PlotColors::tipperMasked(),
                                                          PlotColors::tipperMasked(), 5));

        if(m_errorBars.size() >= 7)
        {
          m_errorBars[4]->setPen(QPen(PlotColors::tipperMasked()));
          m_errorBars[6]->setPen(QPen(PlotColors::tipperMasked()));
        }
      }
    }
  }
  else
  {
    for(int i = 0; i < m_plot->graphCount(); ++i)
      m_plot->graph(i)->removeFromLegend();

    for(unsigned i = 0; i < pointNames.size() && i < static_cast<unsigned>(m_plot->graphCount()); ++i)
    {
      m_plot->graph(i)->setProperty("tipperArrowRole", 0);
      m_plot->graph(i)->setName(pointNames[i]);
      m_plot->graph(i)->setPen(QPen(pointColors[i]));
      m_plot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle,
                                                        pointColors[i], pointColors[i], 5));
      m_plot->graph(i)->addToLegend();

      if(i < m_errorBars.size())
        m_errorBars[i]->setPen(QPen(pointColors[i]));
    }

    for(unsigned i = 0; i < pointNames.size() && i + pointNames.size() < static_cast<unsigned>(m_plot->graphCount()); ++i)
    {
      m_plot->graph(i + pointNames.size())->setProperty("tipperArrowRole", 0);
      m_plot->graph(i + pointNames.size())->setName(pointNames[i]);
      m_plot->graph(i + pointNames.size())->setPen(QPen(PlotColors::masked()));
      m_plot->graph(i + pointNames.size())->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle,
                                                                            PlotColors::masked(), 5));

      if(i + pointNames.size() < m_errorBars.size())
        m_errorBars[i + pointNames.size()]->setPen(QPen(PlotColors::masked()));
    }

    for(unsigned i = 0; i < pointNames.size() && i + pointNames.size()*2 < static_cast<unsigned>(m_plot->graphCount()); ++i)
    {
      m_plot->graph(i + pointNames.size()*2)->setPen(QPen(pointColors[i]));
    }
  }
}

void TipperPlot::update_axis_style()
{
  m_plot->yAxis->setTicks(!m_arrowMode);
  m_plot->yAxis->setTickLabels(!m_arrowMode);
  m_plot->yAxis->setLabel(m_arrowMode ? "" : "Tipper");
}

void TipperPlot::set_arrow_selection_data(const std::vector<std::vector<bool>> &mask,
                                          const dvector &frequencies)
{
  if(m_plot->graphCount() < 7)
    return;

  QVector<double> realActiveX, realActiveY;
  QVector<double> realMaskedX, realMaskedY;
  QVector<double> imagActiveX, imagActiveY;
  QVector<double> imagMaskedX, imagMaskedY;

  for(unsigned i = 0; i < frequencies.size(); ++i)
  {
    const double period = 1.0 / frequencies[i];

    if(mask[0][i] && mask[1][i])
    {
      realActiveX.push_back(period);
      realActiveY.push_back(0.0);
    }
    else
    {
      realMaskedX.push_back(period);
      realMaskedY.push_back(0.0);
    }

    if(mask[2][i] && mask[3][i])
    {
      imagActiveX.push_back(period);
      imagActiveY.push_back(0.0);
    }
    else
    {
      imagMaskedX.push_back(period);
      imagMaskedY.push_back(0.0);
    }
  }

  m_plot->graph(0)->setData(realActiveX, realActiveY);
  m_plot->graph(4)->setData(realMaskedX, realMaskedY);
  m_plot->graph(2)->setData(imagActiveX, imagActiveY);
  m_plot->graph(6)->setData(imagMaskedX, imagMaskedY);
}

void TipperPlot::draw_arrow(double period, double xComponent,
                            double yComponent, const QPen &pen)
{
  QCPItemLine *arrow = new QCPItemLine(m_plot);
  arrow->setSelectable(false);
  arrow->setLayer("overlay");
  arrow->setPen(pen);
  arrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8));

  arrow->start->setType(QCPItemPosition::ptPlotCoords);
  arrow->start->setAxes(m_plot->xAxis, m_plot->yAxis);
  arrow->start->setCoords(period, 0.0);

  arrow->end->setParentAnchor(arrow->start);
  arrow->end->setType(QCPItemPosition::ptAbsolute);
  arrow->end->setCoords(xComponent * tipperArrowPixelsPerUnit,
                        -yComponent * tipperArrowPixelsPerUnit);

  m_arrowItems.push_back(arrow);
}

void TipperPlot::draw_reference_arrow()
{
  QCPItemLine *arrow = new QCPItemLine(m_plot);
  arrow->setSelectable(false);
  arrow->setLayer("overlay");
  arrow->setPen(QPen(PlotColors::tipperReal()));
  arrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8));
  arrow->start->setType(QCPItemPosition::ptAxisRectRatio);
  arrow->start->setAxisRect(m_plot->axisRect());
  arrow->start->setCoords(0.72, 0.12);
  arrow->end->setParentAnchor(arrow->start);
  arrow->end->setType(QCPItemPosition::ptAbsolute);
  arrow->end->setCoords(tipperArrowPixelsPerUnit, 0);
  m_referenceItems.push_back(arrow);

  QCPItemText *label = new QCPItemText(m_plot);
  label->setSelectable(false);
  label->setLayer("overlay");
  label->setPositionAlignment(Qt::AlignLeft | Qt::AlignVCenter);
  label->position->setParentAnchor(arrow->end);
  label->position->setType(QCPItemPosition::ptAbsolute);
  label->position->setCoords(8, 0);
  label->setText("|T| = 1.0");
  label->setColor(PlotColors::tipperReal());
  m_referenceItems.push_back(label);
}

void TipperPlot::draw_tipper_arrows(const std::vector<dvector> &tipper,
                                    const std::vector<std::vector<bool>> *mask,
                                    const dvector &frequencies,
                                    bool predicted)
{
  if(tipper.size() < 4)
    return;

  QPen realPen(PlotColors::tipperReal());
  QPen imagPen(PlotColors::tipperImaginary());
  if(predicted)
  {
    realPen.setStyle(Qt::DashLine);
    imagPen.setStyle(Qt::DashLine);
  }

  QPen maskedPen(PlotColors::tipperMasked());

  for(unsigned i = 0; i < frequencies.size(); ++i)
  {
    const double period = 1.0 / frequencies[i];

    QPen activeRealPen = realPen;
    if(mask != nullptr && (!(*mask)[0][i] || !(*mask)[1][i]))
      activeRealPen = maskedPen;
    draw_arrow(period, tipper[0][i], tipper[1][i], activeRealPen);

    QPen activeImagPen = imagPen;
    if(mask != nullptr && (!(*mask)[2][i] || !(*mask)[3][i]))
      activeImagPen = maskedPen;
    draw_arrow(period, tipper[2][i], tipper[3][i], activeImagPen);
  }
}
