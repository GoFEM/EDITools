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
    if(std::any_of(m_componentVisible.begin(), m_componentVisible.end(), [](bool visible) { return visible; }))
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

void TipperPlot::clear_observed_data()
{
  clear_arrow_items();
  MTDataPlot::clear_observed_data();
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

void TipperPlot::set_view_style(const TFViewStyle &style)
{
  const double scale = style.arrowScale / m_viewStyle.arrowScale;
  MTDataPlot::set_view_style(style);
  if(!m_arrowMode) return;
  for(int i: {0, 2, 4, 6}) {
    const QColor color = view_color(i >= 4 ? 7 : i == 0 ? 5 : 6);
    auto *graph = m_plot->graph(i);
    auto scatter = graph->scatterStyle(); scatter.setPen(QPen(color)); scatter.setBrush(color);
    graph->setScatterStyle(scatter); graph->setPen(QPen(color, style.lineWidth));
    m_errorBars[i]->setPen(QPen(color, style.lineWidth));
  }
  const auto restyle = [&](QCPItemLine *arrow) {
    auto pen = arrow->pen(); pen.setColor(view_color(arrow->property("viewColor").toUInt()));
    pen.setWidthF(style.lineWidth); arrow->setPen(pen);
    if(arrow->end->type() == QCPItemPosition::ptAxisRectRatio)
      arrow->start->setCoords(-style.arrowScale, 0);
    else arrow->end->setCoords(arrow->end->coords() * scale);
  };
  for(auto *arrow: m_arrowItems) restyle(arrow);
  for(auto *item: m_referenceItems) {
    if(auto *arrow = qobject_cast<QCPItemLine *>(item)) restyle(arrow);
    if(auto *label = qobject_cast<QCPItemText *>(item)) label->setColor(view_color(5));
  }
  update_component_legend();
}

std::vector<RealDataType> TipperPlot::get_graph_data_types(const QCPGraph *graph) const
{
  if(m_arrowMode)
  {
    const int idx = graph->property("tipperArrowRole").toInt();
    if(idx == 1 || idx == 2) {
      const unsigned offset = idx == 1 ? 0 : 2;
      const std::array<RealDataType, 4> types{{RealTzx, RealTzy, ImagTzx, ImagTzy}};
      std::vector<RealDataType> result;
      for(unsigned c = offset; c < offset + 2; ++c)
        if(m_componentVisible[c]) result.push_back(types[c]);
      return result;
    }
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

void TipperPlot::apply_component_visibility()
{
  if(!m_arrowMode) {
    MTDataPlot::apply_component_visibility();
    return;
  }
  for(unsigned group = 0; group < 2; ++group) {
    unsigned active = 0;
    for(unsigned c = group * 2; c < group * 2 + 2; ++c)
      if(m_componentVisible[c]) active |= 1u << c;
    if(active) m_arrowComponents[group] = active;
  }
  for(int i = 0; i < m_plot->graphCount(); ++i) {
    auto *graph = m_plot->graph(i);
    const int role = graph->property("tipperArrowRole").toInt();
    const unsigned offset = role == 1 ? 0 : 2;
    const bool visible = (role == 1 || role == 2) &&
                        (m_componentVisible[offset] || m_componentVisible[offset + 1]);
    graph->setVisible(visible);
    graph->setSelectable(visible && i < static_cast<int>(m_errorBars.size())
                         ? QCP::stMultipleDataRanges : QCP::stNone);
    if(!visible) graph->setSelection(QCPDataSelection());
    if(i < 4 && (role == 1 || role == 2)) {
      QString name = role == 1 ? tr("Real") : tr("Imaginary");
      const auto mask = m_arrowComponents[offset / 2];
      if(!(mask & (1u << offset))) name += " Tzy";
      else if(!(mask & (1u << (offset + 1)))) name += " Tzx";
      graph->setName(name);
    }
  }
  for(auto *errorBar: m_errorBars) errorBar->setVisible(false);
  update_component_legend();
}

void TipperPlot::set_legend_component_visible(unsigned component, bool visible)
{
  if(!m_arrowMode) {
    MTDataPlot::set_legend_component_visible(component, visible);
    return;
  }
  auto components = m_componentVisible;
  const auto mask = m_arrowComponents[component / 2];
  for(unsigned c = component; c < component + 2; ++c)
    components[c] = visible && (mask & (1u << c));
  set_component_visibility(components);
  emit componentVisibilityChanged();
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
    if(m_plot->graphCount() >= 3)
    {
      m_plot->graph(0)->setName("Real");
      m_plot->graph(0)->setProperty("tipperArrowRole", 1);
      m_plot->graph(0)->setPen(QPen(PlotColors::tipperReal()));
      m_plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                        PlotColors::tipperReal(),
                                                        PlotColors::tipperReal(), 5));

      m_plot->graph(2)->setName("Imaginary");
      m_plot->graph(2)->setProperty("tipperArrowRole", 2);
      m_plot->graph(2)->setPen(QPen(PlotColors::tipperImaginary()));
      m_plot->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc,
                                                        PlotColors::tipperImaginary(),
                                                        PlotColors::tipperImaginary(), 5));

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
    for(unsigned i = 0; i < pointNames.size() && i < static_cast<unsigned>(m_plot->graphCount()); ++i)
    {
      m_plot->graph(i)->setProperty("tipperArrowRole", 0);
      m_plot->graph(i)->setName(pointNames[i]);
      m_plot->graph(i)->setPen(QPen(pointColors[i]));
      m_plot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle,
                                                        pointColors[i], pointColors[i], 5));

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
  rebuild_component_legend(m_arrowMode ? std::vector<int>{0, 2} : std::vector<int>{0, 1, 2, 3});
  apply_component_visibility();
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

    if((!m_componentVisible[0] || mask[0][i]) && (!m_componentVisible[1] || mask[1][i]))
    {
      realActiveX.push_back(period);
      realActiveY.push_back(0.0);
    }
    else
    {
      realMaskedX.push_back(period);
      realMaskedY.push_back(0.0);
    }

    if((!m_componentVisible[2] || mask[2][i]) && (!m_componentVisible[3] || mask[3][i]))
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
                            double yComponent, const QPen &pen, unsigned colorIndex)
{
  if(!std::isfinite(period) || !std::isfinite(xComponent) || !std::isfinite(yComponent))
    return;

  QCPItemLine *arrow = new QCPItemLine(m_plot);
  arrow->setSelectable(false);
  arrow->setLayer("overlay");
  arrow->setPen(pen);
  arrow->setProperty("viewColor", colorIndex);
  arrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8));

  arrow->start->setType(QCPItemPosition::ptPlotCoords);
  arrow->start->setAxes(m_plot->xAxis, m_plot->yAxis);
  arrow->start->setCoords(period, 0.0);

  arrow->end->setParentAnchor(arrow->start);
  arrow->end->setType(QCPItemPosition::ptAbsolute);
  arrow->end->setCoords(xComponent * m_viewStyle.arrowScale,
                        -yComponent * m_viewStyle.arrowScale);

  m_arrowItems.push_back(arrow);
}

void TipperPlot::draw_reference_arrow()
{
  QCPItemLine *arrow = new QCPItemLine(m_plot);
  arrow->setObjectName("tipperReferenceArrow");
  arrow->setSelectable(false);
  arrow->setLayer("overlay");
  arrow->setPen(QPen(view_color(5), m_viewStyle.lineWidth));
  arrow->setProperty("viewColor", 5u);
  arrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8));
  // Anchor the head at the right edge so resizing or a larger scale cannot
  // push the reference and its label outside the plot.
  arrow->end->setType(QCPItemPosition::ptAxisRectRatio);
  arrow->end->setAxisRect(m_plot->axisRect());
  arrow->end->setCoords(.96, .12);
  arrow->start->setParentAnchor(arrow->end);
  arrow->start->setType(QCPItemPosition::ptAbsolute);
  arrow->start->setCoords(-m_viewStyle.arrowScale, 0);
  m_referenceItems.push_back(arrow);

  QCPItemText *label = new QCPItemText(m_plot);
  label->setSelectable(false);
  label->setLayer("overlay");
  label->setPositionAlignment(Qt::AlignRight | Qt::AlignTop);
  label->position->setParentAnchor(arrow->end);
  label->position->setType(QCPItemPosition::ptAbsolute);
  label->position->setCoords(0, 5);
  label->setText("|T| = 1.0");
  label->setColor(view_color(5));
  m_referenceItems.push_back(label);
}

void TipperPlot::draw_tipper_arrows(const std::vector<dvector> &tipper,
                                    const std::vector<std::vector<bool>> *mask,
                                    const dvector &frequencies,
                                    bool predicted)
{
  if(tipper.size() < 4)
    return;

  QPen realPen(view_color(5), m_viewStyle.lineWidth);
  QPen imagPen(view_color(6), m_viewStyle.lineWidth);
  if(predicted)
  {
    realPen.setStyle(Qt::DashLine);
    imagPen.setStyle(Qt::DashLine);
  }

  QPen maskedPen(view_color(7), m_viewStyle.lineWidth);

  for(unsigned i = 0; i < frequencies.size(); ++i)
  {
    const double period = 1.0 / frequencies[i];

    for(unsigned offset: {0u, 2u}) {
      if(!m_componentVisible[offset] && !m_componentVisible[offset + 1]) continue;
      QPen pen = offset == 0 ? realPen : imagPen;
      unsigned colorIndex = offset == 0 ? 5 : 6;
      if(mask && ((m_componentVisible[offset] && !(*mask)[offset][i]) ||
                  (m_componentVisible[offset + 1] && !(*mask)[offset + 1][i])))
        { pen = maskedPen; colorIndex = 7; }
      draw_arrow(period, m_componentVisible[offset] ? tipper[offset][i] : 0.,
                 m_componentVisible[offset + 1] ? tipper[offset + 1][i] : 0., pen, colorIndex);
    }
  }
}
