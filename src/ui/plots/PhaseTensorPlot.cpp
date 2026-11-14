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

#include "ui/plots/PhaseTensorPlot.h"
#include "core/MTMapData.h"
#include "ui/plots/PlotColorMaps.h"

class PhaseTensorPlot::Ellipses : public QCPLayerable
{
public:
  explicit Ellipses(QCustomPlot *plot): QCPLayerable(plot) { setLayer("main"); }
  struct Glyph { double period; MTMapData::PhaseTensor tensor; bool masked; };
  std::vector<Glyph> observed, predicted;
  QCPColorGradient colors = PlotColorMaps::gradient();
  TFViewStyle style;
protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(true); }
  void draw(QCPPainter *p) override {
    const double radius = std::min(style.ellipseSize / 2., clipRect().height() * .3);
    for(const auto *samples: {&observed, &predicted}) for(const auto &s: *samples) {
      p->save();
      p->translate(mParentPlot->xAxis->coordToPixel(s.period), mParentPlot->yAxis->coordToPixel(0.));
      p->rotate(s.tensor.azimuth);
      const auto color = s.masked ? QColor(QString::fromStdString(style.colors[4])) :
        QColor::fromRgb(colors.color(s.tensor.phiMin, {style.colorMin, style.colorMax}));
      p->setPen(QPen(samples == &predicted ? QColor("#222222") : color,
                    style.lineWidth * (samples == &predicted ? 1.2 : .7)));
      p->setBrush(samples == &predicted ? QBrush(Qt::NoBrush) : QBrush(color));
      p->drawEllipse(QRectF(-radius*s.tensor.axisRatio, -radius, 2*radius*s.tensor.axisRatio, 2*radius));
      p->restore();
    }
  }
};

PhaseTensorPlot::PhaseTensorPlot(QCustomPlot *plot):
  MTDataPlot(plot)
{
  set_layout();
  m_ellipses = new Ellipses(plot);
  m_ellipses->setVisible(false);
  m_colorScale = new QCPColorScale(plot);
  m_colorScale->setType(QCPAxis::atRight);
  m_colorScale->setGradient(m_ellipses->colors);
  m_colorScale->setDataRange({0, 90});
  m_colorScale->setLabel(tr("φmin (°)"));
  m_colorScale->setBarWidth(10);
  m_colorScale->setVisible(false);
  // Keep the scalar layout unchanged until ellipse mode is selected.
}

void PhaseTensorPlot::set_observed_data(MTStationData &data, bool rescaleAxes)
{
  m_data = &data;

  if(m_ellipseMode) {
    set_ellipses(data, false);
    apply_axis_ranges(rescaleAxes);
    m_plot->yAxis->setRange(-1, 1);
    m_plot->replot();
    return;
  }

  std::vector<dvector> phase, phase_err;
  m_data->get_phase_tensor(phase, phase_err);
  const auto mask = m_data->phase_tensor_mask();

  set_graph_data(mask, phase, phase_err, data.frequencies());

  apply_axis_ranges(rescaleAxes);

  m_plot->replot();
}

void PhaseTensorPlot::set_predicted_data(const MTStationData &data, bool rescaleAxes)
{
  if(m_ellipseMode) { set_ellipses(data, true); m_plot->replot(); return; }
  std::vector<dvector> phase, phase_err;
  data.get_phase_tensor(phase, phase_err);

  set_graph_responses(phase, data.frequencies());

  apply_axis_ranges(rescaleAxes);

  m_plot->replot();
}

void PhaseTensorPlot::set_ellipse_mode(bool on)
{
  if(on == m_ellipseMode) return;
  m_ellipseMode = on;
  m_colorScale->setVisible(on);
  if(on) m_plot->plotLayout()->addElement(0, 1, m_colorScale);
  else { m_plot->plotLayout()->take(m_colorScale); m_plot->plotLayout()->simplify(); }
  m_plot->yAxis->setLabel(on ? tr("Phase tensor") : tr("Phase Tensor"));
  m_plot->yAxis->setTicks(!on);
  m_plot->yAxis->setTickLabels(!on);
  m_plot->yAxis->grid()->setVisible(!on);
  m_plot->graph(0)->setName(on ? tr("PT") : QStringLiteral("XX"));
  m_plot->graph(4)->setName(on ? tr("PT") : QStringLiteral("XX"));
  m_plot->setToolTip(on ? tr("Equal major diameters; color is minimum principal phase. Configure size, colors and color limits in the View ribbon. North is up, east right. Outline: calculated response. Select ellipse centers to mask/unmask the whole tensor; Link Z / PT controls impedance linkage.") : QString());
  rebuild_component_legend(on ? std::vector<int>{0} : std::vector<int>{0, 1, 2, 3});
  apply_component_visibility();
}

void PhaseTensorPlot::set_view_style(const TFViewStyle &style)
{
  MTDataPlot::set_view_style(style);
  m_ellipses->style = style;
  auto gradient = PlotColorMaps::gradient(QString::fromStdString(style.colorMap));
  if(style.reverseColorMap) gradient = gradient.inverted();
  m_ellipses->colors = gradient;
  m_colorScale->setGradient(gradient);
  m_colorScale->setDataRange({style.colorMin, style.colorMax});
}

void PhaseTensorPlot::apply_component_visibility()
{
  MTDataPlot::apply_component_visibility();
  if(!m_ellipses) return;
  const bool visible = std::any_of(m_componentVisible.begin(), m_componentVisible.end(), [](bool v) { return v; });
  m_ellipses->setVisible(m_ellipseMode && visible);
  if(!m_ellipseMode) return;
  for(int i = 0; i < m_plot->graphCount(); ++i) {
    auto *graph = m_plot->graph(i);
    graph->setVisible(visible && (i == 0 || i == 4));
    graph->setSelectable(graph->visible() ? QCP::stMultipleDataRanges : QCP::stNone);
    if(!graph->visible()) graph->setSelection(QCPDataSelection());
  }
  for(auto *bars: m_errorBars) bars->setVisible(false);
  update_component_legend();
}

void PhaseTensorPlot::set_legend_component_visible(unsigned component, bool visible)
{
  if(!m_ellipseMode) { MTDataPlot::set_legend_component_visible(component, visible); return; }
  set_component_visibility({visible, visible, visible, visible});
  emit componentVisibilityChanged();
}

std::vector<RealDataType> PhaseTensorPlot::get_graph_data_types(const QCPGraph *graph) const
{
  return m_ellipseMode ? std::vector<RealDataType>{PTxx, PTxy, PTyx, PTyy}
                       : MTDataPlot::get_graph_data_types(graph);
}

void PhaseTensorPlot::set_ellipses(const MTStationData &data, bool predicted)
{
  auto &samples = predicted ? m_ellipses->predicted : m_ellipses->observed;
  samples.clear();
  std::vector<dvector> tensor, errors;
  data.get_phase_tensor(tensor, errors);
  const auto mask = data.phase_tensor_mask();
  QVector<double> active, masked;
  for(unsigned i = 0; i < data.frequencies().size(); ++i) {
    MTMapData::PhaseTensor result;
    if(!MTMapData::phase_tensor({tensor[0][i], tensor[1][i], tensor[2][i], tensor[3][i]}, result)) continue;
    const bool enabled = predicted || (mask[0][i] && mask[1][i] && mask[2][i] && mask[3][i]);
    const double period = 1. / data.frequencies()[i];
    if(!std::isfinite(period) || period <= 0.) continue;
    samples.push_back({period, result, !enabled});
    (enabled ? active : masked).push_back(period);
  }
  if(!predicted) {
    for(int i = 0; i < m_plot->graphCount(); ++i) m_plot->graph(i)->data()->clear();
    m_plot->graph(0)->setData(active, QVector<double>(active.size(), 0.));
    m_plot->graph(4)->setData(masked, QVector<double>(masked.size(), 0.));
    apply_component_visibility();
  }
}

void PhaseTensorPlot::clear_predicted_data()
{
  m_ellipses->predicted.clear();
  MTDataPlot::clear_predicted_data();
}

void PhaseTensorPlot::clear_observed_data()
{
  m_ellipses->observed.clear(); m_ellipses->predicted.clear();
  MTDataPlot::clear_observed_data();
}

void PhaseTensorPlot::showPointToolTip(QMouseEvent *event)
{
  if(!m_ellipseMode) { MTDataPlot::showPointToolTip(event); return; }
  QString text = tr("North ↑ · East → · color: φmin · outline: response\nSelect centers to mask/unmask all PT entries.");
  double nearest = 18.;
  for(const auto &s: m_ellipses->observed) {
    const auto center = QPointF(m_plot->xAxis->coordToPixel(s.period), m_plot->yAxis->coordToPixel(0.));
    const double distance = QLineF(center, event->position()).length();
    if(distance >= nearest) continue;
    nearest = distance;
    text = tr("T: %1 s · φmin: %2° · φmax: %3°\nAzimuth: %4° · skew: %5°%6\nSelect center to mask/unmask all PT entries.")
      .arg(s.period, 0, 'g', 5).arg(s.tensor.phiMin, 0, 'f', 1).arg(s.tensor.phiMax, 0, 'f', 1)
      .arg(s.tensor.azimuth, 0, 'f', 1).arg(s.tensor.skew, 0, 'f', 1).arg(s.masked ? tr(" · masked") : QString());
  }
  m_plot->setToolTip(text);
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
