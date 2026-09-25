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

#include "include/MTDataPlot.h"
#include <QCheckBox>
#include <QSignalBlocker>
#include <QStyleOptionButton>
#include <functional>

namespace {
// Native checkboxes handle mouse, keyboard and accessibility on screen. The
// legend item draws the same controls when QCustomPlot exports without widgets.
class ComponentLegendItem : public QCPPlottableLegendItem
{
public:
  ComponentLegendItem(QCPLegend *legend, QCPGraph *graph, int component,
                      const std::function<void(bool)> &toggle)
    : QCPPlottableLegendItem(legend, graph), m_checkBox(new QCheckBox(legend->parentPlot()))
  {
    m_checkBox->setObjectName("componentVisibility" + QString::number(component));
    m_checkBox->setFont(font());
    setSelectable(false);
    sync();
    QObject::connect(m_checkBox, &QCheckBox::toggled, m_checkBox, toggle);
    QObject::connect(legend->parentPlot(), &QCustomPlot::afterReplot, m_checkBox, [this] {
      m_checkBox->setGeometry(rect());
      m_checkBox->setVisible(realVisibility());
      m_checkBox->raise();
    });
  }

  ~ComponentLegendItem() override { delete m_checkBox.data(); }

  void sync()
  {
    const QSignalBlocker block(m_checkBox);
    m_checkBox->setChecked(mPlottable->visible());
    m_checkBox->setText(mPlottable->name());
    m_checkBox->setAccessibleName(QObject::tr("Show %1").arg(mPlottable->name()));
    m_checkBox->setToolTip(QObject::tr("Show or hide %1 points, error bars and response curves on this plot.")
                             .arg(mPlottable->name()));
    auto palette = m_checkBox->palette();
    const auto *graph = static_cast<QCPGraph *>(mPlottable);
    palette.setColor(QPalette::WindowText, graph->visible() ? graph->scatterStyle().pen().color() : QColor("#808080"));
    m_checkBox->setPalette(palette);
  }

protected:
  QSize minimumOuterSizeHint() const override
  {
    return m_checkBox->sizeHint() + QSize(mMargins.left() + mMargins.right(), mMargins.top() + mMargins.bottom());
  }

  void draw(QCPPainter *painter) override
  {
    if(!painter->modes().testFlag(QCPPainter::pmNoCaching)) return;
    QStyleOptionButton option;
    option.initFrom(m_checkBox);
    option.rect = rect();
    option.text = m_checkBox->text();
    option.state &= ~(QStyle::State_On | QStyle::State_Off | QStyle::State_HasFocus | QStyle::State_MouseOver);
    option.state |= m_checkBox->isChecked() ? QStyle::State_On : QStyle::State_Off;
    painter->save();
    painter->setFont(m_checkBox->font());
    m_checkBox->style()->drawControl(QStyle::CE_CheckBox, &option, painter, m_checkBox);
    painter->restore();
  }

private:
  QPointer<QCheckBox> m_checkBox;
};
}

namespace PlotColors
{
namespace
{
const QColor &componentXX()
{
  static const QColor color("#579BC9");
  return color;
}

const QColor &componentXY()
{
  static const QColor color("#1D5FA2");
  return color;
}

const QColor &componentYX()
{
  static const QColor color("#B45309");
  return color;
}

const QColor &componentYY()
{
  static const QColor color("#D08738");
  return color;
}
}

const std::vector<QColor> &componentColors()
{
  // Tensor rows form two hue families: XX/XY blue, YX/YY orange.
  // Keep the off-diagonal components darker than their diagonal partners.
  static const std::vector<QColor> colors = {
    componentXX(),
    componentXY(),
    componentYX(),
    componentYY()
  };
  return colors;
}

const std::vector<QColor> &tipperScalarColors()
{
  // Preserve the established tipper palette independently of impedance colors.
  static const std::vector<QColor> colors = {
    QColor("#2563EB"),
    QColor("#D97706"),
    QColor("#60A5FA"),
    QColor("#FBBF24")
  };
  return colors;
}

QColor masked()
{
  static const QColor color("#9CA3AF");
  return color;
}

QColor tipperReal()
{
  static const QColor color("#111827");
  return color;
}

QColor tipperImaginary()
{
  static const QColor color("#6B7280");
  return color;
}

QColor tipperMasked()
{
  static const QColor color("#BFC4CC");
  return color;
}
}

MTDataPlot::MTDataPlot(QCustomPlot *plot):
  m_plot(plot), m_associated_plot(nullptr),
  m_yAxisAutoscale(true), m_fixedYRange(plot->yAxis->range())
{
  m_contextMenu = new QMenu(m_plot);
  m_contextMenu->addAction(tr("Mask"), this, &MTDataPlot::maskSelectedData);
  m_contextMenu->addAction(tr("Inverted Mask"), this, &MTDataPlot::invMaskSelectedData);
  m_contextMenu->addAction(tr("Unmask"), this, &MTDataPlot::unmaskSelectedData);

  m_plot->setFocusPolicy(Qt::StrongFocus);

  connect(m_plot, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(showPointToolTip(QMouseEvent*)));
}

MTDataPlot::~MTDataPlot() = default;

void MTDataPlot::set_associated_plot(MTDataPlot &plot)
{
  m_associated_plot = &plot;
}

void MTDataPlot::set_error_bars_visible(bool on)
{
  m_errorBarsVisible = on;
  apply_component_visibility();
}

void MTDataPlot::set_component_visibility(const std::array<bool, 4> &visible)
{
  m_componentVisible = visible;
  apply_component_visibility();
  m_plot->replot();
}

void MTDataPlot::set_legend_component_visible(unsigned component, bool visible)
{
  auto components = m_componentVisible;
  components.at(component) = visible;
  set_component_visibility(components);
  emit componentVisibilityChanged();
}

void MTDataPlot::rebuild_component_legend(const std::vector<int> &graphIndices)
{
  m_plot->legend->clearItems();
  for(int component: graphIndices)
    m_plot->legend->addItem(new ComponentLegendItem(m_plot->legend, m_plot->graph(component), component,
      [this, component](bool visible) { set_legend_component_visible(component, visible); }));
}

void MTDataPlot::update_component_legend()
{
  for(int i = 0; i < m_plot->legend->itemCount(); ++i)
    if(auto *item = dynamic_cast<ComponentLegendItem *>(m_plot->legend->item(i))) item->sync();
}

void MTDataPlot::apply_component_visibility()
{
  for(int i = 0; i < m_plot->graphCount(); ++i) {
    auto *graph = m_plot->graph(i);
    const bool visible = m_componentVisible[i % 4];
    graph->setVisible(visible);
    if(!visible) graph->setSelection(QCPDataSelection());
  }
  for(unsigned i = 0; i < m_errorBars.size(); ++i)
    m_errorBars[i]->setVisible(m_errorBarsVisible && m_componentVisible[i % 4]);
  update_component_legend();
}

void MTDataPlot::set_masking_mode(bool on)
{
  if(on)
    m_plot->setSelectionRectMode(QCP::srmSelect);
  else
    m_plot->setSelectionRectMode(QCP::srmNone);
}

void MTDataPlot::set_y_axis_autoscale(bool on)
{
  m_yAxisAutoscale = on;
}

bool MTDataPlot::y_axis_autoscale() const
{
  return m_yAxisAutoscale;
}

void MTDataPlot::set_y_axis_range(double lower, double upper)
{
  m_fixedYRange = QCPRange(lower, upper);
}

QCPRange MTDataPlot::y_axis_range() const
{
  return m_yAxisAutoscale ? m_plot->yAxis->range() : m_fixedYRange;
}

QCPRange MTDataPlot::fixed_y_axis_range() const
{
  return m_fixedYRange;
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

std::vector<RealDataType> MTDataPlot::get_graph_data_types(const QCPGraph *graph) const
{
  auto itype = m_name2type.find(graph->name().toStdString());
  if(itype == m_name2type.end())
    return {};

  return {itype->second};
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
    bool has_values = false;
    for (unsigned j = 0; j < frequencies.size(); ++j)
    {
      x.push_back(1.0 / frequencies[j]);
      const bool valid = std::isfinite(data[i][j]);
      y.push_back(valid ? data[i][j] : std::numeric_limits<double>::quiet_NaN());
      has_values |= valid;
    }

    if(!has_values) { x.clear(); y.clear(); }
    m_plot->graph(i + data.size()*2)->setData(x, y);
  }
}

void MTDataPlot::clear_predicted_data()
{
  for(int i = static_cast<int>(m_errorBars.size()); i < m_plot->graphCount(); ++i)
    m_plot->graph(i)->data()->clear();
}

void MTDataPlot::set_layout_generic(const std::vector<QString> &data_graph_names,
                                    const std::vector<QString> &masked_graph_names)
{
  const auto &colors = PlotColors::componentColors();
  const QColor maskedColor = PlotColors::masked();

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
    graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, maskedColor, 5));
    graph->setSelectable(QCP::stMultipleDataRanges);
    graph->selectionDecorator()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, maskedColor, 8));
    graph->removeFromLegend();

    QCPScatterStyle style = graph->selectionDecorator()->scatterStyle();
    style.setSize(9);
    graph->selectionDecorator()->setScatterStyle(style, QCPScatterStyle::spSize);
    graph->selectionDecorator()->setUsedScatterProperties(QCPScatterStyle::spSize);

    connect(graph, SIGNAL(selectionChanged(bool)), this, SLOT(dataSelected(bool)));

    m_errorBars.push_back(new QCPErrorBars(m_plot->xAxis, m_plot->yAxis));
    m_errorBars.back()->removeFromLegend();
    m_errorBars.back()->setPen(QPen(maskedColor));
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
  rebuild_component_legend({0, 1, 2, 3});

  m_plot->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(m_plot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(plotContextRequest(QPoint)));
}

void MTDataPlot::apply_axis_ranges(bool rescaleAxes, bool useDefaultYRange,
                                   const QCPRange &defaultYRange)
{
  if(rescaleAxes && std::any_of(m_componentVisible.begin(), m_componentVisible.end(), [](bool visible) { return visible; }))
  {
    m_plot->xAxis->rescale(true);
    QCPRange xrange = m_plot->xAxis->range();
    m_plot->xAxis->setRange(xrange.lower / 2., xrange.upper * 2.);

    if(m_yAxisAutoscale)
    {
      if(useDefaultYRange)
        m_plot->yAxis->setRange(defaultYRange);
      else
        m_plot->yAxis->rescale(true);
    }
  }

  if(!m_yAxisAutoscale)
    m_plot->yAxis->setRange(m_fixedYRange);
}

void MTDataPlot::mask_selected_data(bool on)
{
  QList<QCPGraph*> selectedGraphs = m_plot->selectedGraphs();
  for(auto &graph: selectedGraphs)
  {
    QCPDataSelection selection = graph->selection();
    const std::vector<RealDataType> types = get_graph_data_types(graph);
    if(types.empty())
      continue;

    for (QCPDataRange dataRange: selection.dataRanges())
    {
      auto begin = graph->data()->at(dataRange.begin());
      auto end = graph->data()->at(dataRange.end());
      for (auto it = begin; it != end; ++it)
      {
        for(const auto type: types)
          m_data->set_data_mask(type, 1. / it->key, on, m_linkTensorMasks);
      }
    }
  }

  set_observed_data(*m_data, false);
  if(m_associated_plot)
    m_associated_plot->set_observed_data(*m_associated_plot->m_data, false);
  emit observationsChanged();
}

void MTDataPlot::maskSelectedData()
{
  mask_selected_data(false);
}

void MTDataPlot::invMaskSelectedData()
{
  const auto &frequencies = m_data->frequencies();
  std::map<RealDataType, std::vector<bool>> selected;
  for(auto *graph: m_plot->selectedGraphs()) {
    const auto types = get_graph_data_types(graph);
    for(auto type: types) if(!selected.count(type)) selected[type].assign(frequencies.size(), false);
    for(const auto &range: graph->selection().dataRanges()) {
      for(auto point = graph->data()->at(range.begin()); point != graph->data()->at(range.end()); ++point) {
        const double frequency = 1. / point->key;
        int closestIndex = -1; double closest = 1e-3;
        for(unsigned i = 0; i < frequencies.size(); ++i) {
          const double distance = std::abs(frequency - frequencies[i]) / frequencies[i];
          if(distance < closest) { closest = distance; closestIndex = i; }
        }
        if(closestIndex >= 0) for(auto type: types) selected[type][closestIndex] = true;
      }
    }
  }
  // A masked input wins over an unmask when several selected components share
  // a dependent PT/Z entry. Compute the whole operation before modifying data.
  std::map<RealDataType, std::vector<bool>> related;
  for(const auto &entry: selected) {
    for(unsigned i = 0; i < frequencies.size(); ++i) m_data->set_data_mask(entry.first, frequencies[i], entry.second[i]);
    if(m_linkTensorMasks) for(auto type: MTStationData::linked_mask_types(entry.first)) {
      if(!related.count(type)) related[type].assign(frequencies.size(), true);
      for(unsigned i = 0; i < frequencies.size(); ++i) related[type][i] = related[type][i] && entry.second[i];
    }
  }
  for(const auto &entry: related)
    for(unsigned i = 0; i < frequencies.size(); ++i) m_data->set_data_mask(entry.first, frequencies[i], entry.second[i]);

  set_observed_data(*m_data, false);

  if(m_associated_plot)
    m_associated_plot->set_observed_data(*m_associated_plot->m_data, false);
  emit observationsChanged();
}

void MTDataPlot::unmaskSelectedData()
{
  mask_selected_data(true);
}
