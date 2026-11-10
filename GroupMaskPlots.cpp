#include "GroupMaskPlots.h"
#include "include/MTSurveyData.h"
#include "include/MTDataPlot.h"
#include <QCheckBox>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QToolTip>
#include <QVBoxLayout>
#include <cmath>

GroupMaskPlots::GroupMaskPlots(QWidget *parent, const MTSurveyData &survey, const std::set<std::string> &stations)
  : QWidget(parent)
{
  auto *layout = new QVBoxLayout(this); layout->setContentsMargins(0, 0, 0, 0);
  auto *key = new QHBoxLayout;
  auto *legend = new QLabel(tr("<span style='color:#777'>× Already masked</span> · <span style='color:#bf233a'>⊗ To mask</span> · <span style='color:#14804a'>○ To unmask</span>"), this);
  legend->setWordWrap(true); key->addWidget(legend, 1);
  auto *reset = new QPushButton(tr("Fit view"), this); reset->setObjectName("groupMaskFitView"); key->addWidget(reset);
  layout->addLayout(key);
  auto *grid = new QGridLayout; layout->addLayout(grid, 1);
  const QStringList axes{tr("Apparent resistivity (Ω m)"), tr("Phase (°)"), tr("Tipper"), tr("Phase tensor")};
  const QStringList tensorLabels{"XX", "XY", "YX", "YY"}, tipperLabels{"Re X", "Re Y", "Im X", "Im Y"};
  for(unsigned pane = 0; pane < plots.size(); ++pane) {
    auto *panel = new QWidget(this); auto *column = new QVBoxLayout(panel); column->setContentsMargins(0, 0, 0, 0);
    auto *visibility = new QHBoxLayout; column->addLayout(visibility);
    auto *plot = plots[pane] = new QCustomPlot(panel); plot->setObjectName("groupMaskPlot" + QString::number(pane));
    plot->setMinimumSize(250, 190); column->addWidget(plot, 1);
    plot->xAxis->setLabel(tr("Period (s)")); plot->yAxis->setLabel(axes[pane]);
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    if(pane == 0) {
      plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
      plot->yAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    }
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot->setToolTip(tr("All selected stations, including disabled sites. Thin curves show usable samples; gaps remain at masked or missing data. Large markers show changes from the current mask plan, including Z/PT links. Hover over a point for its station. Drag or scroll to inspect; this does not edit masks. Phases use the stored convention without wrapping."));
    const auto &colors = pane == 2 ? PlotColors::tipperScalarColors() : PlotColors::componentColors();
    for(unsigned c = 0; c < 4; ++c) {
      auto *box = new QCheckBox(pane == 2 ? tipperLabels[c] : tensorLabels[c], panel);
      box->setObjectName(QString("groupMaskVisible%1_%2").arg(pane).arg(c)); box->setChecked(true);
      box->setToolTip(tr("Show or hide this component on this plot. Visibility does not change the mask plan."));
      auto palette = box->palette(); palette.setColor(QPalette::WindowText, colors[c]); box->setPalette(palette);
      visibility->addWidget(box);
      connect(box, &QCheckBox::toggled, this, [this, pane, c](bool on) {
        for(const auto &curve : curves) if(curve.pane == pane && curve.component == c)
          for(auto *graph : {curve.observed, curve.masked, curve.pending}) graph->setVisible(on);
        plots[pane]->replot();
      });
    }
    visibility->addStretch(); grid->addWidget(panel, pane / 2, pane % 2);
    connect(plot, &QCustomPlot::mouseMove, this, [plot](QMouseEvent *event) {
      int index;
      auto *graph = plot->plottableAt<QCPGraph>(event->pos(), false, &index);
      if(!graph || index < 0 || index >= graph->dataCount()) { QToolTip::hideText(); return; }
      const auto point = graph->data()->at(index);
      QToolTip::showText(event->globalPos(), tr("%1\n%2 s · %3").arg(graph->name()).arg(point->key, 0, 'g', 12).arg(point->value, 0, 'g', 8), plot);
    });
  }
  // Cache original curves once. Mask previews only replace the change markers.
  const std::array<RealDataType, 4> z{{RealZxx, RealZxy, RealZyx, RealZyy}},
    t{{RealTzx, RealTzy, RealTzx, RealTzy}}, pt{{PTxx, PTxy, PTyx, PTyy}};
  using Getter = void (MTStationData::*)(std::vector<dvector> &, std::vector<dvector> &) const;
  const std::array<Getter, 4> getters{{&MTStationData::get_apparent_resistivity, &MTStationData::get_phase,
                                     &MTStationData::get_tipper, &MTStationData::get_phase_tensor}};
  for(const auto &name : stations) if(survey.is_station_present(name)) {
    const auto &station = survey.get_station_data(name);
    for(unsigned pane = 0; pane < plots.size(); ++pane) {
      std::vector<dvector> values, errors; (station.*getters[pane])(values, errors);
      const auto masks = pane < 2 ? station.impedance_mask(true) : pane == 2 ? station.tipper_mask(true) : station.phase_tensor_mask(true);
      const auto &colors = pane == 2 ? PlotColors::tipperScalarColors() : PlotColors::componentColors();
      for(unsigned c = 0; c < 4; ++c) {
        const auto type = pane < 2 ? z[c] : pane == 2 ? t[c] : pt[c];
        auto *plot = plots[pane];
        const QString title = QString::fromStdString(name) + " · " + (pane == 2 ? tipperLabels[c] : tensorLabels[c]);
        auto graph = [&](const char *state, const QColor &color, QCPScatterStyle::ScatterShape shape, double size) {
          auto *g = plot->addGraph(); g->setName(title); g->setProperty("station", QString::fromStdString(name));
          g->setProperty("component", c); g->setProperty("maskState", state);
          g->setPen(QPen(color, 1.)); g->setScatterStyle(QCPScatterStyle(shape, color, size)); g->setLineStyle(QCPGraph::lsNone); return g;
        };
        auto *observed = graph("observed", colors[c], QCPScatterStyle::ssCircle, 3.);
        observed->setLineStyle(QCPGraph::lsLine);
        auto *masked = graph("masked", QColor("#777777"), QCPScatterStyle::ssCross, 5.);
        auto *pending = graph("pending", QColor("#bf233a"), QCPScatterStyle::ssCrossCircle, 9.);
        pending->setLayer("overlay");
        Curve curve{name, type, pane, c, observed, masked, pending, {}, {}, {}};
        for(unsigned i = 0; i < station.frequencies().size(); ++i) {
          const double f = station.frequencies()[i], p = 1. / f;
          if(!std::isfinite(f) || f <= 0. || !std::isfinite(p) || p <= 0.) continue;
          const double value = values[c][i];
          const bool valid = std::isfinite(value) && (pane != 0 || value > 0.);
          observed->addData(p, valid && masks[c][i] ? value : NAN);
          if(valid && !masks[c][i]) masked->addData(p, value);
          curve.periods.push_back(p); curve.values.push_back(valid ? value : NAN); curve.frequencies.push_back(f);
        }
        curves.push_back(std::move(curve));
      }
    }
  }
  connect(reset, &QPushButton::clicked, this, [this] { fitAxes(); });
  for(auto *plot : plots) connect(plot->xAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, [this, plot](const QCPRange &range) {
    for(auto *other : plots) if(other != plot) { const QSignalBlocker block(other->xAxis); other->xAxis->setRange(range); other->replot(); }
  });
  fitAxes();
}

void GroupMaskPlots::setEdits(const GroupMaskEdits &edits, bool unmask)
{
  std::map<std::string, std::set<std::pair<RealDataType, double>>> selected;
  for(const auto &entry : edits) selected[entry.first].insert(entry.second.begin(), entry.second.end());
  const QColor color(unmask ? "#14804a" : "#bf233a");
  for(const auto &curve : curves) {
    curve.pending->data()->clear(); curve.pending->setPen(QPen(color, 1.8));
    curve.pending->setScatterStyle(QCPScatterStyle(unmask ? QCPScatterStyle::ssCircle : QCPScatterStyle::ssCrossCircle, QPen(color, 1.8), Qt::NoBrush, 9.));
    curve.pending->setName(curve.observed->name() + (unmask ? tr(" · to unmask") : tr(" · to mask")));
    const auto match = selected.find(curve.station);
    if(match == selected.end()) continue;
    for(int i = 0; i < curve.periods.size(); ++i)
      if(std::isfinite(curve.values[i]) && match->second.count({curve.type, curve.frequencies[i]}))
        curve.pending->addData(curve.periods[i], curve.values[i]);
  }
  for(auto *plot : plots) plot->replot();
}

void GroupMaskPlots::fitAxes()
{
  double low = INFINITY, high = 0.;
  for(const auto &curve : curves) for(double p : curve.periods) { low = std::min(low, p); high = std::max(high, p); }
  for(unsigned pane = 0; pane < plots.size(); ++pane) {
    auto *plot = plots[pane]; const QSignalBlocker block(plot->xAxis);
    plot->yAxis->setRange(pane == 0 ? .1 : -1., pane == 0 ? 100. : 1.);
    plot->rescaleAxes(true);
    plot->yAxis->scaleRange(1.1);
    plot->xAxis->setRange(std::isfinite(low) ? low / 1.1 : .1, high > 0. ? high * 1.1 : 10.);
    plot->replot();
  }
}
