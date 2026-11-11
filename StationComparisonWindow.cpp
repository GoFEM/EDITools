#include "StationComparisonWindow.h"
#include "include/MTSurveyData.h"
#include "include/MapBackground.h"
#include "qcustomplot.h"
#include "include/SurveyGeometry.h"
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QVBoxLayout>
#include <cmath>
#include <numeric>

namespace {
QColor stationColor(int index) {
  return QColor::fromHsv((index * 137) % 360, 180, 175);
}
using SurveyGeometry::located;
using SurveyGeometry::distance;
} // namespace
class StationComparisonWindow::Overlay : public QCPLayerable {
public:
  struct Station { QString name; QPointF position; QColor color; bool checked; };
  std::vector<Station> stations;
  explicit Overlay(QCustomPlot *map) : QCPLayerable(map) { setObjectName("comparisonMapOverlay"); setLayer("overlay"); setVisible(false); }
  QPointF pixel(const Station &station) const { return {mParentPlot->xAxis->coordToPixel(station.position.x()), mParentPlot->yAxis->coordToPixel(station.position.y())}; }
protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(true); }
  void draw(QCPPainter *p) override {
    for(const auto &station: stations) {
      const auto point = pixel(station);
      p->setPen(QPen(station.color)); p->setBrush(station.color);
      const double radius = station.checked ? 4.5 : 2.5;
      p->drawEllipse(point, radius, radius);
      if(station.checked) p->drawText(point + QPointF(7., -5.), station.name);
    }
  }
};
StationComparisonWindow::StationComparisonWindow(QWidget *parent, QCustomPlot *shared) : WorkspaceView(parent), sharedMap(shared) {
  setObjectName("stationComparison");
  setWindowTitle(tr("Neighbouring stations"));
  resize(1180, 760);
  auto *layout = new QVBoxLayout(this);
  auto *controls = new QHBoxLayout;
  layout->addLayout(controls);
  controls->addWidget(new QLabel(tr("Reference"), this));
  anchor = new QComboBox(this);
  anchor->setObjectName("compareAnchor");
  controls->addWidget(anchor);
  count = new QSpinBox(this);
  count->setObjectName("compareCount");
  count->setRange(1, 100);
  count->setValue(5);
  controls->addWidget(count);
  auto *nearest = new QPushButton(tr("Nearest"), this);
  nearest->setObjectName("compareNearest");
  controls->addWidget(nearest);
  nearest->setToolTip(
      tr("Select the reference and this many nearest enabled neighbours, using WGS84 geodesic distance. "
         "Stations without geographic positions can still be checked manually."));
  family = new QComboBox(this);
  family->setObjectName("compareFamily");
  family->addItems({tr("Resistivity / phase"), tr("Tipper"), tr("Phase tensor")});
  controls->addWidget(family);
  component = new QComboBox(this);
  component->setObjectName("compareComponent");
  controls->addWidget(component);
  errors = new QCheckBox(tr("Errors"), this);
  errors->setObjectName("compareErrors");
  controls->addWidget(errors);
  controls->addStretch();
  auto *splitter = new QSplitter(this);
  layout->addWidget(splitter, 1);
  auto *left = new QWidget(splitter);
  auto *leftLayout = new QVBoxLayout(left);
  leftLayout->setContentsMargins(0, 0, 0, 0);
  stations = new QListWidget(left);
  stations->setObjectName("compareStations");
  leftLayout->addWidget(stations, 1);
  stations->setToolTip(tr("One checkbox per station controls its curves. Colours identify stations. "
                          "Distances are from the reference; phases retain the stored convention. Masked, "
                          "missing and disabled samples are omitted without connecting across gaps."));
  map = new QCustomPlot(left);
  map->setObjectName(sharedMap ? "compareExportMap" : "compareMap");
  map->xAxis->ticker()->setTickCount(3);
  map->setMinimumSize(280, 240);
  if(sharedMap) map->hide(); else leftLayout->addWidget(map, 1);
  background = new MapBackground(map);
  map->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  map->setToolTip(tr("Click a station to toggle its curves. Pan and zoom to inspect the survey."));
  auto *right = new QWidget(splitter);
  auto *plots = new QVBoxLayout(right);
  plots->setContentsMargins(0, 0, 0, 0);
  first = new QCustomPlot(right);
  first->setObjectName("compareFirstPlot");
  second = new QCustomPlot(right);
  second->setObjectName("compareSecondPlot");
  for (auto *plot : {first, second}) {
    plots->addWidget(plot, 1);
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot->xAxis->setLabel(tr("Period (s)"));
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
  }
  splitter->setSizes({330, 820});
  splitter->setStretchFactor(1, 1);
  status = new QLabel(this);
  status->setObjectName("compareStatus");
  layout->addWidget(status);
  connect(nearest, &QPushButton::clicked, this, [this] { selectNearest(); });
  connect(anchor, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { rebuild(); });
  connect(stations, &QListWidget::itemChanged, this, [this] { draw(); });
  connect(errors, &QCheckBox::toggled, this, [this] { draw(); });
  connect(component, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { draw(); });
  connect(family, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] {
    const QSignalBlocker block(component);
    component->clear();
    component->addItems(family->currentIndex() == 1 ? QStringList{"X", "Y"}
                                                    : QStringList{"XX", "XY", "YX", "YY"});
    component->setCurrentIndex(family->currentIndex() == 1 ? 0 : 1);
    draw();
  });
  connect(map, &QCustomPlot::plottableClick, this, [this](QCPAbstractPlottable *plot, int, QMouseEvent *) {
    for (int i = 0; i < stations->count(); ++i)
      if (stations->item(i)->data(Qt::UserRole).toString() == plot->name()) {
        auto *item = stations->item(i);
        if (item->flags() & Qt::ItemIsEnabled)
          item->setCheckState(item->checkState() == Qt::Checked ? Qt::Unchecked : Qt::Checked);
        break;
      }
  });
  if(sharedMap) {
    overlay = new Overlay(sharedMap);
    connect(sharedMap, &QCustomPlot::mousePress, this, [this](QMouseEvent *event) { if(sharedActive) sharedPress = event->position().toPoint(); });
    connect(sharedMap, &QCustomPlot::mouseRelease, this, [this](QMouseEvent *event) {
      if(!sharedActive || !overlay || event->button() != Qt::LeftButton || (event->position().toPoint() - sharedPress).manhattanLength() > 4) return;
      if(!sharedMap->axisRect()->rect().contains(event->position().toPoint())) return;
      QString nearest; double best = 12.;
      for(const auto &station: overlay->stations) {
        const auto delta = overlay->pixel(station) - event->position(); const double distance = std::hypot(delta.x(), delta.y());
        if(distance < best) { best = distance; nearest = station.name; }
      }
      for(int i = 0; !nearest.isEmpty() && i < stations->count(); ++i) {
        auto *item = stations->item(i);
        if(item->data(Qt::UserRole).toString() == nearest && (item->flags() & Qt::ItemIsEnabled)) {
          item->setCheckState(item->checkState() == Qt::Checked ? Qt::Unchecked : Qt::Checked);
          stations->setCurrentItem(item); emit stationActivated(nearest); break;
        }
      }
    });
  }
  component->addItems({"XX", "XY", "YX", "YY"});
  component->setCurrentIndex(1);
  addToolGroup(tr("Comparison tools"), tr("Stations and curves"), controls);
  setStateRefresh([this] {
    const auto selected = component->currentText(); const QSignalBlocker block(component);
    component->clear(); component->addItems(family->currentIndex() == 1 ? QStringList{"X", "Y"} : QStringList{"XX", "XY", "YX", "YY"});
    component->setCurrentIndex(std::max(0, component->findText(selected))); rebuild();
  });
  connect(anchor, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { emit stationActivated(anchor->currentText()); });
}
StationComparisonWindow::~StationComparisonWindow() {
  setSharedMapActive(false, sharedCoordinates); delete overlay.data();
}
void StationComparisonWindow::setSharedMapActive(bool active, const SurveyCoordinates &coordinates) {
  if(!sharedMap || !overlay) return;
  if(!active) {
    if(sharedActive) { sharedMap->setInteractions(savedInteractions); overlay->setVisible(false); sharedMap->replot(QCustomPlot::rpQueuedReplot); }
    sharedActive = false; return;
  }
  if(!sharedActive) savedInteractions = sharedMap->interactions();
  sharedActive = true; sharedCoordinates = coordinates;
  sharedMap->setInteraction(QCP::iSelectPlottables, false);
  updateSharedMap();
}
void StationComparisonWindow::updateSharedMap() {
  if(!sharedActive || !overlay) return;
  overlay->stations.clear();
  if(survey) for(int i = 0; i < stations->count(); ++i) {
    const auto *item = stations->item(i); const auto name = item->data(Qt::UserRole).toString();
    const auto &station = survey->get_station_data(name.toStdString());
    if(!located(station.position())) continue;
    try {
      const auto p = sharedCoordinates.transform({station.position()}).front();
      overlay->stations.push_back({name, QPointF(p[1], p[0]), item->foreground().color(), item->checkState() == Qt::Checked});
    } catch(const std::exception &) {}
  }
  overlay->setVisible(true); sharedMap->replot(QCustomPlot::rpQueuedReplot);
}
void StationComparisonWindow::setData(std::shared_ptr<MTSurveyData> data, const QStringList &selected) {
  const bool changed = survey != data;
  survey = std::move(data);
  const QString previous = selected.isEmpty() ? anchor->currentText() : selected.front();
  {
    const QSignalBlocker block(anchor);
    anchor->clear();
    if (survey)
      for (const auto &name : survey->get_stations_names())
        anchor->addItem(QString::fromStdString(name));
    anchor->setCurrentIndex(std::max(0, anchor->findText(previous)));
  }
  if (changed)
    stations->clear();
  rebuild();
  if (selected.size() > 1) {
    {
      const QSignalBlocker block(stations);
      for (int i = 0; i < stations->count(); ++i)
        stations->item(i)->setCheckState(selected.contains(stations->item(i)->data(Qt::UserRole).toString())
                                             ? Qt::Checked
                                             : Qt::Unchecked);
    }
    draw();
  } else if (changed || !selected.isEmpty())
    selectNearest();
}
void StationComparisonWindow::rebuild() {
  QStringList checked;
  for (int i = 0; i < stations->count(); ++i)
    if (stations->item(i)->checkState() == Qt::Checked)
      checked << stations->item(i)->data(Qt::UserRole).toString();
  {
    const QSignalBlocker block(stations);
    stations->clear();
    if (survey && survey->is_station_present(anchor->currentText().toStdString())) {
      const auto &position = survey->get_station_data(anchor->currentText().toStdString()).position();
      std::vector<std::pair<double, std::string>> ordered;
      for (const auto &name : survey->get_stations_names())
        ordered.emplace_back(distance(position, survey->get_station_data(name).position()), name);
      std::sort(ordered.begin(), ordered.end());
      for (const auto &entry : ordered) {
        const auto name = QString::fromStdString(entry.second);
        auto *item = new QListWidgetItem(name + (std::isfinite(entry.first)
                                                     ? tr(" · %1 km").arg(entry.first, 0, 'f', 2)
                                                     : tr(" · no location")),
                                         stations);
        item->setData(Qt::UserRole, name);
        item->setCheckState(checked.contains(name) ? Qt::Checked : Qt::Unchecked);
        if (!survey->is_active(entry.second)) {
          item->setText(item->text() + tr(" · disabled"));
          item->setCheckState(Qt::Unchecked);
          item->setFlags(item->flags() & ~Qt::ItemIsEnabled);
        }
      }
    }
  }
  draw();
}
void StationComparisonWindow::selectNearest() {
  int remaining = count->value();
  const bool referenceLocated =
      survey && survey->is_station_present(anchor->currentText().toStdString()) &&
      located(survey->get_station_data(anchor->currentText().toStdString()).position());
  {
    const QSignalBlocker block(stations);
    for (int i = 0; i < stations->count(); ++i) {
      auto *item = stations->item(i);
      const auto name = item->data(Qt::UserRole).toString();
      const bool reference = name == anchor->currentText();
      const bool select = (item->flags() & Qt::ItemIsEnabled) &&
                          (reference || (referenceLocated && remaining > 0 &&
                                         located(survey->get_station_data(name.toStdString()).position())));
      item->setCheckState(select ? Qt::Checked : Qt::Unchecked);
      if (select && !reference)
        --remaining;
    }
  }
  draw();
}
void StationComparisonWindow::draw() {
  const QSignalBlocker block(stations);
  for (auto *plot : {map, first, second})
    plot->clearPlottables();
  if (!survey || !stations->count()) {
    for (auto *plot : {map, first, second})
      plot->replot();
    status->clear();
    updateSharedMap();
    return;
  }
  const int mode = family->currentIndex(), c = std::max(0, component->currentIndex());
  const RealDataType rho[]{RhoZxx, RhoZxy, RhoZyx, RhoZyy}, phase[]{PhsZxx, PhsZxy, PhsZyx, PhsZyy},
      pt[]{PTxx, PTxy, PTyx, PTyy};
  const RealDataType tipReal[]{RealTzx, RealTzy}, tipImag[]{ImagTzx, ImagTzy};
  const std::array<RealDataType, 2> types{{mode == 0   ? rho[c]
                                           : mode == 1 ? tipReal[c]
                                                       : pt[c],
                                           mode == 0   ? phase[c]
                                           : mode == 1 ? tipImag[c]
                                                       : pt[c]}};
  first->yAxis->setScaleType(mode == 0 ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
  first->yAxis->setTicker(mode == 0 ? QSharedPointer<QCPAxisTicker>(new QCPAxisTickerLog)
                                    : QSharedPointer<QCPAxisTicker>(new QCPAxisTicker));
  first->yAxis->setLabel(mode == 0   ? tr("Apparent resistivity (Ω m)")
                         : mode == 1 ? tr("Real tipper")
                                     : tr("Phase tensor"));
  second->yAxis->setLabel(mode == 0 ? tr("Phase (°)") : tr("Imaginary tipper"));
  second->setVisible(mode != 2);
  auto coordinates = survey->coordinates();
  background->setCoordinates(coordinates);
  background->setTopography(survey->topography());
  map->xAxis->setLabel(coordinates.utm ? tr("East (m)") : tr("Longitude (°)"));
  map->yAxis->setLabel(coordinates.utm ? tr("North (m)") : tr("Latitude (°)"));
  const auto &origin = survey->get_station_data(anchor->currentText().toStdString()).position();
  unsigned visible = 0, samples = 0;
  double eastMin=INFINITY,eastMax=-INFINITY,northMin=INFINITY,northMax=-INFINITY;
  const auto names = survey->get_stations_names();
  for (int i = 0; i < stations->count(); ++i) {
    auto *item = stations->item(i);
    const auto name = item->data(Qt::UserRole).toString().toStdString();
    const auto &station = survey->get_station_data(name);
    const bool on = item->checkState() == Qt::Checked && station.active();
    const auto color = stationColor(int(std::lower_bound(names.begin(), names.end(), name) - names.begin()));
    item->setForeground(on ? color : QColor(Qt::gray));
    if (located(station.position()))
      try {
        auto p = coordinates.transform({station.position()}).front();
        if (!coordinates.utm && located(origin))
          p[1] = origin[1] + std::remainder(p[1] - origin[1], 360.);
        eastMin=std::min(eastMin,p[1]);eastMax=std::max(eastMax,p[1]);
        northMin=std::min(northMin,p[0]);northMax=std::max(northMax,p[0]);
        auto *g = map->addGraph();
        g->setName(QString::fromStdString(name));
        g->setData(QVector<double>{p[1]}, QVector<double>{p[0]});
        g->setLineStyle(QCPGraph::lsNone);
        g->setScatterStyle(
            QCPScatterStyle(QCPScatterStyle::ssDisc, on ? color : QColor(Qt::gray), on ? 8 : 5));
      } catch (const std::exception &) {
      }
    if (!on)
      continue;
    ++visible;
    for (int axis = 0; axis < (mode == 2 ? 1 : 2); ++axis) {
      auto *plot = axis ? second : first;
      auto *graph = plot->addGraph();
      graph->setName(QString::fromStdString(name));
      graph->setPen(QPen(color, 1.2));
      graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
      QVector<double> x, y, e;
      for (unsigned f = 0; f < station.frequencies().size(); ++f) {
        const double frequency = station.frequencies()[f];
        if (!std::isfinite(frequency) || frequency <= 0.)
          continue;
        double value = NAN, error = NAN;
        const bool valid =
            station.scalar_value(types[axis], f, value, error) && !(mode == 0 && axis == 0 && value <= 0.);
        x << 1. / frequency;
        y << (valid ? value : NAN);
        e << (valid && std::isfinite(error) && error >= 0. ? error : NAN);
        samples += valid;
      }
      // QCP sorts graph points by period; error vectors must follow the same order.
      std::vector<int> order(x.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int a, int b) { return x[a] < x[b]; });
      QVector<double> sx, sy, se;
      for (int j : order) {
        sx << x[j];
        sy << y[j];
        se << e[j];
      }
      graph->setData(sx, sy, true);
      if (errors->isChecked()) {
        auto *bars = new QCPErrorBars(plot->xAxis, plot->yAxis);
        bars->setDataPlottable(graph);
        bars->setData(se);
        bars->setPen(QPen(color));
      }
    }
  }
  for (auto *plot : {map, first, second}) {
    if (plot->graphCount()) {
      if(plot==map) {
        const double span=std::max({eastMax-eastMin,northMax-northMin,coordinates.utm?100.:.01});
        plot->xAxis->setRange((eastMin+eastMax)/2.,eastMax>eastMin?eastMax-eastMin:span,Qt::AlignCenter);
        plot->yAxis->setRange((northMin+northMax)/2.,northMax>northMin?northMax-northMin:span,Qt::AlignCenter);
      } else {
        // Stable fallback ranges for a single finite value; repeated redraws must not grow them.
        plot->xAxis->setRange(.1,10.);
        plot->yAxis->setRange(plot->yAxis->scaleType()==QCPAxis::stLogarithmic?.1:-1.,
                             plot->yAxis->scaleType()==QCPAxis::stLogarithmic?10.:1.);
        plot->rescaleAxes();
      }
      plot->xAxis->scaleRange(1.1);
      plot->yAxis->scaleRange(1.1);
    }
    plot->replot();
  }
  status->setText(tr("%1 checked stations · %2 visible scalar samples").arg(visible).arg(samples));
  updateSharedMap();
}
