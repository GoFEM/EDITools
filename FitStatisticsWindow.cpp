#include "include/MapBackground.h"
#include "include/PlotColorMaps.h"
#include "include/FileLabels.h"
#include "FitStatisticsWindow.h"
#include <QJsonArray>
#include "include/ToolTips.h"
#include <QDialogButtonBox>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QMessageBox>
#include <QPrinter>
#include <QPageLayout>
#include <QPageSize>
#include <QPushButton>
#include <QRegularExpression>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTabWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <set>

namespace {
QCPColorScale *colorScale(QCustomPlot *plot)
{
  auto *scale = new QCPColorScale(plot);
  plot->plotLayout()->addElement(1, 1, scale);
  scale->setType(QCPAxis::atRight);
  scale->setGradient(PlotColorMaps::gradient());
  scale->axis()->setLabel("nRMS");
  return scale;
}
void clear(QCustomPlot *plot)
{
  plot->clearPlottables();
  plot->clearItems();
}
void target(QCustomPlot *plot)
{
  auto *line = new QCPItemStraightLine(plot);
  line->point1->setCoords(0., 1.);
  line->point2->setCoords(1., 1.);
  line->setPen(QPen(QColor("#808080"), 1, Qt::DotLine));
}
void finish(QCustomPlot *plot, bool zero = true)
{
  plot->rescaleAxes();
  const double upper = std::max(1.1, plot->yAxis->range().upper * 1.1);
  if(zero) plot->yAxis->setRange(0., upper);
  plot->replot();
}
QString shortName(const QString &path)
{
  const auto file = FileLabels::fileName(path);
  const auto match = QRegularExpression("_predicted_iter(\\d+)\\.txt$").match(file);
  return match.hasMatch() ? QObject::tr("Iteration %1").arg(match.captured(1).toInt()) : file;
}
}

FitStatisticsWindow::FitStatisticsWindow(QWidget *parent, std::function<void()> refresh, bool central)
  : WorkspaceView(parent)
{
  setObjectName("fitStatisticsWindow");
  setWindowTitle(tr("Data fit statistics"));
  resize(1450, 940);
  auto *layout = new QVBoxLayout(this);
  auto *controls = new QHBoxLayout;
  errorSource = new QComboBox(this);
  errorSource->setObjectName("fitErrorSource");
  errorSource->addItems({tr("Response"), tr("Observed")});
  controls->addWidget(ToolTips::label(errorSource, tr("Errors:"), tr("nRMS = √mean(((predicted − observed) / error)²). Choose inversion/response errors or current observation errors, including error floors. Only matching, unmasked data are compared. N is matched / available response scalars; the dotted line marks nRMS = 1.")));
  controls->addWidget(errorSource);
  auto *refreshButton = new QPushButton(tr("Refresh"), this);
  controls->addWidget(refreshButton);

  controls->addWidget(new QLabel(tr("Colors:"), this));
  curveColors = new QComboBox(this);
  curveColors->setObjectName("fitCurveColors");
  curveColors->addItem(tr("Distinct"));
  curveColors->addItems(PlotColorMaps::names());
  controls->addWidget(curveColors);
  controls->addStretch();
  auto *pdf = new QPushButton(tr("Save PDF…"), this);
  pdf->setObjectName("fitSavePdf");
  controls->addWidget(pdf);
  layout->addLayout(controls);
  auto *splitter = new QSplitter(this);
  layout->addWidget(splitter, 1);
  auto *left = new QWidget(splitter);
  auto *leftLayout = new QVBoxLayout(left);
  leftLayout->setContentsMargins(0, 0, 0, 0);
  responseList = new QTreeWidget(left);
  responseList->setObjectName("fitResponses");
  responseList->setHeaderLabels({tr("Response"), tr("nRMS"), tr("N")});
  responseList->headerItem()->setToolTip(1, tr("Root mean square of residuals divided by the selected errors."));
  responseList->headerItem()->setToolTip(2, tr("Matched / available response scalars. Counts may differ between responses because their coverage differs."));
  responseList->setRootIsDecorated(false);
  responseList->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
  leftLayout->addWidget(responseList);
  auto *selection = new QHBoxLayout;
  for(bool checked: {true, false}) {
    auto *button = new QPushButton(checked ? tr("Check all") : tr("Uncheck all"), left);
    selection->addWidget(button);
    connect(button, &QPushButton::clicked, this, [this, checked] {
      {
        const QSignalBlocker block(responseList);
        for(int i = 0; i < responseList->topLevelItemCount(); ++i)
          responseList->topLevelItem(i)->setCheckState(0, checked ? Qt::Checked : Qt::Unchecked);
      }
      updatePlots();
    });
  }
  leftLayout->addLayout(selection);
  auto *tabs = new QTabWidget(splitter);
  tabs->setObjectName("fitTabs");
  auto *overview = new QWidget(tabs);
  auto *grid = new QGridLayout(overview);
  overall = makePlot("fitOverall", tr("Overall misfit"), tr("Response number"), "nRMS");
  periods = makePlot("fitPeriods", tr("Misfit by period"), tr("Period (s)"), "nRMS");
  histogram = makePlot("fitHistogram", tr("Normalized residuals — all values retained"), tr("(Predicted − observed) / error"), tr("Scalar count"));
  components = makePlot("fitComponents", tr("Misfit by component"), tr("Component"), "nRMS");
  periods->xAxis->setScaleType(QCPAxis::stLogarithmic);
  periods->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
  grid->addWidget(overall, 0, 0); grid->addWidget(periods, 0, 1);
  grid->addWidget(histogram, 1, 0); grid->addWidget(components, 1, 1);
  tabs->addTab(overview, tr("Overview"));
  auto *stationTab = new QWidget(tabs);
  auto *stationLayout = new QVBoxLayout(stationTab);
  stations = makePlot("fitStations", tr("Misfit by station"), tr("Station, north to south"), "nRMS");
  stationLayout->addWidget(stations);
  tabs->addTab(stationTab, tr("Stations"));
  auto *spatial = new QWidget(tabs);
  auto *spatialGrid = new QGridLayout(spatial);
  detailA = new QComboBox(spatial); detailA->setObjectName("fitDetailA");
  detailB = new QComboBox(spatial); detailB->setObjectName("fitDetailB");
  spatialGrid->addWidget(detailA, 0, 0); spatialGrid->addWidget(detailB, 0, 1);
  auto *colorControls = new QWidget(spatial);
  auto *colorLayout = new QGridLayout(colorControls);
  colorLayout->setContentsMargins(0, 0, 0, 0);
  colorMap = new QComboBox(colorControls);
  colorMap->setObjectName("fitColorMap");
  colorMap->addItems(PlotColorMaps::names());
  colorLayout->addWidget(ToolTips::label(colorMap, tr("Colormap:"), tr("Choose two checked responses above to compare. Both panels share color limits. Auto fits the checked responses; manual limits clip to the endpoint colors. Grey means no matching observations.")), 0, 0);
  colorLayout->addWidget(colorMap, 0, 1);
  reverseColors = new QCheckBox(tr("Reverse"), colorControls);
  reverseColors->setObjectName("fitReverseColors");
  colorLayout->addWidget(reverseColors, 0, 2);
  stationNames = new QCheckBox(tr("Station names"), colorControls);
  stationNames->setObjectName("fitStationNames");
  stationNames->setChecked(true);
  colorLayout->addWidget(stationNames, 0, 3);
  auto addLimits = [&](int row, const QString &title, const QString &prefix, ColorLimits &limits) {
    colorLayout->addWidget(new QLabel(title, colorControls), row, 0);
    limits.automatic = new QCheckBox(tr("Auto"), colorControls);
    limits.automatic->setObjectName(prefix + "Auto");
    limits.automatic->setChecked(true);
    ToolTips::set(limits.automatic, tr("Fit both panels to the checked responses. Turn off to set shared color limits."));
    colorLayout->addWidget(limits.automatic, row, 1);
    limits.minimum = new QDoubleSpinBox(colorControls);
    limits.maximum = new QDoubleSpinBox(colorControls);
    limits.minimum->setObjectName(prefix + "Min");
    limits.maximum->setObjectName(prefix + "Max");
    for(auto *spin: {limits.minimum, limits.maximum}) {
      spin->setRange(0., 1e12);
      spin->setDecimals(3);
      ToolTips::set(spin, tr("Shared color limit when Auto is off. Values outside the range use the endpoint colors."));
      spin->setSingleStep(.25);
      spin->setKeyboardTracking(false);
      spin->setEnabled(false);
    }
    limits.minimum->setMaximum(1e12 - .001);
    limits.maximum->setMinimum(.001);
    limits.maximum->setValue(1.);
    colorLayout->addWidget(new QLabel(tr("Min:"), colorControls), row, 2);
    colorLayout->addWidget(limits.minimum, row, 3);
    colorLayout->addWidget(new QLabel(tr("Max:"), colorControls), row, 4);
    colorLayout->addWidget(limits.maximum, row, 5);
    const auto fields = limits;
    connect(limits.automatic, &QCheckBox::toggled, this, [this, fields](bool automatic) {
      fields.minimum->setEnabled(!automatic); fields.maximum->setEnabled(!automatic);
      updateDetails();
    });
    connect(limits.minimum, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this, fields](double value) {
      if(value >= fields.maximum->value()) {
        const QSignalBlocker block(fields.maximum);
        fields.maximum->setValue(value + .001);
      }
      updateDetails();
    });
    connect(limits.maximum, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this, fields](double value) {
      if(value <= fields.minimum->value()) {
        const QSignalBlocker block(fields.minimum);
        fields.minimum->setValue(std::max(0., value - .001));
      }
      updateDetails();
    });
  };
  addLimits(1, tr("Map nRMS:"), "fitMapRange", mapLimits);
  addLimits(2, tr("Heatmap nRMS:"), "fitHeatRange", heatLimits);
  colorLayout->setColumnStretch(6, 1);
  spatialGrid->addWidget(colorControls, 1, 0, 1, 2);
  mapA = makePlot("fitMapA", tr("Station nRMS"), "", "");
  mapB = makePlot("fitMapB", tr("Station nRMS"), "", "");
  backgroundA = new MapBackground(mapA); backgroundB = new MapBackground(mapB);
  colorLayout->addWidget(MapBackground::button(this, backgroundA->layers(), [this](int mask) {
    backgroundA->setLayers(mask); backgroundB->setLayers(mask); mapA->replot(); mapB->replot();
  }), 0, 6);
  heatA = makePlot("fitHeatmapA", tr("Station / period nRMS"), tr("Period (s)"), tr("Station"));
  heatB = makePlot("fitHeatmapB", tr("Station / period nRMS"), tr("Period (s)"), tr("Station"));
  mapScaleA = colorScale(mapA); mapScaleB = colorScale(mapB);
  heatScaleA = colorScale(heatA); heatScaleB = colorScale(heatB);
  for(auto *plot: {mapA, mapB})
    connect(plot, &QCustomPlot::afterLayout, this, [plot] {
      // Keep distances to scale on screen and in exported pages, expanding
      // whichever range needs room so a narrow survey is never clipped.
      const double width = plot->axisRect()->width(), height = plot->axisRect()->height();
      if(width <= 0 || height <= 0) return;
      if(plot->xAxis->range().size() / width > plot->yAxis->range().size() / height)
        plot->yAxis->setScaleRatio(plot->xAxis);
      else
        plot->xAxis->setScaleRatio(plot->yAxis);
    });
  auto *spatialTabs = new QTabWidget(spatial);
  spatialTabs->setObjectName("fitSpatialTabs");
  auto *maps = new QWidget(spatialTabs);
  auto *mapLayout = new QHBoxLayout(maps);
  mapLayout->addWidget(mapA); mapLayout->addWidget(mapB);
  if(!central) spatialTabs->addTab(maps, tr("Station maps"));
  else { setCentralMaps({mapA, mapB}, {tr("Response A"), tr("Response B")}, tr("Fit appearance")); maps->hide(); }
  auto *heatmaps = new QWidget(spatialTabs);
  auto *heatLayout = new QHBoxLayout(heatmaps);
  heatLayout->addWidget(heatA); heatLayout->addWidget(heatB);
  spatialTabs->addTab(heatmaps, tr("Station / period heatmaps"));
  spatialGrid->addWidget(spatialTabs, 2, 0, 1, 2);
  spatialGrid->setRowStretch(2, 1);
  tabs->addTab(spatial, tr("Maps and heatmaps"));
  splitter->setSizes({340, 1100});
  summary = new QLabel(this);
  summary->setObjectName("fitSummary");
  summary->setWordWrap(true);
  layout->addWidget(summary);
  connect(responseList, &QTreeWidget::itemChanged, this, [this] { updatePlots(); });
  connect(detailA, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updateDetails(); });
  connect(detailB, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updateDetails(); });
  connect(errorSource, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [refresh] { refresh(); });
  connect(refreshButton, &QPushButton::clicked, this, [refresh] { refresh(); });
  connect(pdf, &QPushButton::clicked, this, [this] { savePdf(); });

  connect(curveColors, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updatePlots(); });
  connect(colorMap, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updateDetails(); });
  connect(reverseColors, &QCheckBox::toggled, this, [this] { updateDetails(); });
  connect(stationNames, &QCheckBox::toggled, this, [this] { updateDetails(); });
  addToolGroup(tr("Fit tools"), tr("Data and output"), controls);
  addToolGroup(tr("Fit appearance"), tr("Spatial colors"), colorLayout);
  spatialGrid->removeWidget(colorControls); delete colorControls;
  auto *details = new QHBoxLayout;
  spatialGrid->removeWidget(detailA); spatialGrid->removeWidget(detailB);
  details->addWidget(new QLabel(tr("Left:"), this)); details->addWidget(detailA);
  details->addWidget(new QLabel(tr("Right:"), this)); details->addWidget(detailB);
  addToolGroup(tr("Fit tools"), tr("Spatial comparison"), details);
  splitter->setObjectName("fitMainSplit");
  setupRanges(); setExportButton("pdf", pdf);
  setStateRefresh([this] { updatePlots(); refreshRanges(); });

}

QCustomPlot *FitStatisticsWindow::makePlot(const QString &name, const QString &title,
                                         const QString &x, const QString &y)
{
  auto *plot = new QCustomPlot(this);
  plot->setObjectName(name);
  plot->setMinimumSize(260, 210);
  plot->plotLayout()->insertRow(0);
  plot->plotLayout()->addElement(0, 0, new QCPTextElement(plot, title, QFont(font().family(), 10, QFont::Bold)));
  plot->xAxis->setLabel(x); plot->yAxis->setLabel(y);
  plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  plot->axisRect()->setMinimumMargins(QMargins(55, 10, 15, 45));
  return plot;
}

void FitStatisticsWindow::setData(const std::shared_ptr<MTSurveyData> &data,
                                 const std::map<std::string, MTSurveyData> &loaded)
{
  survey = data;
  backgroundA->setTopography(survey ? survey->topography() : Topography::State{});
  backgroundB->setTopography(survey ? survey->topography() : Topography::State{});
  std::map<QString, Qt::CheckState> states;
  for(int i = 0; i < responseList->topLevelItemCount(); ++i) {
    const auto *item = responseList->topLevelItem(i);
    states[item->data(0, Qt::UserRole).toString()] = item->checkState(0);
  }
  const QSignalBlocker block(responseList);
  responseList->clear();
  responses.clear(); positions.clear(); stationOrder.clear();
  if(survey) {
    stationOrder = survey->get_stations_names();
    auto locations = survey->geographic_locations();
    auto coordinates = survey->coordinates();
    try {
      if(!coordinates.utm) coordinates = SurveyCoordinates::suggested(locations);
      locations = coordinates.transform(locations);
      renderedCoordinates = coordinates;
      backgroundA->setCoordinates(coordinates, .001); backgroundB->setCoordinates(coordinates, .001);
      mapXLabel = coordinates.centered ? tr("East offset (km)") : tr("Easting (km)");
      mapYLabel = coordinates.centered ? tr("North offset (km)") : tr("Northing (km)");
      for(auto &position: locations) { position[0] /= 1000.; position[1] /= 1000.; }
    } catch(const std::exception &) {
      renderedCoordinates = {};
      backgroundA->setCoordinates(SurveyCoordinates{}); backgroundB->setCoordinates(SurveyCoordinates{});
      mapXLabel = tr("Longitude (°)"); mapYLabel = tr("Latitude (°)");
    }
    for(unsigned i = 0; i < stationOrder.size(); ++i)
      if(std::isfinite(locations[i][0]) && std::isfinite(locations[i][1]))
        positions[stationOrder[i]] = locations[i];
    std::stable_sort(stationOrder.begin(), stationOrder.end(), [this](const std::string &a, const std::string &b) {
      const auto pa = positions.find(a), pb = positions.find(b);
      if(pa == positions.end() || pb == positions.end()) return pa != positions.end() && pb == positions.end();
      return pa->second[0] > pb->second[0];
    });
  }
  for(const auto &entry: loaded) {
    Response response;
    response.path = QString::fromStdString(entry.first);
    response.name = shortName(response.path);
    response.color = QColor::fromHsv((responses.size() * 137 + 215) % 360, 175, 185);
    if(survey) response.fit = FitStatistics::compare(*survey, entry.second,
      errorSource->currentIndex() == 0 ? FitStatistics::ErrorSource::Response : FitStatistics::ErrorSource::Observed);
    auto *item = new QTreeWidgetItem(responseList);
    item->setText(0, QString::number(responses.size() + 1) + ". " + response.name);
    item->setToolTip(0, response.path);
    item->setData(0, Qt::UserRole, response.path);
    item->setForeground(0, response.color);
    item->setCheckState(0, states.count(response.path) ? states.at(response.path) : Qt::Checked);
    item->setText(1, response.fit.total.count ? QString::number(response.fit.total.rms(), 'f', 4) : tr("N/A"));
    item->setData(1, Qt::UserRole, response.fit.total.rms());
    item->setText(2, QString("%1 / %2").arg(response.fit.total.count).arg(response.fit.available));
    responses.push_back(std::move(response));
  }
  updatePlots();
}

std::vector<int> FitStatisticsWindow::enabled() const
{
  std::vector<int> indices;
  for(int i = 0; i < responseList->topLevelItemCount(); ++i)
    if(responseList->topLevelItem(i)->checkState(0) == Qt::Checked) indices.push_back(i);
  return indices;
}

void FitStatisticsWindow::updatePlots()
{
  const auto active = enabled();
  {
    const QSignalBlocker block(responseList);
    auto colors = PlotColorMaps::gradient(curveColors->currentText());
    for(unsigned i = 0; i < responses.size(); ++i) {
      // Keep pale palette endpoints away from line plots on a white background.
      responses[i].color = curveColors->currentIndex() == 0 ? QColor::fromHsv((i * 137 + 215) % 360, 175, 185) :
        QColor::fromRgb(colors.color(responses.size() == 1 ? .5 : .8 * i / (responses.size() - 1), QCPRange(0., 1.)));
      responseList->topLevelItem(i)->setForeground(0, responses[i].color);
    }
  }
  for(auto *plot: {overall, periods, histogram, components, stations}) clear(plot);
  auto overallTicks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  auto componentTicks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  auto stationTicks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  std::vector<std::string> componentNames;
  std::set<std::string> presentComponents;
  double minimum = -1., maximum = 1.;
  for(int i: active) {
    for(const auto &component: responses[i].fit.components) presentComponents.insert(component.first);
    for(double r: responses[i].fit.residuals) { minimum = std::min(minimum, r); maximum = std::max(maximum, r); }
  }
  for(auto type: {RealZxx, RealZxy, RealZyx, RealZyy, RealTzx, RealTzy,
                 PTxx, PTxy, PTyx, PTyy, RhoZxx, RhoZxy, RhoZyx, RhoZyy,
                 PhsZxx, PhsZxy, PhsZyx, PhsZyy}) {
    const auto name = FitStatistics::component_name(type);
    if(presentComponents.erase(name)) componentNames.push_back(name);
  }
  for(unsigned i = 0; i < componentNames.size(); ++i) componentTicks->addTick(i, QString::fromStdString(componentNames[i]));
  for(unsigned i = 0; i < stationOrder.size(); ++i) stationTicks->addTick(i, QString::fromStdString(stationOrder[i]));
  constexpr int bins = 60;
  const double binWidth = (maximum - minimum) / bins;
  QVector<double> overallX, overallY;
  for(unsigned order = 0; order < active.size(); ++order) {
    const auto i = active[order];
    const auto &response = responses[i];
    const auto &fit = response.fit;
    overallTicks->addTick(i, QString::number(i+1));
    if(!fit.total.count) continue;
    overallX.push_back(i); overallY.push_back(fit.total.rms());
    auto *point = overall->addGraph();
    point->setData(QVector<double>{double(i)}, QVector<double>{fit.total.rms()});
    point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, response.color, response.color, 8));
    point->setLineStyle(QCPGraph::lsNone);
    point->setName(response.name);
    auto *line = periods->addGraph();
    QVector<double> x, y;
    for(const auto &period: fit.periods) { x.push_back(period.first); y.push_back(period.second.rms()); }
    line->setData(x, y); line->setName(response.name); line->setPen(QPen(response.color, 1.5));
    line->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 3));
    QVector<double> counts(bins, 0), centers;
    for(double r: fit.residuals) ++counts[std::min(bins - 1, std::max(0, int((r - minimum) / binWidth)))];
    for(int bin = 0; bin < bins; ++bin) centers.push_back(minimum + (bin + .5) * binWidth);
    auto *distribution = histogram->addGraph();
    distribution->setData(centers, counts); distribution->setName(response.name);
    distribution->setLineStyle(QCPGraph::lsStepCenter); distribution->setPen(QPen(response.color, 1.5));
    QColor fill = response.color; fill.setAlpha(25); distribution->setBrush(fill);
    auto *bars = new QCPBars(components->xAxis, components->yAxis);
    bars->setName(response.name); bars->setPen(Qt::NoPen); bars->setBrush(response.color);
    const double width = .8 / active.size(); bars->setWidth(width);
    x.clear(); y.clear();
    for(unsigned c = 0; c < componentNames.size(); ++c) {
      const auto found = fit.components.find(componentNames[c]);
      if(found != fit.components.end()) { x.push_back(c - .4 + (order + .5) * width); y.push_back(found->second.rms()); }
    }
    bars->setData(x, y);
    x.clear(); y.clear();
    for(unsigned s = 0; s < stationOrder.size(); ++s) {
      x.push_back(s);
      const auto found = fit.stations.find(stationOrder[s]);
      y.push_back(found == fit.stations.end() ? std::numeric_limits<double>::quiet_NaN() : found->second.rms());
    }
    auto *stationLine = stations->addGraph(); stationLine->setName(response.name);
    stationLine->setData(x, y); stationLine->setPen(QPen(response.color, 1.5));
    stationLine->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 4));
  }
  auto *progression = overall->addGraph();
  progression->setData(overallX, overallY); progression->setPen(QPen(QColor("#777777"), 1, Qt::DashLine));
  overall->xAxis->setTicker(overallTicks);
  components->xAxis->setTicker(componentTicks);
  stations->xAxis->setTicker(stationTicks); stations->xAxis->setTickLabelRotation(75);
  stations->xAxis->setTickLabelFont(QFont(font().family(), 8));
  for(auto *plot: {overall, periods, components, stations}) { target(plot); finish(plot); }
  finish(histogram);
  components->xAxis->setRange(-.5, std::max(.5, double(componentNames.size()) - .5));
  components->replot();
  const auto periodRange = periods->xAxis->range();
  periods->xAxis->setRange(periodRange.lower / 1.1, periodRange.upper * 1.1);
  periods->replot();
  overall->xAxis->setRange(-.5, std::max(.5, double(responses.size()) - .5));
  overall->replot();
  applyAxisRanges();
  for(auto *combo: {detailA, detailB}) {
    const auto previous = combo->currentData().toString();
    const QSignalBlocker block(combo);
    combo->clear();
    for(int i: active) combo->addItem(QString::number(i+1) + ". " + responses[i].name, responses[i].path);
    const int index = combo->findData(previous);
    combo->setCurrentIndex(index >= 0 ? index : (combo == detailA ? 0 : combo->count() - 1));
  }
  updateDetails();
  summary->setText(!survey || responses.empty() ? tr("Load observations and responses.") :
    tr("%1/%2 responses checked · %3 stations")
      .arg(active.size()).arg(responses.size()).arg(stationOrder.size()));
}

void FitStatisticsWindow::updateDetails()
{
  for(const auto &limits: {mapLimits, heatLimits}) {
    setToolVisible(limits.minimum, !limits.automatic->isChecked()); setToolVisible(limits.maximum, !limits.automatic->isChecked());
    limits.minimum->setEnabled(!limits.automatic->isChecked()); limits.maximum->setEnabled(!limits.automatic->isChecked());
  }

  double maximum = 1., stationMaximum = 1.;
  for(int i: enabled()) {
    for(const auto &station: responses[i].fit.stations)
      stationMaximum = std::max(stationMaximum, station.second.rms());
    for(const auto &station: responses[i].fit.station_periods)
      for(const auto &period: station.second) maximum = std::max(maximum, period.second.rms());
  }
  auto selected = [this](QComboBox *combo) {
    for(unsigned i = 0; i < responses.size(); ++i)
      if(responses[i].path == combo->currentData().toString()) return int(i);
    return -1;
  };
  auto range = [](const ColorLimits &limits, double maximum) {
    if(limits.automatic->isChecked()) {
      const QSignalBlocker minBlock(limits.minimum), maxBlock(limits.maximum);
      limits.minimum->setValue(0.); limits.maximum->setValue(maximum);
      return QCPRange(0., maximum);
    }
    return QCPRange(limits.minimum->value(), limits.maximum->value());
  };
  const auto a = selected(detailA), b = selected(detailB);
  const auto mapRange = range(mapLimits, stationMaximum), heatRange = range(heatLimits, maximum);
  drawMap(mapA, mapScaleA, a, mapRange); drawMap(mapB, mapScaleB, b, mapRange);
  drawHeatmap(heatA, heatScaleA, a, heatRange); drawHeatmap(heatB, heatScaleB, b, heatRange);
}

void FitStatisticsWindow::setupRanges()
{
  auto *fields = new QHBoxLayout;
  rangePlot = new QComboBox(this); rangePlot->setObjectName("fitRangePlot");
  rangePlot->addItems({tr("Overall"), tr("By period"), tr("Histogram"), tr("By component"), tr("By station")});
  rangeAxis = new QComboBox(this); rangeAxis->setObjectName("fitRangeAxis"); rangeAxis->addItems({tr("X axis"), tr("Y axis")});
  rangeAuto = new QCheckBox(tr("Auto range"), this); rangeAuto->setObjectName("fitRangeAuto");
  rangeMin = new QDoubleSpinBox(this); rangeMin->setObjectName("fitRangeMin");
  rangeMax = new QDoubleSpinBox(this); rangeMax->setObjectName("fitRangeMax");
  for(auto *spin: {rangeMin, rangeMax}) { spin->setRange(-1e18, 1e18); spin->setDecimals(8); spin->setKeyboardTracking(false); }
  fields->addWidget(rangePlot); fields->addWidget(rangeAxis); fields->addWidget(rangeAuto);
  fields->addWidget(new QLabel(tr("Min:"), this)); fields->addWidget(rangeMin);
  fields->addWidget(new QLabel(tr("Max:"), this)); fields->addWidget(rangeMax);
  auto *reset = new QPushButton(tr("Fit all axes"), this); reset->setObjectName("fitAxesAutoscale"); fields->addWidget(reset);
  rangeError = new QLabel(this); rangeError->setObjectName("fitAxisRangeError"); fields->addWidget(rangeError);
  addToolGroup(tr("Fit appearance"), tr("Axis ranges"), fields);
  connect(rangePlot, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { refreshRanges(); });
  connect(rangeAxis, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { refreshRanges(); });
  const auto apply = [this] {
    const std::array<QCustomPlot *, 5> plots{{overall, periods, histogram, components, stations}};
    auto *plot = plots[rangePlot->currentIndex()]; const int index = rangeAxis->currentIndex();
    auto &option = axisOptions[plot->objectName()][index];
    const auto *axis = index == 0 ? plot->xAxis : plot->yAxis;
    const double low = rangeMin->value(), high = rangeMax->value();
    if(!rangeAuto->isChecked() && (!QCPRange::validRange(low, high) || low >= high || (axis->scaleType() == QCPAxis::stLogarithmic && low <= 0))) {
      rangeError->setText(tr("Enter ordered limits; log axes require positive values.")); return;
    }
    rangeError->clear(); option.automatic = rangeAuto->isChecked(); option.range = QCPRange(low, high);
    updatePlots(); refreshRanges();
  };
  connect(rangeAuto, &QCheckBox::toggled, this, [apply](bool) { apply(); });
  for(auto *spin: {rangeMin, rangeMax}) connect(spin, &QDoubleSpinBox::editingFinished, this, apply);
  connect(reset, &QPushButton::clicked, this, [this] { axisOptions.clear(); updatePlots(); refreshRanges(); });
  refreshRanges();
}
void FitStatisticsWindow::refreshRanges()
{
  const std::array<QCustomPlot *, 5> plots{{overall, periods, histogram, components, stations}};
  auto *plot = plots[rangePlot->currentIndex()]; const int index = rangeAxis->currentIndex();
  const auto &option = axisOptions[plot->objectName()][index];
  const auto range = option.automatic ? (index == 0 ? plot->xAxis : plot->yAxis)->range() : option.range;
  const QSignalBlocker a(rangeAuto), b(rangeMin), c(rangeMax);
  rangeAuto->setChecked(option.automatic); rangeMin->setValue(range.lower); rangeMax->setValue(range.upper);
  setToolVisible(rangeMin, !option.automatic); setToolVisible(rangeMax, !option.automatic);
}
QJsonObject FitStatisticsWindow::saveViewState() const
{
  auto state = WorkspaceView::saveViewState(); QJsonObject axes;
  for(const auto &entry: axisOptions) {
    QJsonArray ranges;
    for(const auto &option: entry.second) ranges.append(QJsonObject{{"auto", option.automatic}, {"min", option.range.lower}, {"max", option.range.upper}});
    axes[entry.first] = ranges;
  }
  state["axes"] = axes; return state;
}
void FitStatisticsWindow::restoreViewState(const QJsonObject &state)
{
  const auto axes = state["axes"].toObject(); axisOptions.clear();
  for(auto it = axes.begin(); it != axes.end(); ++it) {
    const auto ranges = it.value().toArray(); if(ranges.size() != 2) continue;
    for(unsigned i = 0; i < 2; ++i) {
      const auto value = ranges[i].toObject(); const double lo = value["min"].toDouble(), hi = value["max"].toDouble();
      if(QCPRange::validRange(lo, hi) && lo < hi && !(it.key() == periods->objectName() && i == 0 && lo <= 0))
        axisOptions[it.key()][i] = {value["auto"].toBool(true), QCPRange(lo, hi)};
    }
  }
  WorkspaceView::restoreViewState(state); refreshRanges();
}

void FitStatisticsWindow::applyAxisRanges()
{
  for(auto *plot: {overall, periods, histogram, components, stations}) {
    const auto &options = axisOptions[plot->objectName()];
    if(!options[0].automatic) plot->xAxis->setRange(options[0].range);
    if(!options[1].automatic) plot->yAxis->setRange(options[1].range);
    plot->replot();
  }
}

QCPColorGradient FitStatisticsWindow::spatialGradient() const
{
  auto colors = PlotColorMaps::gradient(colorMap->currentText());
  return reverseColors->isChecked() ? colors.inverted() : colors;
}

void FitStatisticsWindow::drawMap(QCustomPlot *plot, QCPColorScale *scale, int response, const QCPRange &range)
{
  clear(plot);
  auto colors = spatialGradient();
  scale->setGradient(colors);
  scale->setDataRange(range);
  plot->xAxis->setLabel(mapXLabel); plot->yAxis->setLabel(mapYLabel);
  const FitStatistics::Result *fit = response < 0 ? nullptr : &responses[response].fit;
  for(const auto &position: positions) {
    QColor color("#bbbbbb");
    if(fit) {
      const auto found = fit->stations.find(position.first);
      if(found != fit->stations.end()) color = QColor::fromRgb(colors.color(found->second.rms(), range));
    }
    auto *point = plot->addGraph();
    point->setData(QVector<double>{position.second[1]}, QVector<double>{position.second[0]});
    point->setLineStyle(QCPGraph::lsNone);
    point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor("#555555"), color, 9));
    point->setName(QString::fromStdString(position.first));
    if(!stationNames->isChecked()) continue;
    auto *label = new QCPItemText(plot);
    label->position->setCoords(position.second[1], position.second[0]);
    label->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom);
    label->setText(QString::fromStdString(position.first).section('_', -1));
    label->setFont(QFont(font().family(), 7));
    label->setPadding(QMargins(3, 0, 0, 2));
  }
  plot->rescaleAxes();
  plot->xAxis->scaleRange(1.15); plot->yAxis->scaleRange(1.15);
  plot->replot();
}

void FitStatisticsWindow::drawHeatmap(QCustomPlot *plot, QCPColorScale *scale, int response, const QCPRange &range)
{
  clear(plot);
  const auto colors = spatialGradient();
  scale->setGradient(colors);
  scale->setDataRange(range);
  plot->axisRect()->setBackground(QColor("#dddddd"));
  std::set<double> periodSet;
  for(int i: enabled()) for(const auto &period: responses[i].fit.periods) periodSet.insert(period.first);
  const std::vector<double> periodValues(periodSet.begin(), periodSet.end());
  auto xTicks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  auto yTicks = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  plot->xAxis->setTicker(xTicks); plot->yAxis->setTicker(yTicks);
  const auto nx = periodValues.size(), ny = stationOrder.size();
  if(response >= 0 && nx && ny) {
    auto *map = new QCPColorMap(plot->xAxis, plot->yAxis);
    // With one station or period, use a transparent neighbour so the lone
    // cell remains centred on its tick and occupies the whole visible cell.
    map->data()->setSize(std::max(2, int(nx)), std::max(2, int(ny)));
    map->data()->setRange(QCPRange(0., std::max(1., double(nx)-1)), QCPRange(0., std::max(1., double(ny)-1)));
    map->data()->fill(0.); map->data()->fillAlpha(0);
    for(unsigned x = 0; x < nx; ++x)
      if(x % std::max(1u, unsigned((nx + 5) / 6)) == 0) xTicks->addTick(x, QString::number(periodValues[x], 'g', 3));
    for(unsigned s = 0; s < ny; ++s) {
      const auto y = ny - s - 1;
      yTicks->addTick(y, QString::fromStdString(stationOrder[s]).section('_', -1));
      const auto station = responses[response].fit.station_periods.find(stationOrder[s]);
      if(station == responses[response].fit.station_periods.end()) continue;
      for(unsigned x = 0; x < nx; ++x) {
        const auto period = station->second.find(periodValues[x]);
        if(period == station->second.end()) continue;
        map->data()->setCell(x, y, period->second.rms());
        map->data()->setAlpha(x, y, 255);
      }
    }
    map->setInterpolate(false); map->setTightBoundary(false);
    map->setColorScale(scale); map->setGradient(colors); map->setDataRange(range);
  }
  plot->xAxis->setRange(-.5, std::max(.5, double(nx)-.5));
  plot->yAxis->setRange(-.5, std::max(.5, double(ny)-.5));
  plot->yAxis->setTickLabelFont(QFont(font().family(), 7));
  plot->xAxis->grid()->setVisible(false); plot->yAxis->grid()->setVisible(false);
  plot->replot();
}

void FitStatisticsWindow::savePdf()
{
  QString path = QFileDialog::getSaveFileName(this, tr("Save fit statistics"), exportName("pdf"), tr("PDF (*.pdf)"));
  if(path.isEmpty()) return;
  if(!path.endsWith(".pdf", Qt::CaseInsensitive)) path += ".pdf";
  QPrinter printer(QPrinter::HighResolution);
  printer.setOutputFormat(QPrinter::PdfFormat); printer.setOutputFileName(path);
  printer.setPageSize(QPageSize(QPageSize::A4)); printer.setPageOrientation(QPageLayout::Landscape);
  printer.setResolution(120);
  QCPPainter painter;
  if(!painter.begin(&printer)) {
    QMessageBox::warning(this, tr("Save fit statistics"), tr("Cannot write the PDF file."));
    return;
  }
  const auto pageRect = printer.pageLayout().paintRectPixels(printer.resolution());
  const int width = pageRect.width(), height = pageRect.height();
  const std::vector<std::vector<QCustomPlot *>> pages{{overall, periods, histogram, components}, {stations}, {mapA, mapB}, {heatA, heatB}};
  for(unsigned page = 0; page < pages.size(); ++page) {
    if(page) printer.newPage();
    painter.setPen(Qt::black);
    painter.drawText(QRect(0, 0, width, 24), Qt::AlignLeft, tr("Data fit statistics — %1").arg(errorSource->currentText()));
    painter.setFont(QFont(font().family(), 8));
    int legendX = 0, legendY = 28;
    const int lineHeight = painter.fontMetrics().height() + 5;
    for(int i: enabled()) {
      const auto text = QString("%1: %2 (nRMS %3, N %4)").arg(i+1).arg(responses[i].name)
        .arg(responses[i].fit.total.rms(), 0, 'f', 4).arg(responses[i].fit.total.count);
      const int itemWidth = painter.fontMetrics().horizontalAdvance(text) + 24;
      if(legendX && legendX + itemWidth > width) { legendX = 0; legendY += lineHeight; }
      painter.setPen(responses[i].color);
      painter.drawText(QRect(legendX, legendY, itemWidth, lineHeight), Qt::AlignLeft | Qt::AlignVCenter, text);
      legendX += itemWidth;
    }
    int plotTop = legendY + lineHeight + 12;
    painter.setPen(Qt::black);
    if(page >= 2) {
      painter.drawText(QRect(0, plotTop, width, 24), Qt::AlignCenter, detailA->currentText() + "       |       " + detailB->currentText());
      plotTop += 28;
    }
    const int columns = pages[page].size() == 1 ? 1 : 2;
    const int rows = pages[page].size() > 2 ? 2 : 1;
    for(unsigned i = 0; i < pages[page].size(); ++i) {
      painter.save();
      painter.translate((i % columns) * width / columns, plotTop + (i / columns) * (height - plotTop) / rows);
      pages[page][i]->toPainter(&painter, width / columns, (height - plotTop) / rows);
      painter.restore();
    }
  }
  painter.end();
}
