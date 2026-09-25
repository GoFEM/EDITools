#include "include/MapBackground.h"
#include "include/PlotColorMaps.h"
#include "include/FileLabels.h"
#include "PeriodMapWindow.h"
#include "include/HelpButton.h"
#include <QFileDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QPushButton>
#include <QSignalBlocker>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>

namespace {
QColor arrowColor(bool imaginary) { return imaginary ? QColor("#c2410c") : QColor("#111827"); }
}

PeriodMapWindow::PeriodMapWindow(QWidget *parent, std::function<void()> changed)
  : QDialog(parent, Qt::Window), observationsChanged(std::move(changed))
{
  setObjectName("periodMapWindow");
  setWindowTitle(tr("Induction / phase tensor maps"));
  resize(1200, 1000);
  auto *layout = new QVBoxLayout(this);
  auto *top = new QHBoxLayout;
  top->addWidget(new QLabel(tr("Data:"), this));
  dataset = new QComboBox(this); dataset->setObjectName("mapDataset");
  dataset->setMinimumWidth(230);
  top->addWidget(dataset, 1);
  top->addWidget(new QLabel(tr("Period (s):"), this));
  period = new QComboBox(this); period->setObjectName("mapPeriod");
  period->setMinimumWidth(125);
  top->addWidget(period);
  for(int direction: {-1, 1}) {
    auto *button = new QPushButton(direction < 0 ? tr("Previous") : tr("Next"), this);
    button->setObjectName(direction < 0 ? "mapPreviousPeriod" : "mapNextPeriod");
    top->addWidget(button);
    connect(button, &QPushButton::clicked, this, [this, direction] {
      if(period->count()) period->setCurrentIndex(std::max(0, std::min(period->count() - 1, period->currentIndex() + direction)));
    });
  }
  auto *fit = new QPushButton(tr("Fit view"), this); fit->setObjectName("mapFitView"); top->addWidget(fit);
  auto *pdf = new QPushButton(tr("Save PDF…"), this); pdf->setObjectName("mapSavePdf"); top->addWidget(pdf);
  layout->addLayout(top);

  auto *layers = new QHBoxLayout;
  auto addLayer = [&](QCheckBox *&box, const QString &text, const QString &name, bool checked) {
    box = new QCheckBox(text, this); box->setObjectName(name); box->setChecked(checked); layers->addWidget(box);
  };
  addLayer(tensors, tr("Phase tensors"), "mapPhaseTensors", true);
  addLayer(realArrows, tr("Real arrows"), "mapRealArrows", true);
  addLayer(imagArrows, tr("Imaginary arrows"), "mapImagArrows", false);
  addLayer(names, tr("Station names"), "mapStationNames", false);
  layers->addStretch();
  layers->addWidget(UiHelp::label(this, tr("Tolerance (%):"), tr("Use each station's closest period within this percentage of the selected period. No interpolation or extrapolation. Missing, incomplete and masked data are omitted; valid zero vectors are counted but have no arrowhead."), "mapPeriodHelp"));
  tolerance = new QDoubleSpinBox(this); tolerance->setObjectName("mapPeriodTolerance");
  tolerance->setRange(0., 25.); tolerance->setDecimals(3); tolerance->setValue(.1); tolerance->setKeyboardTracking(false);
  tolerance->setToolTip(tr("Use the closest available period within this percentage. No interpolation or extrapolation."));
  layers->addWidget(tolerance);
  layout->addLayout(layers);

  auto *controls = new QGridLayout;
  controls->addWidget(UiHelp::label(this, tr("Convention:"), tr("Parkinson reverses the stored tipper sign; Wiese keeps it. The sign applies to both real and imaginary arrows. Tzx is north and Tzy is east, accounting for the projection's direction of true north."), "mapConventionHelp"), 0, 0);
  convention = new QComboBox(this); convention->setObjectName("mapArrowConvention");
  convention->addItems({tr("Parkinson (−T)"), tr("Wiese (+T)")});
  convention->setToolTip(tr("The sign applies to both real and imaginary arrows. Stored Tzx is north and Tzy is east."));
  controls->addWidget(convention, 0, 1);
  auto sizeControl = [&](const QString &name, double value) {
    auto *spin = new QDoubleSpinBox(this); spin->setObjectName(name);
    spin->setRange(.001, 100000.); spin->setDecimals(3); spin->setValue(value);
    spin->setSingleStep(.1); spin->setKeyboardTracking(false); return spin;
  };
  ellipseSize = sizeControl("mapEllipseSize", .5);
  arrowSize = sizeControl("mapArrowSize", 1.);
  controls->addWidget(UiHelp::label(this, tr("Ellipse (km):"), tr("Major-axis diameter. Ellipses have equal major diameters; their minor/major ratio is |Φmin/Φmax|. Orientation is clockwise from geographic north. Singular or undefined tensors are omitted."), "mapEllipseHelp"), 0, 2);
  controls->addWidget(ellipseSize, 0, 3);
  controls->addWidget(UiHelp::label(this, tr("Arrow scale (km):"), tr("Arrow length for a tipper magnitude |T| = 1. The map reference arrow shows |T| = 0.5."), "mapArrowHelp"), 0, 4);
  controls->addWidget(arrowSize, 0, 5);
  controls->addWidget(UiHelp::label(this, tr("Color by:"), tr("Color represents maximum phase φmax, minimum phase φmin, or skew β, in degrees. Min and Max set the color range; values outside it use endpoint colors."), "mapColorHelp"), 1, 0);
  colorBy = new QComboBox(this); colorBy->setObjectName("mapColorBy");
  colorBy->addItems({tr("φmax (°)"), tr("φmin (°)"), tr("Skew β (°)")});
  controls->addWidget(colorBy, 1, 1);
  controls->addWidget(new QLabel(tr("Colormap:"), this), 1, 2);
  colorMap = new QComboBox(this); colorMap->setObjectName("mapColorMap");
  colorMap->addItems(PlotColorMaps::names());
  controls->addWidget(colorMap, 1, 3);
  reverseColors = new QCheckBox(tr("Reverse"), this); reverseColors->setObjectName("mapReverseColors");
  controls->addWidget(reverseColors, 1, 4);
  auto *limits = new QHBoxLayout;
  colorMin = new QDoubleSpinBox(this); colorMin->setObjectName("mapColorMin");
  colorMax = new QDoubleSpinBox(this); colorMax->setObjectName("mapColorMax");
  for(auto *spin: {colorMin, colorMax}) {
    spin->setRange(-180., 180.); spin->setDecimals(2); spin->setKeyboardTracking(false);
  }
  colorMin->setMaximum(179.99); colorMax->setMinimum(-179.99); colorMax->setValue(90.);
  limits->addWidget(new QLabel(tr("Min:"), this)); limits->addWidget(colorMin);
  limits->addWidget(new QLabel(tr("Max:"), this)); limits->addWidget(colorMax);
  controls->addLayout(limits, 1, 5);
  layout->addLayout(controls);

  auto *maskControls = new QHBoxLayout;
  selectMode = new QCheckBox(tr("Select stations"), this); selectMode->setObjectName("mapSelectStations");
  selectMode->setToolTip(tr("Click a station, arrow or ellipse; drag a box to select several stations. Ctrl-click toggles a station; Ctrl-drag adds stations. Turn off to pan."));
  maskControls->addWidget(selectMode);
  maskControls->addWidget(UiHelp::button(this, tr("Station selection"), selectMode->toolTip(), "mapSelectionHelp"));
  maskControls->addWidget(UiHelp::label(this, tr("Mask group:"), tr("Mask / Unmask edits observed data at the displayed period, within the tolerance, even while viewing a computed response. Computed responses are unchanged. Disabled stations and unmatched observed periods are skipped. Masked stations keep selectable center dots so they can be unmasked.\nTippers affects both directions and real/imaginary parts. Phase tensor affects all four tensor components and, when Data → Link Z / PT masks is on, all impedance components too. Impedance + tensor always affects both tensors and derived resistivity and phase."), "mapMaskHelp"));
  maskGroup = new QComboBox(this); maskGroup->setObjectName("mapMaskGroup");
  maskGroup->addItems({tr("Tippers"), tr("Phase tensor"), tr("Impedance + tensor")});
  maskGroup->setToolTip(tr("Tippers: both directions, real and imaginary. Phase tensor: all four tensor components, plus impedance when Link Z / PT masks is on. Impedance + phase tensor: all impedance and tensor components, including derived resistivity and phase."));
  maskControls->addWidget(maskGroup);
  maskButton = new QPushButton(tr("Mask"), this); maskButton->setObjectName("mapMask");
  unmaskButton = new QPushButton(tr("Unmask"), this); unmaskButton->setObjectName("mapUnmask");
  clearSelection = new QPushButton(tr("Clear selection"), this); clearSelection->setObjectName("mapClearSelection");
  for(auto *button: {maskButton, unmaskButton, clearSelection}) maskControls->addWidget(button);
  maskControls->addStretch();
  layout->addLayout(maskControls);
  selectionSummary = new QLabel(this); selectionSummary->setObjectName("mapSelectionSummary");
  selectionSummary->setWordWrap(false); maskControls->addWidget(selectionSummary);

  plot = new QCustomPlot(this); plot->setObjectName("periodMapPlot");
  background = new MapBackground(plot);
  layers->insertWidget(4, MapBackground::button(this, 0, [this](int mask) { background->setLayers(mask); plot->replot(); }));
  plot->setMinimumSize(650, 450);
  plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  plot->plotLayout()->insertRow(0);
  title = new QCPTextElement(plot, "", QFont(font().family(), 12, QFont::Bold));
  plot->plotLayout()->addElement(0, 0, title);
  colorScale = new QCPColorScale(plot); plot->plotLayout()->addElement(1, 1, colorScale);
  colorScale->setType(QCPAxis::atRight);
  note = new QCPTextElement(plot, "", QFont(font().family(), 9));
  plot->plotLayout()->addElement(2, 0, note);
  plot->xAxis->setLabel(tr("Easting (km)")); plot->yAxis->setLabel(tr("Northing (km)"));
  plot->legend->setVisible(true); plot->legend->setBrush(Qt::white);
  plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight | Qt::AlignTop);
  layout->addWidget(plot, 1);
  summary = new QLabel(this); summary->setObjectName("mapSummary"); summary->setWordWrap(true);
  layout->addWidget(summary);

  connect(selectMode, &QCheckBox::toggled, this, [this](bool selecting) {
    plot->setSelectionRectMode(selecting ? QCP::srmCustom : QCP::srmNone);
    plot->setInteraction(QCP::iRangeDrag, !selecting);
    plot->setCursor(selecting ? Qt::CrossCursor : Qt::ArrowCursor);
  });
  connect(plot, &QCustomPlot::mousePress, this, [this](QMouseEvent *event) {
    selectionStart = event->pos(); selectionDragged = false;
    selectionStarted = selectMode->isChecked() && event->button() == Qt::LeftButton &&
                       plot->axisRect()->rect().contains(selectionStart);
  });
  connect(plot, &QCustomPlot::mouseRelease, this, [this](QMouseEvent *event) {
    if(selectionStarted && !selectionDragged && event->button() == Qt::LeftButton)
      selectStations(QRect(event->pos(), event->pos()), event->modifiers().testFlag(Qt::ControlModifier), false);
  });
  connect(plot->selectionRect(), &QCPSelectionRect::accepted, this, [this](const QRect &rect, QMouseEvent *event) {
    if(selectionStarted && event->button() == Qt::LeftButton)
      selectStations(rect.normalized(), event->modifiers().testFlag(Qt::ControlModifier), true);
  });
  connect(maskButton, &QPushButton::clicked, this, [this] { setSelectedMasks(false); });
  connect(unmaskButton, &QPushButton::clicked, this, [this] { setSelectedMasks(true); });
  connect(clearSelection, &QPushButton::clicked, this, [this] { selectedStations.clear(); updateSelection(); });

  connect(plot, &QCustomPlot::afterLayout, this, [this] {
    const double width = plot->axisRect()->width(), height = plot->axisRect()->height();
    if(width <= 0. || height <= 0.) return;
    if(plot->xAxis->range().size() / width > plot->yAxis->range().size() / height)
      plot->yAxis->setScaleRatio(plot->xAxis);
    else plot->xAxis->setScaleRatio(plot->yAxis);
    if(referenceArrow) {
      const auto start = referenceArrow->start->pixelPosition();
      referenceArrow->end->setCoords(plot->xAxis->pixelToCoord(start.x()) + .5 * arrowSize->value(),
                                     plot->yAxis->pixelToCoord(start.y()));
    }
  });
  connect(plot, &QCustomPlot::mouseMove, this, [this](QMouseEvent *event) {
    if((event->pos() - selectionStart).manhattanLength() > 3) selectionDragged = true;
    QString tooltip;
    double distance = 20.;
    for(const auto &site: sites) {
      const double dx = plot->xAxis->coordToPixel(site.position.x()) - event->pos().x();
      const double dy = plot->yAxis->coordToPixel(site.position.y()) - event->pos().y();
      if(std::hypot(dx, dy) < distance) { distance = std::hypot(dx, dy); tooltip = site.details; }
    }
    plot->setToolTip(tooltip);
  });
  connect(dataset, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updatePeriods(); updateMap(); });
  connect(period, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updateMap(); });
  for(auto *box: {names, tensors, realArrows, imagArrows, reverseColors})
    connect(box, &QCheckBox::toggled, this, [this] { updateMap(); });
  for(auto *spin: {tolerance, ellipseSize, arrowSize})
    connect(spin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this] { updateMap(); });
  for(auto *combo: {convention, colorMap})
    connect(combo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updateMap(); });
  connect(colorBy, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { resetColorRange(); updateMap(); });
  connect(colorMin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double value) {
    if(value >= colorMax->value()) { const QSignalBlocker block(colorMax); colorMax->setValue(value + .01); }
    updateMap();
  });
  connect(colorMax, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double value) {
    if(value <= colorMin->value()) { const QSignalBlocker block(colorMin); colorMin->setValue(value - .01); }
    updateMap();
  });
  connect(fit, &QPushButton::clicked, this, [this] { fitView(); });
  connect(pdf, &QPushButton::clicked, this, [this] { savePdf(); });
}

const MTSurveyData *PeriodMapWindow::source() const
{
  const auto key = dataset->currentData().toString();
  const auto found = responses.find(key);
  return key.isEmpty() ? survey.get() : (found == responses.end() ? nullptr : found->second);
}

void PeriodMapWindow::setData(const std::shared_ptr<MTSurveyData> &data,
                              const std::map<std::string, MTSurveyData> &loaded)
{
  const auto selected = dataset->currentData().toString();
  const auto previousSites = sites;
  if(survey != data) selectedStations.clear();
  survey = data; responses.clear(); sites.clear(); locationError.clear(); coordinatesDescription.clear(); originDescription.clear();
  {
    const QSignalBlocker block(dataset);
    dataset->clear(); dataset->addItem(tr("Observed data"), QString());
    for(const auto &entry: loaded) {
      const auto key = QString::fromStdString(entry.first);
      responses[key] = &entry.second;
      dataset->addItem(FileLabels::fileName(key), key);
      dataset->setItemData(dataset->count() - 1, key, Qt::ToolTipRole);
    }
    dataset->setCurrentIndex(std::max(0, dataset->findData(selected)));
  }
  if(survey) {
    const auto stationNames = survey->get_stations_names();
    const auto geographic = survey->geographic_locations();
    std::vector<std::array<double, 3>> points;
    for(unsigned i = 0; i < geographic.size(); ++i) {
      if(!std::isfinite(geographic[i][0]) || !std::isfinite(geographic[i][1]) || !std::isfinite(geographic[i][2]) ||
         geographic[i][0] < -80. || geographic[i][0] > 84. || std::abs(geographic[i][1]) > 180.) continue;
      sites.push_back({stationNames[i], {}, {}, {}, {}});
      points.push_back(geographic[i]);
    }
    try {
      if(!points.empty()) {
        auto coordinates = survey->coordinates();
        if(!coordinates.utm) coordinates = SurveyCoordinates::suggested(points);
        const auto count = points.size();
        // A short step toward true north gives each station's local projected
        // direction, accounting for UTM meridian convergence.
        for(unsigned i = 0; i < count; ++i) {
          auto north = points[i]; north[0] += north[0] > 83.99998 ? -1e-5 : 1e-5; points.push_back(north);
        }
        const auto projected = coordinates.transform(points);
        for(unsigned i = 0; i < count; ++i) {
          sites[i].position = QPointF(projected[i][1] / 1000., projected[i][0] / 1000.);
          QPointF north(projected[i + count][1] - projected[i][1], projected[i + count][0] - projected[i][0]);
          if(points[i + count][0] < points[i][0]) north = -north;
          const double length = std::hypot(north.x(), north.y());
          if(!std::isfinite(length) || length == 0.) throw std::runtime_error("Cannot determine the local north direction.");
          sites[i].north = north / length;
          sites[i].east = QPointF(sites[i].north.y(), -sites[i].north.x());
        }
        background->setCoordinates(coordinates, .001);
        coordinatesDescription = tr("WGS84 / UTM %1%2%3").arg(coordinates.zone)
          .arg(coordinates.north ? "N" : "S").arg(coordinates.centered ? tr("; survey-centered offsets") : QString());
        originDescription = tr("UTM origin: E = %1 m, N = %2 m")
          .arg(coordinates.origin_easting, 0, 'f', 2).arg(coordinates.origin_northing, 0, 'f', 2);
        plot->xAxis->setLabel(coordinates.centered ? tr("East offset (km)") : tr("Easting (km)"));
        plot->yAxis->setLabel(coordinates.centered ? tr("North offset (km)") : tr("Northing (km)"));
      }
    } catch(const std::exception &error) {
      sites.clear(); locationError = QString::fromUtf8(error.what());
    }
  }
  geometryChanged = sites.size() != previousSites.size();
  for(auto it = selectedStations.begin(); it != selectedStations.end();) {
    if(std::none_of(sites.begin(), sites.end(), [&](const Site &site) { return site.name == *it; }))
      it = selectedStations.erase(it);
    else ++it;
  }
  for(unsigned i = 0; !geometryChanged && i < sites.size(); ++i)
    geometryChanged = sites[i].name != previousSites[i].name || sites[i].position != previousSites[i].position;
  if(geometryChanged && !sites.empty()) {
    std::vector<double> nearest;
    for(unsigned i = 0; i < sites.size(); ++i) {
      double distance = std::numeric_limits<double>::infinity();
      for(unsigned j = 0; j < sites.size(); ++j) if(i != j) {
        const auto delta = sites[i].position - sites[j].position;
        const double length = std::hypot(delta.x(), delta.y());
        if(length > 0.) distance = std::min(distance, length);
      }
      if(std::isfinite(distance)) nearest.push_back(distance);
    }
    std::sort(nearest.begin(), nearest.end());
    const double spacing = nearest.empty() ? 1. : nearest[nearest.size() / 2];
    const QSignalBlocker ellipseBlock(ellipseSize), arrowBlock(arrowSize);
    ellipseSize->setValue(spacing * .6); arrowSize->setValue(spacing * 1.5);
  }
  updatePeriods();
  updateMap(geometryChanged);
}

void PeriodMapWindow::updatePeriods()
{
  const double previous = period->currentData().isValid() ? period->currentData().toDouble() : 1.;
  const QSignalBlocker block(period);
  period->clear();
  if(!source()) return;
  double best = std::numeric_limits<double>::infinity();
  int selected = 0;
  for(double value: source()->get_unique_periods()) {
    if(!std::isfinite(value) || value <= 0.) continue;
    period->addItem(QString::number(value, 'g', 8), value);
    const double difference = std::abs(std::log(value / previous));
    if(difference < best) { best = difference; selected = period->count() - 1; }
  }
  period->setCurrentIndex(selected);
}

void PeriodMapWindow::resetColorRange()
{
  const QSignalBlocker minBlock(colorMin), maxBlock(colorMax);
  colorMin->setValue(colorBy->currentIndex() == 2 ? -10. : (colorBy->currentIndex() == 1 ? -90. : 0.));
  colorMax->setValue(colorBy->currentIndex() == 2 ? 10. : 90.);
}

void PeriodMapWindow::extendBounds(const QPointF &point)
{
  if(!haveBounds) { boundsX = QCPRange(point.x(), point.x()); boundsY = QCPRange(point.y(), point.y()); haveBounds = true; }
  boundsX.lower = std::min(boundsX.lower, point.x()); boundsX.upper = std::max(boundsX.upper, point.x());
  boundsY.lower = std::min(boundsY.lower, point.y()); boundsY.upper = std::max(boundsY.upper, point.y());
}

void PeriodMapWindow::drawEllipse(const Site &site, const MTMapData::PhaseTensor &tensor, double actualPeriod,
                                   QCPColorGradient &colors, const QCPRange &range)
{
  const double angle = tensor.azimuth * std::acos(-1.) / 180.;
  const auto major = site.north * std::cos(angle) + site.east * std::sin(angle);
  const auto minor = site.east * std::cos(angle) - site.north * std::sin(angle);
  QVector<double> order, x, y;
  for(int i = 0; i <= 96; ++i) {
    const double theta = i * 2. * std::acos(-1.) / 96.;
    const auto point = site.position + .5 * ellipseSize->value() *
      (major * std::cos(theta) + minor * (tensor.axisRatio * std::sin(theta)));
    order.push_back(i); x.push_back(point.x()); y.push_back(point.y()); extendBounds(point);
  }
  auto *ellipse = new QCPCurve(plot->xAxis, plot->yAxis);
  ellipse->setObjectName("phaseTensor_" + QString::fromStdString(site.name));
  ellipse->setProperty("stationName", QString::fromStdString(site.name));
  ellipse->setProperty("phiMin", tensor.phiMin); ellipse->setProperty("phiMax", tensor.phiMax);
  ellipse->setProperty("skew", tensor.skew); ellipse->setProperty("azimuth", tensor.azimuth);
  ellipse->setProperty("axisRatio", tensor.axisRatio); ellipse->setProperty("period", actualPeriod);
  ellipse->setData(order, x, y); ellipse->setPen(QPen(QColor("#333333"), .8));
  const double value = colorBy->currentIndex() == 0 ? tensor.phiMax : (colorBy->currentIndex() == 1 ? tensor.phiMin : tensor.skew);
  ellipse->setBrush(QColor::fromRgb(colors.color(value, range)));
  ellipse->setSelectable(QCP::stNone); ellipse->removeFromLegend();
}

void PeriodMapWindow::drawArrow(const Site &site, const std::array<double, 2> &vector, bool imaginary, double actualPeriod)
{
  if(std::hypot(vector[0], vector[1]) == 0.) return; // A valid zero vector has no arrowhead.
  const double sign = convention->currentIndex() == 0 ? -1. : 1.;
  const auto end = site.position + sign * arrowSize->value() * (site.north * vector[0] + site.east * vector[1]);
  if(!std::isfinite(end.x()) || !std::isfinite(end.y())) return;
  auto *arrow = new QCPItemLine(plot);
  arrow->setLayer("axes");
  arrow->setObjectName((imaginary ? "inductionImag_" : "inductionReal_") + QString::fromStdString(site.name));
  arrow->setProperty("stationName", QString::fromStdString(site.name));
  arrow->setProperty("north", vector[0]); arrow->setProperty("east", vector[1]); arrow->setProperty("period", actualPeriod);
  arrow->setSelectable(false);
  arrow->start->setCoords(site.position); arrow->end->setCoords(end);
  arrow->setPen(QPen(arrowColor(imaginary), 1.6));
  arrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8));
  extendBounds(end);
}

void PeriodMapWindow::updateMap(bool fit)
{
  referenceArrow = nullptr;
  selectionGraph = nullptr;
  plot->clearPlottables(); plot->clearItems(); haveBounds = false;
  auto colors = PlotColorMaps::gradient(colorMap->currentText());
  if(reverseColors->isChecked()) colors = colors.inverted();
  const QCPRange range(colorMin->value(), colorMax->value());
  colorScale->setGradient(colors); colorScale->setDataRange(range);
  colorScale->axis()->setLabel(colorBy->currentText());
  colorScale->setVisible(tensors->isChecked());
  const auto *data = source();
  const double requested = period->currentData().toDouble();
  unsigned nTensor = 0, nReal = 0, nImag = 0, matchingPeriods = 0;
  QVector<double> stationX, stationY;
  for(auto &site: sites) {
    extendBounds(site.position);
    stationX.push_back(site.position.x()); stationY.push_back(site.position.y());
    site.details = QString::fromStdString(site.name);
    const auto *station = data && data->is_station_present(site.name) ? &data->get_station_data(site.name) : nullptr;
    const int index = station ? MTMapData::nearest_period(station->frequencies(), requested, tolerance->value() / 100.) : -1;
    if(index >= 0 && survey->is_active(site.name) && station->active()) {
      ++matchingPeriods;
      const double actual = 1. / station->frequencies()[index];
      site.details += tr("\nPeriod: %1 s").arg(actual, 0, 'g', 9);
      MTMapData::PhaseTensor tensor;
      if(MTMapData::phase_tensor(*station, index, tensor)) {
        site.details += tr("\nφmin: %1°; φmax: %2°\nβ: %3°; azimuth: %4°")
          .arg(tensor.phiMin, 0, 'f', 2).arg(tensor.phiMax, 0, 'f', 2).arg(tensor.skew, 0, 'f', 2).arg(tensor.azimuth, 0, 'f', 2);
        if(tensors->isChecked()) { drawEllipse(site, tensor, actual, colors, range); ++nTensor; }
      }
      for(bool imaginary: {false, true}) {
        std::array<double, 2> vector;
        if(MTMapData::induction_vector(*station, index, imaginary, vector)) {
          site.details += tr("\n%1 Tzx: %2; Tzy: %3").arg(imaginary ? tr("Im") : tr("Re"))
            .arg(vector[0], 0, 'g', 5).arg(vector[1], 0, 'g', 5);
          if((imaginary ? imagArrows : realArrows)->isChecked()) {
            drawArrow(site, vector, imaginary, actual);
            if(imaginary) ++nImag; else ++nReal;
          }
        }
      }
    } else site.details += tr("\nNo active data at this period within the tolerance.");
    if(names->isChecked()) {
      auto *label = new QCPItemText(plot); label->setObjectName("periodMapLabel_" + QString::fromStdString(site.name));
      label->position->setCoords(site.position); label->setText(QString::fromStdString(site.name));
      label->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom); label->setPadding(QMargins(4, 0, 0, 3));
      label->setFont(QFont(font().family(), 8)); label->setSelectable(false);
    }
  }
  auto *stations = plot->addGraph(); stations->setName(tr("Station")); stations->setData(stationX, stationY);
  stations->setObjectName("mapStations");
  stations->setLineStyle(QCPGraph::lsNone); stations->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, QColor("#888888"), 3));
  stations->setSelectable(QCP::stNone); stations->removeFromLegend();
  selectionGraph = plot->addGraph(); selectionGraph->setObjectName("mapSelectedStations");
  selectionGraph->setLayer("axes"); selectionGraph->setLineStyle(QCPGraph::lsNone);
  selectionGraph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor("#2563eb"), 2), Qt::NoBrush, 13));
  selectionGraph->setSelectable(QCP::stNone); selectionGraph->removeFromLegend();
  plot->legend->setVisible(realArrows->isChecked() || imagArrows->isChecked());
  for(bool imaginary: {false, true}) if((imaginary ? imagArrows : realArrows)->isChecked()) {
    auto *key = plot->addGraph(); key->setName(imaginary ? tr("Imaginary induction") : tr("Real induction"));
    key->setPen(QPen(arrowColor(imaginary), 1.6)); key->setSelectable(QCP::stNone);
  }
  if(realArrows->isChecked() || imagArrows->isChecked()) {
    referenceArrow = new QCPItemLine(plot); referenceArrow->setObjectName("mapReferenceArrow");
    referenceArrow->start->setType(QCPItemPosition::ptAxisRectRatio); referenceArrow->start->setCoords(.07, .92);
    referenceArrow->setHead(QCPLineEnding(QCPLineEnding::esSpikeArrow, 6, 8)); referenceArrow->setSelectable(false);
    auto *label = new QCPItemText(plot); label->position->setParentAnchor(referenceArrow->start);
    label->position->setType(QCPItemPosition::ptAbsolute); label->position->setCoords(0., -14.);
    label->setText("|T| = 0.5"); label->setPositionAlignment(Qt::AlignLeft | Qt::AlignVCenter); label->setSelectable(false);
  }
  title->setText(tr("%1 — T = %2 s").arg(dataset->currentText()).arg(requested, 0, 'g', 8));
  note->setText(tr("%1   |   %2")
                .arg(coordinatesDescription + "\n" + originDescription).arg(convention->currentText()));
  summary->setText(!locationError.isEmpty() ? tr("Cannot project station locations: %1").arg(locationError) :
    (sites.empty() ? tr("Load stations with geographic locations.") :
     tr("%1/%2 stations · %3 tensors · %4 real / %5 imaginary vectors")
       .arg(matchingPeriods).arg(sites.size()).arg(nTensor).arg(nReal).arg(nImag)));
  plot->setProperty("phaseTensorCount", nTensor); plot->setProperty("realVectorCount", nReal); plot->setProperty("imagVectorCount", nImag);
  updateSelection(false);
  if(fit) fitView(); else plot->replot();
}

int PeriodMapWindow::observedIndex(const std::string &name) const
{
  if(!survey || !survey->is_station_present(name) || !survey->is_active(name)) return -1;
  return MTMapData::nearest_period(survey->get_station_data(name).frequencies(),
                                 period->currentData().toDouble(), tolerance->value() / 100.);
}

void PeriodMapWindow::selectStations(const QRect &rectangle, bool multiple, bool drag)
{
  std::set<std::string> found;
  double closest = 12.;
  for(const auto &site: sites) {
    const QPointF pixel(plot->xAxis->coordToPixel(site.position.x()), plot->yAxis->coordToPixel(site.position.y()));
    if(drag) {
      if(rectangle.contains(pixel.toPoint())) found.insert(site.name);
    } else {
      const auto offset = pixel - rectangle.topLeft();
      const double distance = std::hypot(offset.x(), offset.y());
      if(distance < closest) { closest = distance; found = {site.name}; }
    }
  }
  // Stations remain selectable after masking; glyphs can also be clicked away from their centers.
  if(!drag && found.empty()) {
    auto hit = [&](QCPLayerable *glyph) {
      const auto name = glyph->property("stationName");
      if(!name.isValid()) return;
      const double distance = glyph->selectTest(rectangle.topLeft(), false);
      if(distance >= 0. && distance < closest) { closest = distance; found = {name.toString().toStdString()}; }
    };
    for(int i = 0; i < plot->plottableCount(); ++i) hit(plot->plottable(i));
    for(int i = 0; i < plot->itemCount(); ++i) hit(plot->item(i));
  }
  if(!multiple) selectedStations = found;
  else for(const auto &name: found) {
    if(!drag && selectedStations.count(name)) selectedStations.erase(name);
    else selectedStations.insert(name);
  }
  updateSelection();
}

void PeriodMapWindow::updateSelection(bool replot)
{
  QVector<double> x, y;
  unsigned eligible = 0;
  for(const auto &site: sites) if(selectedStations.count(site.name)) {
    x.push_back(site.position.x()); y.push_back(site.position.y());
    if(observedIndex(site.name) >= 0) ++eligible;
  }
  if(selectionGraph) selectionGraph->setData(x, y);
  maskButton->setEnabled(eligible > 0); unmaskButton->setEnabled(eligible > 0);
  clearSelection->setEnabled(!selectedStations.empty());
  selectionSummary->setText(tr("%1 selected · %2 editable")
                            .arg(selectedStations.size()).arg(eligible));
  if(replot) plot->replot();
}

void PeriodMapWindow::setSelectedMasks(bool enabled)
{
  bool changed = false;
  for(const auto &name: selectedStations) {
    const int index = observedIndex(name);
    if(index < 0) continue;
    auto &station = survey->get_station_data(name);
    const double frequency = station.frequencies()[index];
    std::vector<RealDataType> types = maskGroup->currentIndex() == 0 ?
      std::vector<RealDataType>{RealTzx, RealTzy} : std::vector<RealDataType>{PTxx, PTxy, PTyx, PTyy};
    if(maskGroup->currentIndex() == 2) types.insert(types.end(), {RealZxx, RealZxy, RealZyx, RealZyy});
    for(auto type: types) station.set_data_mask(type, frequency, enabled, linkTensorMasks);
    changed = true;
  }
  if(!changed) return;
  if(observationsChanged) observationsChanged();
  updateMap();
}

void PeriodMapWindow::fitView()
{
  if(haveBounds) {
    const double span = std::max({boundsX.size(), boundsY.size(), ellipseSize->value(), arrowSize->value(), .1});
    plot->xAxis->setRange(boundsX.lower - .08 * span, boundsX.upper + .08 * span);
    plot->yAxis->setRange(boundsY.lower - .08 * span, boundsY.upper + .08 * span);
  } else {
    plot->xAxis->setRange(-1., 1.); plot->yAxis->setRange(-1., 1.);
  }
  plot->replot();
}

PeriodMapWindow::DisplaySettings PeriodMapWindow::displaySettings() const
{
  DisplaySettings s;
  s.period = period->currentData().isValid() ? period->currentData().toDouble() : 1.;
  s.tolerancePercent = tolerance->value(); s.ellipseKm = ellipseSize->value(); s.arrowKm = arrowSize->value();
  s.colorMin = colorMin->value(); s.colorMax = colorMax->value(); s.convention = convention->currentIndex();
  s.mapLayers = background->layers();
  s.colorBy = colorBy->currentIndex(); s.colormap = colorMap->currentText(); s.names = names->isChecked();
  s.real = realArrows->isChecked(); s.imaginary = imagArrows->isChecked(); s.reverse = reverseColors->isChecked();
  return s;
}

void PeriodMapWindow::setDisplaySettings(const DisplaySettings &s)
{
  const QList<QObject *> controls{tolerance, ellipseSize, arrowSize, colorMin, colorMax, convention, colorBy, colorMap, names, realArrows, imagArrows, reverseColors};
  std::vector<std::unique_ptr<QSignalBlocker>> blockers;
  for(auto *control: controls) blockers.emplace_back(new QSignalBlocker(control));
  background->setLayers(s.mapLayers);
  for(auto *action: findChild<QToolButton *>("mapLayers")->menu()->actions()) if(action->isCheckable()) {
    const QSignalBlocker block(action); action->setChecked(s.mapLayers & action->data().toInt());
  }
  tolerance->setValue(s.tolerancePercent); ellipseSize->setValue(s.ellipseKm); arrowSize->setValue(s.arrowKm);
  colorMin->setValue(s.colorMin); colorMax->setValue(s.colorMax); convention->setCurrentIndex(s.convention);
  colorBy->setCurrentIndex(s.colorBy); colorMap->setCurrentText(s.colormap); names->setChecked(s.names);
  realArrows->setChecked(s.real); imagArrows->setChecked(s.imaginary); reverseColors->setChecked(s.reverse);
  updateMap(true);
}

QString PeriodMapWindow::renderReport(QCPPainter &painter, const QRect &rect, double seconds,
                                    bool phaseTensors, bool induction, const QString &responseKey)
{
  const auto s = displaySettings();
  {
    const QSignalBlocker dataBlock(dataset), periodBlock(period), tensorBlock(tensors), realBlock(realArrows), imagBlock(imagArrows);
    dataset->setCurrentIndex(std::max(0, dataset->findData(responseKey)));
    int index = period->findData(seconds);
    if(index < 0) { period->addItem(QString::number(seconds, 'g', 12), seconds); index = period->count() - 1; }
    period->setCurrentIndex(index); tensors->setChecked(phaseTensors);
    realArrows->setChecked(induction && s.real); imagArrows->setChecked(induction && s.imaginary);
  }
  updateMap(true);
  auto reportFont = [](int size) { QFont f("Sans Serif"); f.setPixelSize(size); return f; };
  note->setText(coordinatesDescription + "\n" + originDescription + (induction ? " | " + convention->currentText() : QString()));
  for(auto *axis: {plot->xAxis, plot->yAxis, colorScale->axis()}) {
    axis->setLabelFont(reportFont(16)); axis->setTickLabelFont(reportFont(14));
  }
  plot->legend->setFont(reportFont(14)); title->setFont(reportFont(18)); note->setFont(reportFont(13));
  painter.save(); painter.translate(rect.topLeft()); plot->toPainter(&painter, rect.width(), rect.height()); painter.restore();
  // Layer selection is per report page; retain real/imaginary choices for the next page.
  const QSignalBlocker realBlock(realArrows), imagBlock(imagArrows);
  realArrows->setChecked(s.real); imagArrows->setChecked(s.imaginary);
  return summary->text();
}

void PeriodMapWindow::savePdf()
{
  auto path = QFileDialog::getSaveFileName(this, tr("Save period map"), "period_map.pdf", tr("PDF (*.pdf)"));
  if(path.isEmpty()) return;
  if(!path.endsWith(".pdf", Qt::CaseInsensitive)) path += ".pdf";
  const auto x = plot->xAxis->range(), y = plot->yAxis->range();
  if(!plot->savePdf(path, 1100, 900)) QMessageBox::warning(this, tr("Save period map"), tr("Cannot write the PDF file."));
  plot->xAxis->setRange(x); plot->yAxis->setRange(y); plot->replot();
}
