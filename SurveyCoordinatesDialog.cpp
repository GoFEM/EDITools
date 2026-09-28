#include "SurveyCoordinatesDialog.h"
#include "include/MTSurveyData.h"
#include "include/ToolTips.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <set>
#include <stdexcept>

namespace {
QString zoneName(int zone, bool north) {
  return QString::number(zone) + (north ? "N" : "S");
}
QString errorText(double fraction) {
  return QString::number(fraction * 100., 'f', 4);
}
} // namespace
SurveyCoordinatesDialog::SurveyCoordinatesDialog(MTSurveyData &data, QWidget *parent)
    : QDialog(parent), survey(data), geographic(data.geographic_locations()),
      bounds(UTM::extent(geographic)) {
  setObjectName("utmDialog");
  setWindowTitle(tr("Survey UTM coordinates"));
  resize(1040, 780);
  suggested = SurveyCoordinates::suggested(geographic);
  const auto settings = survey.coordinates().utm ? survey.coordinates() : suggested;
  auto *layout = new QVBoxLayout(this);
  auto *controls = new QHBoxLayout;
  zone = new QSpinBox(this);
  zone->setRange(1, 60);
  zone->setValue(settings.zone);
  zone->setObjectName("utmZone");
  hemisphere = new QComboBox(this);
  hemisphere->addItems({tr("North"), tr("South")});
  hemisphere->setCurrentIndex(settings.north ? 0 : 1);
  hemisphere->setObjectName("utmHemisphere");
  controls->addWidget(
      ToolTips::label(zone, tr("Zone"),
                      tr("One WGS84/UTM projection for every station and DEM. Choose a zone here, click a "
                         "comparison row, or double-click a zone on the map. Crossing a nominal zone "
                         "boundary does not introduce a coordinate discontinuity.")));
  controls->addWidget(zone);
  controls->addWidget(hemisphere);
  auto *automatic =
      new QPushButton(tr("Survey centre (%1)").arg(zoneName(suggested.zone, suggested.north)), this);
  automatic->setObjectName("utmSuggested");
  controls->addWidget(automatic);
  ToolTips::set(automatic, tr("Use the zone containing the midpoint of the survey's latitude and longitude "
                              "extent, using the shortest longitude arc across the date line. All stations "
                              "define the extent, including disabled ones. Station density does not bias "
                              "this default; an existing chosen projection is retained when reopening."));
  controls->addStretch();
  center = new QCheckBox(tr("Center UTM"), this);
  center->setObjectName("utmCenter");
  center->setChecked(settings.centered);
  controls->addWidget(center);
  ToolTips::set(center, tr("Subtract the midpoint of the projected station bounds from all map/export "
                           "coordinates. This changes the coordinate origin only, not projection distortion. "
                           "Original latitude, longitude and elevation are retained."));
  layout->addLayout(controls);
  crossing = new QLabel(this);
  crossing->setObjectName("utmCrossing");
  crossing->setWordWrap(true);
  layout->addWidget(crossing);
  auto *splitter = new QSplitter(Qt::Horizontal, this);
  layout->addWidget(splitter, 1);
  auto *mapPanel = new QWidget(splitter);
  auto *mapLayout = new QVBoxLayout(mapPanel);
  mapLayout->setContentsMargins(0, 0, 0, 0);
  map = new QCustomPlot(mapPanel);
  map->setObjectName("utmZoneMap");
  map->setMinimumSize(430, 270);
  mapLayout->addWidget(map, 1);
  map->addLayer("zones", map->layer("grid"), QCustomPlot::limBelow);
  map->xAxis->setLabel(tr("Longitude (°)"));
  map->yAxis->setLabel(tr("Latitude (°)"));
  map->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  map->axisRect()->setupFullAxesBox();
  map->legend->setVisible(true);
  map->axisRect()->insetLayout()->take(map->legend);
  map->plotLayout()->insertRow(0);
  map->plotLayout()->addElement(0, 0, map->legend);
  map->plotLayout()->setRowStretchFactor(0, .001);
  map->legend->setFillOrder(QCPLayoutGrid::foColumnsFirst);
  map->legend->setWrap(2);
  map->legend->setBorderPen(Qt::NoPen);
  for (const auto &color : {QColor("#2563a6"), QColor("#c66b21")}) {
    auto *g = map->addGraph();
    g->setLineStyle(QCPGraph::lsNone);
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, color, color, 7));
    g->setSelectable(QCP::stNone);
  }
  map->graph(0)->setName(tr("Selected zone"));
  map->graph(1)->setName(tr("Outside zone"));
  auto *mapControls = new QHBoxLayout;
  mapLayout->addLayout(mapControls);
  names = new QCheckBox(tr("Station names"), this);
  names->setObjectName("utmStationNames");
  mapControls->addWidget(names);
  auto *fit = new QPushButton(tr("Fit survey"), this);
  fit->setObjectName("utmFitSurvey");
  mapControls->addWidget(fit);
  mapControls->addStretch();
  auto *mapHint = new QLabel(tr("Double-click a zone to choose it"), this);
  mapControls->addWidget(mapHint);
  ToolTips::set(
      map,
      tr("Geographic overview of nominal UTM zones, including Norway and Svalbard exceptions. The selected "
         "zone is blue; its central meridian is dashed. Longitude can extend beyond ±180° to keep date-line "
         "surveys continuous. Pan/zoom to inspect boundaries; double-click to choose a zone."));
  auto *right = new QWidget(splitter);
  auto *rightLayout = new QVBoxLayout(right);
  rightLayout->setContentsMargins(0, 0, 0, 0);
  candidates = new QTableWidget(this);
  candidates->setObjectName("utmCandidates");
  candidates->setColumnCount(3);
  candidates->setHorizontalHeaderLabels({tr("Zone"), tr("Stations (%)"), tr("Extent (%)")});
  candidates->verticalHeader()->hide();
  candidates->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  candidates->setEditTriggers(QAbstractItemView::NoEditTriggers);
  candidates->setSelectionBehavior(QAbstractItemView::SelectRows);
  candidates->setSelectionMode(QAbstractItemView::SingleSelection);
  rightLayout->addWidget(candidates, 1);
  ToolTips::set(
      candidates,
      tr("Maximum absolute local scale error for each candidate projection. Stations uses the actual station "
         "locations. Extent samples a 33 × 33 grid over the geographic bounding rectangle plus the stations; "
         "it is an estimate, not a rigorous bound. 0.1% equals 1 metre per kilometre. Scale refers to "
         "distances on the reference ellipsoid, not elevation-adjusted ground distances. Candidates cover "
         "the survey extent, neighbouring zones, and your current choice."));
  auto *threshold = new QHBoxLayout;
  rightLayout->addLayout(threshold);
  limit = new QDoubleSpinBox(this);
  limit->setObjectName("utmDistortionLimit");
  limit->setDecimals(3);
  limit->setRange(.001, 100.);
  limit->setValue(.1);
  limit->setSuffix(" %");
  threshold->addWidget(ToolTips::label(
      limit, tr("Reference"),
      tr("Highlight projections whose sampled maximum scale error exceeds this reference. The default 0.1% "
         "is 1 metre per kilometre. This is a comparison aid, not an accuracy guarantee or an application "
         "limit; you can still choose a valid projection.")));
  threshold->addWidget(limit);
  splitter->setStretchFactor(0, 3);
  splitter->setStretchFactor(1, 2);
  splitter->setSizes({640, 360});
  distortion = new QLabel(this);
  distortion->setObjectName("utmDistortion");
  distortion->setWordWrap(true);
  layout->addWidget(distortion);
  ToolTips::set(distortion,
                tr("Scale error is grid distance divided by ellipsoidal distance minus one. Positive values "
                   "stretch distance; negative values shorten it. Convergence is the angle between true "
                   "north and grid north at a station. Coordinate centering does not change either quantity. "
                   "Impedance/tipper components are not rotated by this dialog."));
  auto *origin = new QHBoxLayout;
  layout->addLayout(origin);
  easting = new QLineEdit(this);
  northing = new QLineEdit(this);
  easting->setObjectName("utmOriginEasting");
  northing->setObjectName("utmOriginNorthing");
  easting->setReadOnly(true);
  northing->setReadOnly(true);
  origin->addWidget(new QLabel(tr("Origin E (m)"), this));
  origin->addWidget(easting);
  origin->addWidget(new QLabel(tr("N (m)"), this));
  origin->addWidget(northing);
  copy = new QPushButton(tr("Copy origin"), this);
  copy->setObjectName("copyUtmOrigin");
  origin->addWidget(copy);
  table = new QTableWidget(this);
  table->setObjectName("utmStations");
  table->setColumnCount(8);
  table->setHorizontalHeaderLabels({tr("Station"), tr("Zone"), tr("Easting (m)"), tr("Northing (m)"),
                                    tr("X / North (m)"), tr("Y / East (m)"), tr("Scale (m/km)"),
                                    tr("Convergence (°)")});
  table->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  table->horizontalHeader()->setStretchLastSection(true);
  table->setMinimumHeight(125);
  table->setMaximumHeight(210);
  layout->addWidget(table);
  ToolTips::set(table,
                tr("Zone is the station's nominal zone and hemisphere. Easting/Northing use the chosen "
                   "common projection. X/Y subtract the displayed origin. Scale gives signed local distance "
                   "error in metres per kilometre; convergence gives true-north/grid-north angle in degrees. "
                   "Values do not include terrain-height corrections."));
  status = new QLabel(this);
  status->setObjectName("utmStatus");
  status->setWordWrap(true);
  layout->addWidget(status);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
  apply = buttons->button(QDialogButtonBox::Ok);
  apply->setText(tr("Use these coordinates"));
  layout->addWidget(buttons);
  connect(zone, qOverload<int>(&QSpinBox::valueChanged), this, [this] { updatePreview(); });
  connect(hemisphere, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { updatePreview(); });
  connect(center, &QCheckBox::toggled, this, [this] { updatePreview(); });
  connect(limit, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this] { updatePreview(); });
  connect(automatic, &QPushButton::clicked, this, [this] {
    const QSignalBlocker a(zone), b(hemisphere);
    zone->setValue(suggested.zone);
    hemisphere->setCurrentIndex(suggested.north ? 0 : 1);
    updatePreview();
  });
  connect(candidates, &QTableWidget::cellClicked, this,
          [this](int row, int) { zone->setValue(candidates->item(row, 0)->data(Qt::UserRole).toInt()); });
  connect(fit, &QPushButton::clicked, this, [this] { fitSurvey(); });
  connect(names, &QCheckBox::toggled, this, [this] {
    updateZoneOverlay();
    map->replot();
  });
  connect(map->xAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this,
          [this] { updateZoneOverlay(); });
  connect(map->yAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this,
          [this] { updateZoneOverlay(); });
  connect(map, &QCustomPlot::mouseDoubleClick, this, [this](QMouseEvent *event) {
    if (event->button() != Qt::LeftButton || !map->axisRect()->rect().contains(event->pos()))
      return;
    const double lat = map->yAxis->pixelToCoord(event->pos().y()),
                 lon = map->xAxis->pixelToCoord(event->pos().x());
    if (lat < -80. || lat > 84.)
      return;
    const QSignalBlocker a(zone), b(hemisphere);
    zone->setValue(UTM::zoneAt(lat, lon));
    hemisphere->setCurrentIndex(lat >= 0. ? 0 : 1);
    updatePreview();
  });
  connect(copy, &QPushButton::clicked, this,
          [this] { QApplication::clipboard()->setText(QString::fromStdString(preview.description())); });
  connect(buttons, &QDialogButtonBox::accepted, this, [this] {
    survey.set_coordinates(preview);
    accept();
  });
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  fitSurvey();
  updatePreview();
}

const UTM::Analysis &SurveyCoordinatesDialog::analysis(int selectedZone) {
  const int key = selectedZone + (hemisphere->currentIndex() ? 60 : 0);
  auto found = analyses.find(key);
  if (found == analyses.end()) {
    SurveyCoordinates c;
    c.utm = true;
    c.zone = selectedZone;
    c.north = hemisphere->currentIndex() == 0;
    found = analyses.emplace(key, UTM::analyze(geographic, c)).first;
  }
  return found->second;
}
void SurveyCoordinatesDialog::fitSurvey() {
  const double dx = std::max(4., (bounds.east - bounds.west) * .6),
               dy = std::max(.5, (bounds.north - bounds.south) * .6);
  map->xAxis->setRange(bounds.longitude() - dx, bounds.longitude() + dx);
  map->yAxis->setRange(std::max(-80., bounds.latitude() - dy), std::min(84., bounds.latitude() + dy));
  updateZoneOverlay();
  map->replot();
}
void SurveyCoordinatesDialog::updateZoneOverlay() {
  map->clearItems();
  const auto x = map->xAxis->range(), y = map->yAxis->range();
  // Limit overlays while zoomed far beyond the globe.
  UTM::Extent view{std::max(x.lower, bounds.longitude() - 180.), std::min(x.upper, bounds.longitude() + 180.),
                   std::max(-80., y.lower), std::min(84., y.upper)};
  if (view.east <= view.west || view.north <= view.south)
    return;
  for (const auto &r : UTM::regions(view)) {
    const bool selected = r.zone == zone->value() && r.north == (hemisphere->currentIndex() == 0);
    auto *area = new QCPItemRect(map);
    area->setObjectName(QString("utmRegion%1%2").arg(r.zone).arg(r.north ? 'N' : 'S'));
    area->setLayer("zones");
    area->setSelectable(false);
    area->topLeft->setCoords(r.west, r.top);
    area->bottomRight->setCoords(r.east, r.south);
    area->setBrush(selected ? QColor(72, 141, 203, 65) : QColor(115, 130, 145, r.zone % 2 ? 15 : 30));
    area->setPen(QPen(QColor("#8b98a5"), .8));
    auto *label = new QCPItemText(map);
    label->setLayer("zones");
    label->setSelectable(false);
    label->position->setCoords((r.west + r.east) / 2., r.top - .06 * (r.top - r.south));
    label->setText(zoneName(r.zone, r.north));
    label->setColor(selected ? QColor("#18588e") : QColor("#596775"));
    label->setFont(QFont(font().family(), 10, QFont::Bold));
    label->setBrush(QColor(255, 255, 255, 170));
  }
  auto *meridian = new QCPItemLine(map);
  meridian->setLayer("zones");
  meridian->setSelectable(false);
  const double lon = bounds.unwrap(zone->value() * 6. - 183.);
  meridian->start->setCoords(lon, -80.);
  meridian->end->setCoords(lon, 84.);
  meridian->setPen(QPen(QColor("#3478ad"), 1, Qt::DashLine));
  auto *middle = new QCPItemText(map);
  middle->setLayer("zones");
  middle->setSelectable(false);
  middle->position->setCoords(bounds.longitude(), bounds.latitude());
  middle->setText("+");
  middle->setColor(QColor("#c63d44"));
  middle->setFont(QFont(font().family(), 16, QFont::Bold));
  if (names->isChecked()) {
    const auto stationNames = survey.get_stations_names();
    for (size_t i = 0; i < geographic.size(); ++i) {
      auto *label = new QCPItemText(map);
      label->setSelectable(false);
      label->position->setCoords(bounds.unwrap(geographic[i][1]), geographic[i][0]);
      label->setText(QString::fromStdString(stationNames[i]));
      label->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom);
      label->setPadding(QMargins(5, 0, 0, 3));
    }
  }
}
void SurveyCoordinatesDialog::updatePreview() {
  try {
    std::set<int> choices{suggested.zone, zone->value(), suggested.zone == 1 ? 60 : suggested.zone - 1,
                          suggested.zone == 60 ? 1 : suggested.zone + 1};
    auto cover = bounds;
    cover.west -= 1e-8;
    cover.east += 1e-8;
    cover.south = std::max(-80., cover.south - 1e-8);
    cover.north = std::min(84., cover.north + 1e-8);
    for (const auto &r : UTM::regions(cover))
      choices.insert(r.zone);
    const auto &a = analysis(zone->value());
    QStringList crossed;
    for (const auto &entry : a.stationZones)
      crossed << QString::number(entry.first);
    crossing->setText(tr("Station zones: %1 · %2 outside selected zone%3")
                          .arg(crossed.join(", "))
                          .arg(a.outsideZone)
                          .arg(a.oppositeHemisphere
                                   ? tr(" · %1 in other hemisphere").arg(a.oppositeHemisphere)
                                   : QString()));
    std::vector<int> ordered(choices.begin(), choices.end());
    std::sort(ordered.begin(), ordered.end(),
              [&](int l, int r) { return bounds.unwrap(l * 6. - 183.) < bounds.unwrap(r * 6. - 183.); });
    const QSignalBlocker block(candidates);
    candidates->setRowCount(int(ordered.size()));
    for (size_t i = 0; i < ordered.size(); ++i) {
      const int z = ordered[i];
      const auto &candidate = analysis(z);
      const QStringList cells{zoneName(z, hemisphere->currentIndex() == 0) +
                                  (z == suggested.zone ? tr(" · centre") : QString()),
                              candidate.valid() ? errorText(candidate.stationMaxError) : tr("N/A"),
                              candidate.valid() ? errorText(candidate.extentMaxError) : tr("N/A")};
      for (int c = 0; c < 3; ++c) {
        auto *item = new QTableWidgetItem(cells[c]);
        item->setData(Qt::UserRole, z);
        if (!candidate.valid() || candidate.extentMaxError * 100. > limit->value())
          item->setForeground(QColor("#a45216"));
        candidates->setItem(i, c, item);
      }
      if (z == zone->value())
        candidates->selectRow(i);
    }
    QVector<double> x[2], y[2];
    for (const auto &p : geographic) {
      const int group =
          UTM::zoneAt(p[0], p[1]) == zone->value() && ((p[0] >= 0.) == (hemisphere->currentIndex() == 0)) ? 0
                                                                                                          : 1;
      x[group].push_back(bounds.unwrap(p[1]));
      y[group].push_back(p[0]);
    }
    for (int i = 0; i < 2; ++i)
      map->graph(i)->setData(x[i], y[i]);
    updateZoneOverlay();
    map->replot();
    if (!a.valid())
      throw std::runtime_error(
          "Selected zone cannot represent the complete survey extent. Choose another zone.");
    preview = SurveyCoordinates::calculate(geographic, zone->value(), hemisphere->currentIndex() == 0,
                                           center->isChecked());
    easting->setText(QString::number(preview.origin_easting, 'g', 17));
    northing->setText(QString::number(preview.origin_northing, 'g', 17));
    const auto positions = preview.transform(geographic);
    const auto stationNames = survey.get_stations_names();
    table->setRowCount(int(stationNames.size()));
    for (size_t i = 0; i < stationNames.size(); ++i) {
      const auto &p = positions[i];
      const auto &f = a.stations[i];
      QStringList cells{QString::fromStdString(stationNames[i]),
                        zoneName(UTM::zoneAt(geographic[i][0], geographic[i][1]), geographic[i][0] >= 0.)};
      for (double value : {p[1] + preview.origin_easting, p[0] + preview.origin_northing, p[0], p[1],
                           (f.scale - 1.) * 1000., f.convergence})
        cells << QString::number(value, 'f', 3);
      for (int c = 0; c < cells.size(); ++c)
        table->setItem(i, c, new QTableWidgetItem(cells[c]));
    }
    distortion->setText(tr("Max scale error: stations %1% · sampled extent %2% (%3 m/km)\nSigned range: %4 "
                           "to %5 m/km · max |convergence|: %6°")
                            .arg(errorText(a.stationMaxError), errorText(a.extentMaxError))
                            .arg(a.extentMaxError * 1000., 0, 'f', 3)
                            .arg((a.minimumScale - 1.) * 1000., 0, 'f', 3)
                            .arg((a.maximumScale - 1.) * 1000., 0, 'f', 3)
                            .arg(a.maxConvergence, 0, 'f', 3));
    const bool above = a.extentMaxError * 100. > limit->value();
    status->setStyleSheet(above ? "color: #a45216" : "");
    status->setText(
        above ? tr("Above the %1% reference. You can compare zones or keep this choice.")
                    .arg(limit->value(), 0, 'f', 3)
              : tr("Within the %1% reference at sampled locations.").arg(limit->value(), 0, 'f', 3));
    apply->setEnabled(true);
    copy->setEnabled(true);
  } catch (const std::exception &e) {
    status->setText(QString::fromUtf8(e.what()));
    status->setStyleSheet("color: #a45216");
    distortion->clear();
    apply->setEnabled(false);
    copy->setEnabled(false);
    table->setRowCount(0);
    easting->clear();
    northing->clear();
  }
}
