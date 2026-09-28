#include "TopographyDialog.h"
#include "include/PlotColorMaps.h"
#include <QCheckBox>
#include <QComboBox>
#include <QCryptographicHash>
#include <QDate>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QEventLoop>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QProgressDialog>
#include <QPushButton>
#include <QSaveFile>
#include <QStandardPaths>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <QUrlQuery>
#include <QVBoxLayout>
#include <QApplication>
#include <cmath>
#include <stdexcept>

namespace {
QString text(const std::string &s) {
  return QString::fromStdString(s);
}
QString num(double v) {
  return QString::number(v, 'g', 17);
}
QString cachedDownload(const QUrl &url, QWidget *parent, QProgressDialog &progress, int attempt = 0) {
  const auto dir = QDir(QStandardPaths::writableLocation(QStandardPaths::CacheLocation)).filePath("dem");
  if (!QDir().mkpath(dir))
    throw std::runtime_error("Cannot create DEM cache.");
  const auto path = QDir(dir).filePath(
      QString::fromLatin1(QCryptographicHash::hash(url.toEncoded(), QCryptographicHash::Sha256).toHex()) +
      ".tif");
  if (QFileInfo(path).size() > 0)
    return path;
  QSaveFile file(path);
  if (!file.open(QIODevice::WriteOnly))
    throw std::runtime_error("Cannot write DEM cache.");
  QNetworkAccessManager network(parent);
  QNetworkRequest request(url);
  request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::NoLessSafeRedirectPolicy);
  request.setRawHeader("User-Agent", "EDITools-DEM/1.0 (+https://github.com/GoFEM/EDITools)");
  auto *reply = network.get(request);
  QEventLoop loop;
  QTimer timeout;
  timeout.setSingleShot(true);
  bool failed = false;
  qint64 bytes = 0;
  QByteArray prefix;
  QObject::connect(reply, &QNetworkReply::readyRead, &loop, [&] {
    const auto chunk = reply->readAll();
    if (prefix.size() < 512)
      prefix += chunk.left(512 - prefix.size());
    bytes += chunk.size();
    if (bytes > 512ll * 1024 * 1024 || file.write(chunk) != chunk.size()) {
      failed = true;
      reply->abort();
    }
  });
  QObject::connect(reply, &QNetworkReply::downloadProgress, &loop, [&](qint64 received, qint64 total) {
    progress.setLabelText(QObject::tr("Downloading DEM: %1 / %2 MB")
                              .arg(received / 1048576.)
                              .arg(total > 0 ? QString::number(total / 1048576., 'f', 1) : QString("?")));
  });
  QObject::connect(&progress, &QProgressDialog::canceled, reply, &QNetworkReply::abort);
  QObject::connect(&timeout, &QTimer::timeout, reply, &QNetworkReply::abort);
  QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
  timeout.start(180000);
  loop.exec();
  if (progress.wasCanceled())
    throw std::runtime_error("DEM download canceled; survey unchanged.");
  const int status = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
  if ((status == 502 || status == 503 || status == 504) && attempt < 2) {
    file.cancelWriting();
    progress.setLabelText(QObject::tr("Provider busy; retrying DEM request…"));
    QEventLoop pause;
    QTimer::singleShot(1000, &pause, &QEventLoop::quit);
    pause.exec();
    if (progress.wasCanceled())
      throw std::runtime_error("DEM download canceled; survey unchanged.");
    return cachedDownload(url, parent, progress, attempt + 1);
  }
  if (failed || reply->error() != QNetworkReply::NoError ||
      reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt() != 200)
    throw std::runtime_error(QString("DEM download failed (HTTP %1): %2")
                                 .arg(reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt())
                                 .arg(reply->errorString())
                                 .toStdString());
  const bool tiff =
      prefix.startsWith(QByteArray("II\x2a\0", 4)) || prefix.startsWith(QByteArray("MM\0\x2a", 4)) ||
      prefix.startsWith(QByteArray("II\x2b\0", 4)) || prefix.startsWith(QByteArray("MM\0\x2b", 4));
  if (!tiff)
    throw std::runtime_error("Provider returned no elevation raster. Reduce the area or check coverage.");
  if (!file.commit())
    throw std::runtime_error("Cannot finish DEM cache file.");
  return path;
}
void publish(const QString &temporary, const QString &destination) {
  QFile in(temporary);
  QSaveFile out(destination);
  if (!in.open(QIODevice::ReadOnly) || !out.open(QIODevice::WriteOnly))
    throw std::runtime_error("Cannot open DEM export destination.");
  while (!in.atEnd()) {
    const auto data = in.read(1024 * 1024);
    if (data.isEmpty() && in.error() != QFile::NoError)
      throw std::runtime_error("Cannot read DEM export.");
    if (out.write(data) != data.size())
      throw std::runtime_error("Cannot write DEM export.");
  }
  if (!out.commit())
    throw std::runtime_error("Cannot finish DEM export.");
}
void sampleStations(QWidget *parent, MTSurveyData &survey, const std::function<void()> &changed) {
  const auto grid = survey.topography().grid;
  if (!grid)
    return;
  QDialog preview(parent);
  preview.setObjectName("demStationPreview");
  preview.setWindowTitle(QObject::tr("Station elevations"));
  preview.resize(680, 420);
  auto *layout = new QVBoxLayout(&preview);
  auto *label =
      new QLabel(QObject::tr("DEM datum: %1 · elevations in metres").arg(text(grid->datum)), &preview);
  label->setWordWrap(true);
  layout->addWidget(label);
  auto *table = new QTableWidget(&preview);
  table->setObjectName("demStationSamples");
  table->setColumnCount(4);
  table->setHorizontalHeaderLabels({"Station", "Current (m)", "DEM (m)", "Difference (m)"});
  table->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  layout->addWidget(table);
  const auto names = survey.get_stations_names();
  std::vector<double> x, y;
  for (const auto &name : names) {
    const auto p = survey.get_station_data(name).position();
    x.push_back(p[1]);
    y.push_back(p[0]);
  }
  Topography::transform(4326, grid->epsg, x, y);
  table->setRowCount(int(names.size()));
  std::vector<double> elevations;
  unsigned count = 0;
  for (size_t i = 0; i < names.size(); ++i) {
    const double old = survey.get_station_data(names[i]).position()[2], z = grid->sample(x[i], y[i]);
    elevations.push_back(z);
    auto *item = new QTableWidgetItem(text(names[i]));
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(std::isfinite(z) ? Qt::Checked : Qt::Unchecked);
    if (!std::isfinite(z))
      item->setFlags(item->flags() & ~Qt::ItemIsEnabled);
    else
      ++count;
    table->setItem(i, 0, item);
    for (int c = 1; c < 4; ++c) {
      const double v = c == 1 ? old : c == 2 ? z : z - old;
      table->setItem(i, c,
                     new QTableWidgetItem(std::isfinite(v) ? QString::number(v, 'f', 2) : QString("—")));
    }
  }
  auto *datum = new QCheckBox(QObject::tr("Vertical datum checked"), &preview);
  datum->setObjectName("demDatumChecked");
  datum->setToolTip(
      QObject::tr("Confirm the DEM heights use the vertical reference required by your survey. No geoid or "
                  "vertical datum conversion is performed. Applying replaces only checked, covered station "
                  "elevations; transfer functions and masks are unchanged."));
  layout->addWidget(datum);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Apply | QDialogButtonBox::Cancel, &preview);
  layout->addWidget(buttons);
  buttons->button(QDialogButtonBox::Apply)->setEnabled(false);
  QObject::connect(datum, &QCheckBox::toggled, &preview,
                   [=](bool on) { buttons->button(QDialogButtonBox::Apply)->setEnabled(on && count); });
  QObject::connect(buttons->button(QDialogButtonBox::Apply), &QPushButton::clicked, &preview,
                   &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &preview, &QDialog::reject);
  if (preview.exec() != QDialog::Accepted)
    return;
  for (size_t i = 0; i < names.size(); ++i)
    if (table->item(i, 0)->checkState() == Qt::Checked && std::isfinite(elevations[i])) {
      auto &s = survey.get_station_data(names[i]);
      auto p = s.position();
      p[2] = elevations[i];
      s.set_position(p);
    }
  changed();
}
} // namespace

TopographyDialog::TopographyDialog(QWidget *parent, MTSurveyData &survey, std::function<void()> changed)
    : QDialog(parent) {
  setObjectName("topographyDialog");
  setWindowTitle(tr("Topography"));
  resize(590, 620);
  auto *layout = new QVBoxLayout(this);
  auto *form = new QFormLayout;
  layout->addLayout(form);
  auto *source = new QComboBox(this);
  source->setObjectName("demSource");
  source->addItems({tr("Local GeoTIFF"), tr("Copernicus GLO-30"), tr("USGS 3DEP")});
  form->addRow(tr("Source"), source);
  source->setToolTip(tr("Local GeoTIFFs or explicit online downloads. Copernicus is a surface model "
                        "including vegetation/buildings; 3DEP provides US terrain. Downloaded data is cached "
                        "locally and the clipped grid is saved in the project."));
  auto spin = [&](const char *name, double min, double max, double value, int decimals = 1) {
    auto *s = new QDoubleSpinBox(this);
    s->setObjectName(name);
    s->setDecimals(decimals);
    s->setRange(min, max);
    s->setValue(value);
    return s;
  };
  auto *margin = spin("demMargin", 0, 500, 5);
  auto *spacing = spin("demSpacing", 1, 10000, 30);
  auto *areaRow = new QWidget(this);
  auto *areaLayout = new QHBoxLayout(areaRow);
  areaLayout->setContentsMargins(0, 0, 0, 0);
  areaLayout->addWidget(margin);
  areaLayout->addWidget(new QLabel(tr("km margin"), this));
  auto *fit = new QPushButton(tr("Survey extent"), this);
  fit->setObjectName("demSurveyExtent");
  areaLayout->addWidget(fit);
  form->addRow(tr("Area"), areaRow);
  form->addRow(tr("Cell size (m)"), spacing);
  spacing->setToolTip(tr("Requested output grid spacing, not a promise of source resolution. Coarser grids "
                         "reduce download, display and project size. At most 16 million cells."));
  const auto initial = Topography::surveyArea(survey.geographic_locations(), survey.coordinates(),
                                              margin->value() * 1000, spacing->value());
  auto *bounds = new QGroupBox(tr("UTM bounds · EPSG:%1").arg(initial.epsg), this);
  auto *boundsLayout = new QGridLayout(bounds);
  const char *ids[]{"demWest", "demSouth", "demEast", "demNorth"};
  const QStringList labels{tr("West"), tr("South"), tr("East"), tr("North")};
  std::array<QDoubleSpinBox *, 4> extent;
  const double values[]{initial.west, initial.south, initial.east, initial.north};
  for (int i = 0; i < 4; ++i) {
    extent[i] = spin(ids[i], -20000000, 20000000, values[i], 2);
    boundsLayout->addWidget(new QLabel(labels[i], this), i / 2, (i % 2) * 2);
    boundsLayout->addWidget(extent[i], i / 2, (i % 2) * 2 + 1);
  }
  bounds->setToolTip(tr("Absolute UTM metres. Downloads and exports use this area. XYZ exports subtract the "
                        "survey origin; GeoTIFF retains standard absolute coordinates."));
  layout->addWidget(bounds);
  auto *metadata = new QFormLayout;
  layout->addLayout(metadata);
  auto *units = new QComboBox(this);
  units->setObjectName("demUnits");
  units->addItems({tr("Metres"), tr("Feet"), tr("US survey feet"), tr("Read file units")});
  metadata->addRow(tr("Input units"), units);
  auto *datum = new QComboBox(this);
  datum->setObjectName("demDatum");
  datum->setEditable(true);
  datum->addItems({"Unknown", "EGM2008", "NAVD88", "WGS84 ellipsoid"});
  metadata->addRow(tr("Vertical datum"), datum);
  datum->setToolTip(tr("Declare the input height reference. Horizontal reprojection and unit conversion "
                       "preserve this datum; no vertical datum conversion is performed."));
  auto *credit = new QLineEdit(tr("User-provided DEM"), this);
  credit->setObjectName("demCredit");
  metadata->addRow(tr("Credit"), credit);
  credit->setToolTip(tr("Public source attribution included in maps and numeric exports. Enter a dataset "
                        "name or licence notice; avoid local paths or private information."));
  auto *buttons = new QHBoxLayout;
  layout->addLayout(buttons);
  auto button = [&](const QString &label, const char *name) {
    auto *b = new QPushButton(label, this);
    b->setObjectName(name);
    buttons->addWidget(b);
    return b;
  };
  auto *load = button(tr("Load…"), "demLoad"), *remove = button(tr("Remove"), "demRemove"),
       *sample = button(tr("Station heights…"), "demSampleStations");
  auto *licence = new QCheckBox(tr("Accept Copernicus DEM licence"), this);
  licence->setObjectName("demAcceptLicence");
  licence->hide();
  layout->addWidget(licence);
  auto *licenceLink = new QLabel(tr("<a "
                                    "href=\"https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/"
                                    "Data/DEM/resources/license/License-COPDEM-30.pdf\">Dataset licence</a>"),
                                 this);
  licenceLink->setOpenExternalLinks(true);
  licenceLink->hide();
  layout->addWidget(licenceLink);
  connect(licence, &QCheckBox::toggled, this,
          [=](bool on) { load->setEnabled(source->currentIndex() != 1 || on); });
  auto *display = new QGroupBox(tr("Display"), this);
  auto *displayLayout = new QGridLayout(display);
  layout->addWidget(display);
  auto check = [&](const QString &label, const char *id, bool value, int column) {
    auto *c = new QCheckBox(label, this);
    c->setObjectName(id);
    c->setChecked(value);
    displayLayout->addWidget(c, 0, column);
    return c;
  };
  const auto style = survey.topography().style;
  auto *visible = check(tr("Show"), "demVisible", style.visible, 0),
       *colors = check(tr("Colours"), "demColors", style.colors, 1),
       *shade = check(tr("Hillshade"), "demShade", style.shade, 2),
       *contours = check(tr("Contours"), "demContours", style.contours, 3);
  auto *palette = new QComboBox(this);
  palette->setObjectName("demPalette");
  palette->addItem("Terrain");
  palette->addItems(PlotColorMaps::names());
  palette->setCurrentText(text(style.palette));
  displayLayout->addWidget(palette, 1, 0, 1, 2);
  auto *opacity = spin("demOpacity", 0, 100, style.opacity * 100, 0);
  opacity->setSuffix(" %");
  displayLayout->addWidget(new QLabel(tr("Opacity"), this), 1, 2);
  displayLayout->addWidget(opacity, 1, 3);
  auto *lower = spin("demLower", -100000, 100000, style.lower),
       *upper = spin("demUpper", -100000, 100000, style.upper),
       *interval = spin("demContourInterval", 1, 10000, style.interval);
  displayLayout->addWidget(new QLabel(tr("Range (m)"), this), 2, 0);
  displayLayout->addWidget(lower, 2, 1);
  displayLayout->addWidget(upper, 2, 2);
  auto *autoRange = new QPushButton(tr("Auto"), this);
  displayLayout->addWidget(autoRange, 2, 3);
  displayLayout->addWidget(new QLabel(tr("Contour interval (m)"), this), 3, 0, 1, 2);
  displayLayout->addWidget(interval, 3, 2, 1, 2);
  auto *exports = new QHBoxLayout;
  layout->addLayout(exports);
  auto *tif = new QPushButton(tr("GeoTIFF…"), this);
  tif->setObjectName("demExportTiff");
  auto *xyz = new QPushButton(tr("XYZ…"), this);
  xyz->setObjectName("demExportXYZ");
  auto *sign = new QComboBox(this);
  sign->setObjectName("demHeightSign");
  sign->addItems({tr("XYZ: elevation (+up)"), tr("XYZ: depth (+down)")});
  exports->addWidget(tif);
  exports->addWidget(xyz);
  exports->addWidget(sign);
  auto *status = new QLabel(this);
  status->setObjectName("demStatus");
  status->setWordWrap(true);
  layout->addWidget(status);
  auto *close = new QDialogButtonBox(QDialogButtonBox::Close, this);
  layout->addWidget(close);
  connect(close, &QDialogButtonBox::rejected, this, &QDialog::reject);
  auto refresh = [=, &survey] {
    const auto g = survey.topography().grid;
    status->setToolTip(g ? text(g->notice) : QString());
    for (auto *b : {remove, sample, tif, xyz})
      b->setEnabled(bool(g));
    display->setEnabled(bool(g));
    if (!g)
      status->setText(tr("No DEM loaded."));
    else {
      size_t valid = 0;
      for (float z : g->values)
        if (std::isfinite(z))
          ++valid;
      status->setText(tr("%1 × %2 · %3 m cells · %4% coverage\n%5 · %6")
                          .arg(g->width)
                          .arg(g->height)
                          .arg(g->spacing)
                          .arg(100. * valid / g->values.size(), 0, 'f', 1)
                          .arg(text(g->source), text(g->datum)));
    }
  };
  auto saveStyle = [=, &survey] {
    if (upper->value() <= lower->value())
      return;
    auto &s = survey.topography().style;
    s.visible = visible->isChecked();
    s.colors = colors->isChecked();
    s.shade = shade->isChecked();
    s.contours = contours->isChecked();
    s.lower = lower->value();
    s.upper = upper->value();
    s.opacity = opacity->value() / 100.;
    s.interval = interval->value();
    s.palette = palette->currentText().toStdString();
    changed();
  };
  for (auto *c : {visible, colors, shade, contours})
    connect(c, &QCheckBox::toggled, this, saveStyle);
  for (auto *s : {lower, upper, opacity, interval})
    connect(s, qOverload<double>(&QDoubleSpinBox::valueChanged), this, saveStyle);
  connect(palette, &QComboBox::currentTextChanged, this, saveStyle);
  auto resetRange = [=, &survey] {
    if (!survey.topography().grid)
      return;
    const auto r = survey.topography().grid->range();
    const QSignalBlocker a(lower), b(upper);
    lower->setValue(r[0]);
    upper->setValue(r[1] > r[0] ? r[1] : r[0] + 1.);
    saveStyle();
  };
  connect(autoRange, &QPushButton::clicked, this, resetRange);
  auto error = [=](const std::exception &e) {
    QMessageBox::warning(this, tr("Topography"), QString::fromUtf8(e.what()));
  };
  connect(fit, &QPushButton::clicked, this, [=, &survey] {
    try {
      const auto r = Topography::surveyArea(survey.geographic_locations(), survey.coordinates(),
                                            margin->value() * 1000, spacing->value());
      extent[0]->setValue(r.west);
      extent[1]->setValue(r.south);
      extent[2]->setValue(r.east);
      extent[3]->setValue(r.north);
    } catch (const std::exception &e) {
      error(e);
    }
  });
  connect(source, qOverload<int>(&QComboBox::currentIndexChanged), this, [=](int index) {
    load->setText(index ? tr("Download") : tr("Load…"));
    units->setEnabled(!index);
    datum->setEnabled(!index);
    credit->setEnabled(!index);
    licence->setVisible(index == 1);
    licenceLink->setVisible(index == 1);
    load->setEnabled(index != 1 || licence->isChecked());
    units->setCurrentIndex(0);
    spacing->setValue(index == 2 ? 10 : 30);
    datum->setCurrentText(index == 1   ? "EGM2008"
                          : index == 2 ? "USGS 3DEP source datum (varies)"
                                       : "Unknown");
    credit->setText(index == 1   ? tr("Copernicus GLO-30 (2021)")
                    : index == 2 ? tr("USGS 3DEP · public domain")
                                 : tr("User-provided DEM"));
  });
  auto request = [=] {
    auto r = initial;
    r.west = extent[0]->value();
    r.south = extent[1]->value();
    r.east = extent[2]->value();
    r.north = extent[3]->value();
    r.spacing = spacing->value();
    const double factors[]{1., .3048, 1200. / 3937., 0.};
    r.unitScale = factors[units->currentIndex()];
    r.source = credit->text().trimmed().toStdString();
    r.datum = datum->currentText().trimmed().toStdString();
    r.validate();
    return r;
  };
  connect(load, &QPushButton::clicked, this, [=, &survey] {
    try {
      auto r = request();
      QStringList paths;
      if (!source->currentIndex()) {
        paths = QFileDialog::getOpenFileNames(this, tr("Load DEM"), {}, tr("GeoTIFF (*.tif *.tiff)"));
        if (paths.empty())
          return;
      }
      QProgressDialog progress(tr("Preparing DEM…"), tr("Cancel"), 0, 100, this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0);
      progress.setAutoClose(false);
      progress.show();
      if (source->currentIndex() == 1) {
        r.notice =
            "produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH "
            "2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved. The "
            "organisations in charge of the Copernicus programme by law or by delegation do not incur any "
            "liability for any use of the Copernicus WorldDEM-30. Licence: "
            "https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Data/DEM/resources/license/"
            "License-COPDEM-30.pdf";
        const auto b = Topography::geographicBounds(r);
        const int west = int(std::floor(b[0])), east = int(std::ceil(b[2])), south = int(std::floor(b[1])),
                  north = int(std::ceil(b[3]));
        if ((east - west) * (north - south) > 16)
          throw std::runtime_error("Download at most 16 Copernicus tiles at once; reduce the area.");
        for (int lat = south; lat < north; ++lat)
          for (int lon = west; lon < east; ++lon) {
            if (progress.wasCanceled())
              throw std::runtime_error("DEM download canceled; survey unchanged.");
            const QString tile = QString("Copernicus_DSM_COG_10_%1%2_00_%3%4_00_DEM")
                                     .arg(lat < 0 ? 'S' : 'N')
                                     .arg(std::abs(lat), 2, 10, QChar('0'))
                                     .arg(lon < 0 ? 'W' : 'E')
                                     .arg(std::abs(lon), 3, 10, QChar('0'));
            paths << cachedDownload(
                QUrl("https://copernicus-dem-30m.s3.amazonaws.com/" + tile + "/" + tile + ".tif"), this,
                progress);
          }
      } else if (source->currentIndex() == 2) {
        const int w = int(std::ceil((r.east - r.west) / r.spacing)),
                  h = int(std::ceil((r.north - r.south) / r.spacing));
        if (w > 8000 || h > 8000)
          throw std::runtime_error("USGS downloads allow at most 8000 rows or columns; increase cell size.");
        QUrl url(
            "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage");
        QUrlQuery q;
        q.addQueryItem("f", "image");
        q.addQueryItem("bbox", QStringList{num(r.west), num(r.north - h * r.spacing),
                                           num(r.west + w * r.spacing), num(r.north)}
                                   .join(','));
        q.addQueryItem("bboxSR", QString::number(r.epsg));
        q.addQueryItem("imageSR", QString::number(r.epsg));
        q.addQueryItem("size", QString("%1,%2").arg(w).arg(h));
        q.addQueryItem("format", "tiff");
        q.addQueryItem("pixelType", "F32");
        q.addQueryItem("noData", "-999999");
        q.addQueryItem("adjustAspectRatio", "false");
        q.addQueryItem("interpolation", "RSP_BilinearInterpolation");
        q.addQueryItem("renderingRule", "{\"rasterFunction\":\"None\"}");
        QByteArray encoded;
        for (const auto &item : q.queryItems()) {
          if (!encoded.isEmpty())
            encoded += '&';
          encoded += QUrl::toPercentEncoding(item.first) + '=' + QUrl::toPercentEncoding(item.second);
        }
        url.setQuery(QString::fromLatin1(encoded));
        paths << cachedDownload(url, this, progress);
      }
      std::vector<std::string> files;
      for (const auto &path : paths)
        files.push_back(path.toStdString());
      progress.setLabelText(tr("Clipping and projecting DEM…"));
      auto grid = Topography::load(files, r, [&](double value) {
        progress.setValue(int(value * 100));
        QApplication::processEvents();
        return !progress.wasCanceled();
      });
      if (progress.wasCanceled())
        return;
      survey.topography().grid = std::move(grid);
      visible->setChecked(true);
      resetRange();
      refresh();
      changed();
    } catch (const std::exception &e) {
      error(e);
    }
  });
  connect(remove, &QPushButton::clicked, this, [=, &survey] {
    survey.topography().grid.reset();
    refresh();
    changed();
  });
  connect(sample, &QPushButton::clicked, this, [=, &survey] {
    try {
      sampleStations(this, survey, changed);
    } catch (const std::exception &e) {
      error(e);
    }
  });
  auto exportGrid = [=, &survey](bool isTiff) {
    try {
      if (!survey.topography().grid)
        return;
      auto r = request();
      const auto file = QFileDialog::getSaveFileName(this, tr("Export topography"),
                                                     isTiff ? "topography.tif" : "topography.xyz",
                                                     isTiff ? tr("GeoTIFF (*.tif)") : tr("XYZ (*.xyz)"));
      if (file.isEmpty())
        return;
      QProgressDialog progress(tr("Preparing export…"), tr("Cancel"), 0, 100, this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0);
      progress.setAutoClose(false);
      auto grid = Topography::reproject(*survey.topography().grid, r, [&](double value) {
        progress.setValue(int(value * 100));
        QApplication::processEvents();
        return !progress.wasCanceled();
      });
      if (progress.wasCanceled())
        return;
      QTemporaryDir temporary;
      if (!temporary.isValid())
        throw std::runtime_error("Cannot create temporary export.");
      const auto path = temporary.filePath(isTiff ? "dem.tif" : "dem.xyz");
      const auto c = survey.coordinates().utm ? survey.coordinates()
                                              : SurveyCoordinates::suggested(survey.geographic_locations());
      if (isTiff)
        Topography::writeGeoTiff(path.toStdString(), *grid, c);
      else {
        Topography::writeXYZ(path.toStdString(), *grid, c, sign->currentIndex() == 1);
      }
      if (progress.wasCanceled())
        return;
      publish(path, file);
      status->setText(tr("Topography exported."));
    } catch (const std::exception &e) {
      error(e);
    }
  };
  connect(tif, &QPushButton::clicked, this, [=] { exportGrid(true); });
  connect(xyz, &QPushButton::clicked, this, [=] { exportGrid(false); });
  if (survey.topography().grid) {
    const auto g = survey.topography().grid;
    if (g->epsg == initial.epsg) {
      extent[0]->setValue(g->west);
      extent[1]->setValue(g->north - g->height * g->spacing);
      extent[2]->setValue(g->west + g->width * g->spacing);
      extent[3]->setValue(g->north);
      spacing->setValue(g->spacing);
    }
  }
  refresh();
}
