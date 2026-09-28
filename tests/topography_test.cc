#include "include/Topography.h"
#include "include/MapBackground.h"
#include "include/TopographyLayer.h"
#include "include/SurveyReport.h"
#include "TopographyDialog.h"
#include "PeriodMapWindow.h"
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QProcess>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
void near(double a, double b, double tolerance = 1e-4) {
  check(std::isfinite(a) && std::abs(a - b) <= tolerance, "Unexpected DEM value");
}
template <class T> T *control(QWidget &w, const char *name) {
  auto *p = w.findChild<T *>(name);
  check(p, "Missing topography control");
  return p;
}
template <class F> void rejects(F f) {
  bool failed = false;
  try {
    f();
  } catch (const std::exception &) {
    failed = true;
  }
  check(failed, "Invalid DEM operation succeeded");
}
Topography::Grid fixture() {
  Topography::Grid g;
  g.width = g.height = 5;
  g.epsg = 32611;
  g.west = 500000;
  g.north = 4100000;
  g.spacing = 10;
  g.source = "Synthetic terrain";
  g.datum = "Synthetic vertical reference";
  g.notice = "Synthetic licence notice";
  for (int r = 0; r < 5; ++r)
    for (int c = 0; c < 5; ++c)
      g.values.push_back(100 + 10 * c + 20 * r);
  g.values[12] = NAN;
  return g;
}
MTSurveyData survey(const QString &dir, const Topography::Grid &g) {
  std::vector<double> x{g.west + 15}, y{g.north - 15};
  Topography::transform(g.epsg, 4326, x, y);
  const auto file = dir + "/A.edi";
  std::ofstream f(file.toStdString());
  f.precision(17);
  f << ">HEAD\nDATAID=A\nLAT=" << y[0] << "\nLONG=" << x[0]
    << "\nELEV=42\nEMPTY=1e32\n>=MTSECT\n>FREQ // 1\n1\n";
  for (auto name : {"ZXXR", "ZXXI", "ZXYR", "ZXYI", "ZYXR", "ZYXI", "ZYYR", "ZYYI"})
    f << '>' << name << " // 1\n1\n";
  for (auto name : {"ZXX.VAR", "ZXY.VAR", "ZYX.VAR", "ZYY.VAR"})
    f << '>' << name << " // 1\n0.01\n";
  f << ">END\n";
  f.close();
  MTSurveyData s;
  s.load_from_edi({file.toStdString()});
  auto c = SurveyCoordinates::suggested(s.geographic_locations());
  c.centered = true;
  c.origin_easting = 500000;
  c.origin_northing = 4100000;
  s.set_coordinates(c);
  return s;
}
void chooseFile(QWidget &dialog, const char *button, const QString &path) {
  bool handled = false;
  QTimer timer;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    auto *f = qobject_cast<QFileDialog *>(QApplication::activeModalWidget());
    if (!f)
      return;
    f->setDirectory(QFileInfo(path).absolutePath());
    auto *edit = f->findChild<QLineEdit *>("fileNameEdit");
    if (!edit)
      return;
    edit->setText(QFileInfo(path).fileName());
    handled = true;
    timer.stop();
    QMetaObject::invokeMethod(f, "accept", Qt::QueuedConnection);
  });
  timer.start(10);
  control<QPushButton>(dialog, button)->click();
  check(handled, "DEM file dialog missing");
}
void tests(const QString &dir) {
  QTimer messages;
  QObject::connect(&messages, &QTimer::timeout, [] {
    auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget());
    if (box) {
      std::cerr << "Unexpected dialog: " << box->text().toStdString() << '\n';
      box->accept();
    }
  });
  messages.start(20);
  auto g = fixture();
  g.validate();
  near(g.sample(500005, 4099995), 100);
  near(g.sample(500010, 4099990), 115);
  check(std::isnan(g.sample(500025, 4099975)), "NoData was filled");
  check(std::isnan(g.sample(499999, 4099995)), "DEM extrapolated");
  Topography::Request r;
  r.epsg = g.epsg;
  r.west = g.west;
  r.north = g.north;
  r.east = g.west + 50;
  r.south = g.north - 50;
  r.spacing = 10;
  r.source = g.source;
  r.datum = g.datum;
  r.notice = g.notice;
  const auto tif = dir + "/terrain.tif";
  Topography::writeGeoTiff(tif.toStdString(), g);
  auto loaded = Topography::load({tif.toStdString()}, r);
  near(loaded->sample(500010, 4099990), 115);
  check(std::isnan(loaded->values[12]), "GeoTIFF lost NoData");
  auto *ds = static_cast<GDALDataset *>(GDALOpen(tif.toUtf8().constData(), GA_ReadOnly));
  check(ds, "GeoTIFF cannot reopen");
  check(std::string(ds->GetMetadataItem("VERTICAL_DATUM")) == g.datum, "Missing datum");
  check(std::string(ds->GetMetadataItem("DEM_LICENCE_NOTICE")) == g.notice, "Missing licence");
  double gt[6];
  ds->GetGeoTransform(gt);
  near(gt[0], g.west);
  near(gt[3], g.north);
  GDALClose(ds);
  auto feet = r;
  feet.unitScale = .3048;
  auto metric = Topography::load({tif.toStdString()}, feet);
  near(metric->values[0], 30.48);
  // Import a geographic raster into UTM, with explicit scale/offset and units.
  const auto geographic = dir + "/geographic.tif";
  auto *geo = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(geographic.toUtf8().constData(), 10,
                                                                       10, 1, GDT_Float32, nullptr);
  check(geo, "Cannot create geographic fixture");
  OGRSpatialReference crs;
  crs.importFromEPSG(4326);
  crs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  geo->SetSpatialRef(&crs);
  double affine[]{-117.001, .0002, 0., 37.001, 0., -.0002};
  geo->SetGeoTransform(affine);
  geo->GetRasterBand(1)->Fill(100);
  geo->GetRasterBand(1)->SetScale(2);
  geo->GetRasterBand(1)->SetOffset(10);
  geo->GetRasterBand(1)->SetUnitType("ft");
  GDALClose(geo);
  const auto region = Topography::surveyArea({{{37., -117., 0.}}}, SurveyCoordinates{}, 20, 10);
  auto projected = Topography::load({geographic.toStdString()}, region);
  near(projected->range()[0], 210 * .3048);
  near(projected->range()[1], 210 * .3048);
  auto clip = r;
  clip.west += 10;
  clip.east -= 10;
  auto cropped = Topography::reproject(g, clip);
  check(cropped->width == 3 && cropped->height == 5, "Clipping failed");
  near(cropped->values[0], 110);
  check(cropped->notice == g.notice, "Clipping lost licence");
  auto huge = r;
  huge.spacing = .001;
  rejects([&] { huge.validate(); });
  rejects([&] { Topography::load({tif.toStdString()}, r, [](double) { return false; }); });
  auto empty = r;
  empty.west += 1000;
  empty.east += 1000;
  rejects([&] { Topography::load({tif.toStdString()}, empty); });
  auto s = survey(dir, g);
  s.topography().grid = loaded;
  s.topography().style.lower = 100;
  s.topography().style.upper = 220;
  const auto xyz = dir + "/terrain.xyz";
  Topography::writeXYZ(xyz.toStdString(), g, s.coordinates(), true);
  std::ifstream f(xyz.toStdString());
  std::string line;
  int rows = 0;
  while (std::getline(f, line))
    if (!line.empty() && line[0] != '#') {
      std::istringstream row(line);
      double n, e, z;
      row >> n >> e >> z;
      if (!rows) {
        near(n, -5);
        near(e, 5);
        near(z, -100);
      }
      ++rows;
    }
  check(rows == 24, "XYZ exported missing cells");
  std::stringstream saved;
  {
    boost::archive::binary_oarchive a(saved);
    a << s;
  }
  MTSurveyData restored;
  {
    boost::archive::binary_iarchive a(saved);
    a >> restored;
  }
  check(restored.topography().grid && restored.topography().grid->notice == g.notice,
        "Project lost topography");
  near(restored.topography().grid->sample(500010, 4099990), 115);
  QFile raster(tif);
  check(raster.open(QIODevice::ReadOnly), "Cannot inspect GeoTIFF metadata");
  check(!raster.readAll().contains(dir.toUtf8()), "GeoTIFF leaked a local path");
  QCustomPlot plot;
  MapBackground background(&plot);
  background.setTopography(s.topography());
  background.setCoordinates(s.coordinates());
  plot.xAxis->setRange(0, 50);
  plot.yAxis->setRange(-50, 0);
  plot.resize(500, 400);
  plot.replot();
  const auto terrain = plot.toPixmap(500, 400);
  background.setTopography({});
  plot.replot();
  check(terrain.toImage() != plot.toPixmap(500, 400).toImage(), "DEM missing from rendered map");
  auto data = std::make_shared<MTSurveyData>(s);
  PeriodMapWindow maps(nullptr);
  maps.setData(data, {});
  check(maps.findChild<TopographyLayer *>("topographyLayer"), "Period map lacks DEM layer");
  unsigned changes = 0;
  TopographyDialog dialog(nullptr, s, [&] { ++changes; });
  dialog.show();
  QApplication::processEvents();
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if (!artifacts.isEmpty()) {
    dialog.grab().save(artifacts + "/topography-dialog.png");
    terrain.save(artifacts + "/topography-map.png");
  }
  control<QCheckBox>(dialog, "demVisible")->setChecked(false);
  check(!s.topography().style.visible, "Visibility did not reach survey");
  control<QCheckBox>(dialog, "demVisible")->setChecked(true);
  control<QDoubleSpinBox>(dialog, "demOpacity")->setValue(40);
  near(s.topography().style.opacity, .4);
  const auto masks = s.get_station_data("A").impedance_mask();
  const auto oldHeight = s.get_station_data("A").position()[2];
  QTimer::singleShot(0, [&] {
    auto *p = QApplication::activeModalWidget();
    check(p, "Station preview missing");
    auto *t = control<QTableWidget>(*p, "demStationSamples");
    check(t->rowCount() == 1, "Station preview rows wrong");
    near(t->item(0, 2)->text().toDouble(), 130, .02);
    auto *b = p->findChild<QDialogButtonBox *>();
    check(!b->button(QDialogButtonBox::Apply)->isEnabled(), "Datum confirmation absent");
    qobject_cast<QDialog *>(p)->reject();
  });
  control<QPushButton>(dialog, "demSampleStations")->click();
  near(s.get_station_data("A").position()[2], oldHeight);
  QTimer::singleShot(0, [&] {
    auto *p = QApplication::activeModalWidget();
    control<QCheckBox>(*p, "demDatumChecked")->setChecked(true);
    p->findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Apply)->click();
  });
  control<QPushButton>(dialog, "demSampleStations")->click();
  near(s.get_station_data("A").position()[2], 130, .02);
  check(s.get_station_data("A").impedance_mask() == masks, "Sampling heights changed masks");
  chooseFile(dialog, "demExportXYZ", dir + "/gui.xyz");
  check(QFileInfo(dir + "/gui.xyz").size() > 0, "GUI XYZ export failed");
  chooseFile(dialog, "demExportTiff", dir + "/gui.tif");
  auto *exported = static_cast<GDALDataset *>(GDALOpen((dir + "/gui.tif").toUtf8().constData(), GA_ReadOnly));
  check(exported, "GUI GeoTIFF export failed");
  check(exported->GetMetadataItem("SURVEY_ORIGIN_EASTING_M"), "Missing survey origin");
  near(std::stod(exported->GetMetadataItem("SURVEY_ORIGIN_EASTING_M")), 500000);
  near(std::stod(exported->GetMetadataItem("SURVEY_ORIGIN_NORTHING_M")), 4100000);
  GDALClose(exported);
  control<QPushButton>(dialog, "demRemove")->click();
  check(!s.topography().grid, "Remove did not clear terrain");
  chooseFile(dialog, "demLoad", tif);
  check(bool(s.topography().grid), "GUI GeoTIFF import failed");
  SurveyReport::Options options;
  options.stationPages = false;
  options.topography = true;
  const auto report = dir + "/report.pdf";
  check(SurveyReport::write_pdf(report, s, nullptr, options), "Report failed");
  QProcess text;
  text.start("pdftotext", {report, "-"});
  if (text.waitForFinished() && text.exitCode() == 0) {
    const auto body = text.readAllStandardOutput();
    check(body.contains("Synthetic terrain"), "Report omitted DEM credit");
    check(body.contains("Synthetic licence notice"), "Report omitted DEM licence");
    check(!body.contains(dir.toUtf8()), "Report leaked local path");
  }
  check(changes > 0, "Topography changes not notified");
}
void online(const QString &dir, int provider) {
  auto g = fixture();
  auto s = survey(dir, g);
  TopographyDialog dialog(nullptr, s, [] {});
  dialog.show();
  control<QComboBox>(dialog, "demSource")->setCurrentIndex(provider);
  control<QDoubleSpinBox>(dialog, "demMargin")->setValue(.2);
  control<QDoubleSpinBox>(dialog, "demSpacing")->setValue(100);
  control<QPushButton>(dialog, "demSurveyExtent")->click();
  if (provider == 1)
    control<QCheckBox>(dialog, "demAcceptLicence")->setChecked(true);
  QString error;
  QTimer watchdog;
  QObject::connect(&watchdog, &QTimer::timeout, [&] {
    auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget());
    if (box) {
      error = box->text();
      box->accept();
    }
  });
  watchdog.start(50);
  control<QPushButton>(dialog, "demLoad")->click();
  if (!error.isEmpty())
    throw std::runtime_error(error.toStdString());
  check(bool(s.topography().grid), "Online DEM absent");
  auto range = s.topography().grid->range();
  std::cout << "Provider " << provider << ": " << s.topography().grid->width << " x "
            << s.topography().grid->height << ", heights " << range[0] << " to " << range[1] << " m\n";
}
} // namespace
int main(int argc, char **argv) {
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("EDIToolsTests");
  QCoreApplication::setApplicationName("Topography");
  QTemporaryDir dir;
  try {
    check(dir.isValid(), "No test directory");
    if (argc > 1 && std::string(argv[1]) == "--online") {
      if (argc > 2)
        online(dir.path(), std::stoi(argv[2]));
      else {
        online(dir.path(), 2);
        online(dir.path(), 1);
      }
    } else
      tests(dir.path());
    std::cout << "Topography checks passed.\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
