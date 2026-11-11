#include "ProfilePlots.h"
#include "ProfileWindow.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QTabWidget>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <fstream>
#include <iostream>
#include <numeric>

namespace {
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
bool near(double a, double b) { return std::abs(a - b) < 1e-8; }
template <class T> T *control(QWidget &w, const char *name) {
  auto *p = w.findChild<T *>(name);
  check(p, "Missing profile control");
  return p;
}
void calculations() {
  using namespace ProfileAnalysis;
  std::vector<Point> line{{{0, 0}}, {{10, 0}}, {{10, 10}}};
  std::vector<Station> sites{
      {"A", {{2, 1}}}, {"B", {{11, 5}}}, {"C", {{10, 11}}}, {"Outside", {{5, 3}}}, {"Invalid", {{NAN, 0}}}};
  auto chosen = select(sites, line, 1.);
  check(chosen.size() == 3 && chosen[0].name == "A" && near(chosen[0].along, 2) && near(chosen[0].across, 1),
        "Straight projection wrong");
  check(near(chosen[1].along, 15) && near(chosen[1].across, -1) && near(chosen[2].along, 20),
        "Bend or endpoint projection wrong");
  std::reverse(line.begin(), line.end());
  auto reverse = select(sites, line, 1.);
  check(reverse.size() == 3 && reverse.back().name == "A" && near(reverse.back().along, 18) &&
            near(reverse.back().across, -1),
        "Reversing profile changes membership or offsets");
  check(select(sites, {{{0, 0}}, {{0, 0}}}, 10).empty(), "Degenerate line selected stations");
  check(select({{"Cross", {{0, 0}}}}, {{{-1, -1}}, {{1, 1}}, {{-1, 1}}, {{1, -1}}}, 0).size() == 1,
        "Crossing duplicated station");
  check(near(medianSpacing({{"A", {{0, 0}}}, {"B", {{0, 0}}}, {"C", {{3, 0}}}, {"D", {{10, 0}}}}), 3),
        "Duplicate positions bias spacing");
  check(medianSpacing({{"A", {{0, 0}}}}) == 0, "Single station spacing undefined");
  const auto bounds = cells({0, .1, 2}, .3);
  check(near(bounds[0][1], .05) && bounds[1][1] < bounds[2][0], "Large gaps filled");
  const auto coincident = cells({0, 0, 0, 2}, 1.);
  check(coincident[0] == coincident[1] && coincident[1] == coincident[2] && near(coincident[1][1], .5),
        "Coincident stations lose their cells");
  auto roseResult = rose({{"A", 179}, {"B", 1}, {"C", NAN}}, 10, true);
  check(roseResult.samples == 2 && roseResult.stations == 2 &&
            (roseResult.mean < 1e-8 || roseResult.mean > 180 - 1e-8) && roseResult.concentration > .99,
        "Axial wrap or invalid direction handling wrong");
  roseResult = rose({{"A", 0}, {"A", 0}, {"A", 0}, {"B", 90}}, 10, true);
  check(near(roseResult.percent[0], 50) && near(roseResult.percent[9], 50) && std::isnan(roseResult.mean),
        "Equal station weighting wrong");
  check(near(rose({{"A", 0}, {"A", 0}, {"A", 0}, {"B", 90}}, 10, false).percent[0], 75),
        "Sample weighting wrong");
}
void click(QCustomPlot *map, double x, double y) {
  QPointF point(map->xAxis->coordToPixel(x), map->yAxis->coordToPixel(y));
  QMouseEvent down(QEvent::MouseButtonPress, point, point, Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
  QMouseEvent up(QEvent::MouseButtonRelease, point, point, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(map, &down);
  QApplication::sendEvent(map, &up);
  QApplication::processEvents();
}
void save(QWidget &window, const char *button, const QString &path) {
  bool handled = false;
  QTimer timer;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    for (auto *w : QApplication::topLevelWidgets())
      if (auto *dialog = qobject_cast<QFileDialog *>(w))
        if (dialog->isVisible()) {
          dialog->setDirectory(QFileInfo(path).absolutePath());
          control<QLineEdit>(*dialog, "fileNameEdit")->setText(path);
          handled = true;
          timer.stop();
          QMetaObject::invokeMethod(dialog, "accept", Qt::QueuedConnection);
          return;
        }
  });
  timer.start(10);
  control<QPushButton>(window, button)->click();
  check(handled, "Export dialog not opened");
  check(QFileInfo(path).size() > 100, "Export empty");
}
void dialog(const QString &directory) {
  const auto input = directory + "/synthetic.gofem";
  {
    std::ofstream out(input.toStdString());
    for (const char *name : {"A", "B", "C", "D"})
      for (double frequency : {1., .1, .001}) {
        const std::vector<std::pair<const char *, double>> values{
            {"RealZxx", .001}, {"ImagZxx", .001}, {"RealZxy", .01},  {"ImagZxy", .02},
            {"RealZyx", -.01}, {"ImagZyx", -.01}, {"RealZyy", .001}, {"ImagZyy", .002},
            {"RealTzx", .1},   {"ImagTzx", .05},  {"RealTzy", .2},   {"ImagTzy", .1}};
        for (const auto &v : values)
          out << v.first << ' ' << frequency << " Plane_wave " << name << ' ' << v.second << " .0001\n";
      }
  }
  auto survey = std::make_shared<MTSurveyData>();
  survey->load_from_gofem(input.toStdString());
  survey->get_station_data("A").set_position({40., -111., 0.});
  survey->get_station_data("B").set_position({40.01, -111., 0.});
  survey->get_station_data("C").set_position({40.02, -111., 0.});
  survey->get_station_data("D").set_position({40.01, -110.9, 0.});
  std::shared_ptr<MTSurveyData> extracted;
  ProfileWindow window(nullptr, [&](auto data) { extracted = std::move(data); });
  const std::map<std::string, MTSurveyData> responses{{"calculated", *survey}};
  window.setData(survey, responses);
  window.show();
  QApplication::processEvents();
  auto *map = control<QCustomPlot>(window, "profileMap");
  check(map->graph(0)->dataCount() == 4, "Survey map incomplete");
  auto *distance = control<QDoubleSpinBox>(window, "profileDistance");
  check(distance->value() > 1. && distance->value() < 1.2, "Default corridor is not median nearest spacing");
  distance->setValue(.5);
  const auto projection = SurveyCoordinates::suggested(survey->geographic_locations());
  const auto points = projection.transform({{39.999, -111., 0.}, {40.021, -111., 0.}});
  for (const auto &p : points)
    click(map, p[1] / 1000., p[0] / 1000.);
  auto *table = control<QTableWidget>(window, "profileStations");
  check(table->rowCount() == 3, "Drawing corridor selects wrong stations");
  auto *xy = control<QCustomPlot>(window, "profileSection0"),
       *yx = control<QCustomPlot>(window, "profileSection1"),
       *pt = control<QCustomPlot>(window, "profileTensorPlot");
  check(xy->property("profileSampleCount").toUInt() == 9 && yx->property("profileSampleCount").toUInt() == 9,
        "Sections lost samples");
  const auto yRange = xy->yAxis->range();
  check(xy->yAxis->rangeReversed(), "Period direction wrong");
  xy->xAxis->setRange(.2, 1.5);
  check(near(yx->xAxis->range().lower, .2) && near(pt->xAxis->range().upper, 1.5),
        "Section axes are not linked");
  survey->get_station_data("B").set_data_mask(RealZxy, .1, false);
  window.setData(survey, responses);
  check(table->rowCount() == 3 && xy->property("profileSampleCount").toUInt() == 8 &&
            yx->property("profileSampleCount").toUInt() == 9,
        "Refresh loses profile or propagates masks incorrectly");
  check(near(xy->yAxis->range().upper, yRange.upper), "Masked sample changes section bounds");
  control<QComboBox>(window, "profileDataset")->setCurrentIndex(1);
  check(xy->property("profileSampleCount").toUInt() == 9, "Calculated response not plotted");
  control<QComboBox>(window, "profileVertical")->setCurrentIndex(1);
  check(!xy->yAxis->rangeReversed() && near(xy->yAxis->range().lower, 1. / yRange.upper),
        "Frequency view wrong");
  control<QComboBox>(window, "profilePeriodMin")->setCurrentIndex(1);
  check(xy->property("profileSampleCount").toUInt() == 6, "Period selection ignored");
  control<QCheckBox>(window, "profileDirection0")->setChecked(false);
  check(control<QCustomPlot>(window, "profileRosePlot")->property("profileDirectionCount").toUInt() == 6,
        "Tipper rose missing");
  control<QPushButton>(window, "profileReverse")->click();
  check(table->item(0, 0)->text() == "C", "Reverse order wrong");
  control<QPushButton>(window, "profileExtract")->click();
  check(extracted && extracted->get_stations_names().size() == 3 && survey->get_stations_names().size() == 4,
        "Extraction modifies original survey");
  check(!extracted->get_station_data("B").impedance_mask()[1][1], "Extraction loses masks");
  save(window, "profileCsv", directory + "/profile.csv");
  save(window, "profileSavePdf", directory + "/profile.pdf");
  QFile csv(directory + "/profile.csv");
  csv.open(QIODevice::ReadOnly);
  const auto bytes = csv.readAll();
  check(bytes.contains("origin_easting_m=") && !bytes.contains(directory.toUtf8()),
        "CSV missing origin or leaks path");
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if (!artifacts.isEmpty()) {
    window.grab().save(artifacts + "/profile-sections.png");
    control<QTabWidget>(window, "profileTabs")->setCurrentIndex(1);
    QApplication::processEvents();
    window.grab().save(artifacts + "/profile-tensors.png");
    control<QTabWidget>(window, "profileTabs")->setCurrentIndex(2);
    QApplication::processEvents();
    window.grab().save(artifacts + "/profile-roses.png");
    QFile::copy(directory + "/profile.pdf", artifacts + "/profile.pdf");
  }
  control<QPushButton>(window, "profileClear")->click();
  check(table->rowCount() == 0 && !control<QPushButton>(window, "profileExtract")->isEnabled(),
        "Clear leaves selection");
}
} // namespace
int main(int argc, char **argv) {
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  try {
    calculations();
    QTemporaryDir directory;
    check(directory.isValid(), "Temporary directory failed");
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
    dialog(directory.path());
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cout << "Profile geometry, sections, masks, roses and export checks passed.\n";
}
