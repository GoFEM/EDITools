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
void strikeEstimation() {
  using namespace ProfileAnalysis;
  const auto crossing = estimateStrike({{"A", 89.}, {"B", 1.}, {"C", NAN}}, true);
  check(crossing.samples == 2 && crossing.stations == 2 &&
            (crossing.degrees < 1e-8 || crossing.degrees > 90. - 1e-8) && crossing.concentration > .99,
        "Strike wrap or invalid samples mishandled");
  const auto perpendicular = estimateStrike({{"A", 10.}, {"B", 100.}}, true);
  check(near(perpendicular.degrees, 10.) && near(perpendicular.concentration, 1.),
        "Perpendicular PT axes must infer the same strike family");
  check(std::isnan(estimateStrike({{"A", 0.}, {"B", 45.}}, true).degrees),
        "Undefined strike reported as a direction");
  const std::vector<Direction> uneven{{"A", 10.}, {"A", 10.}, {"A", 10.}, {"B", 40.}};
  check(near(estimateStrike(uneven, true).degrees, 25.) && estimateStrike(uneven, false).degrees < 25.,
        "Strike station weighting ignored");
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
  auto *rotate = control<QCheckBox>(window, "profileRotateStrike");
  auto *strikeMode = control<QComboBox>(window, "profileStrikeMode");
  auto *strikeAngle = control<QDoubleSpinBox>(window, "profileStrikeAngle");
  check(!rotate->isChecked() && strikeMode->currentIndex() == 0 && strikeAngle->isReadOnly() &&
            (std::abs(strikeAngle->value()) < .01 || std::abs(strikeAngle->value() - 90.) < .01),
        "Automatic observed strike default wrong");
  const auto mapPosition = map->graph(1)->data()->at(1)->value;
  rotate->setChecked(true);
  control<QPushButton>(window, "profileSwapStrike")->click();
  check(xy->property("profileSampleCount").toUInt() == 9, "Automatic axis swap loses data");
  strikeMode->setCurrentIndex(1);
  strikeAngle->setValue(45.);
  survey->get_station_data("B").set_data_mask(RealZxx, .1, false);
  window.setData(survey, responses);
  check(xy->property("profileSampleCount").toUInt() == 8 &&
            yx->property("profileSampleCount").toUInt() == 8 &&
            near(map->graph(1)->data()->at(1)->value, mapPosition) && table->rowCount() == 3,
        "Rotation ignores required masks or moves the corridor");
  auto *ptDirections = control<QCheckBox>(window, "profileDirection0");
  ptDirections->setChecked(false);
  check(control<QCustomPlot>(window, "profileRosePlot")->property("profileDirectionCount").toUInt() == 9,
        "Impedance rotation masks leaked into tipper directions");
  ptDirections->setChecked(true);
  strikeAngle->setValue(90.);
  check(xy->property("profileSampleCount").toUInt() == 9 && yx->property("profileSampleCount").toUInt() == 9,
        "Quarter-turn masks unrelated components");
  rotate->setChecked(false);
  check(!survey->get_station_data("B").impedance_mask()[0][1], "Display rotation changed source masks");
  survey->get_station_data("B").set_data_mask(RealZxx, .1, true);
  strikeMode->setCurrentIndex(0);
  control<QPushButton>(window, "profileSwapStrike")->click();
  window.setData(survey, responses);
  rotate->setChecked(true);
  auto *minShape = control<QDoubleSpinBox>(window, "profileMinEllipticity");
  const double oldShape = minShape->value();
  minShape->setValue(1.);
  check(xy->property("profileSampleCount").toUInt() == 0 &&
            !control<QPushButton>(window, "profileSavePdf")->isEnabled(),
        "Unavailable automatic strike silently uses unrotated data");
  minShape->setValue(oldShape);
  rotate->setChecked(false);
  check(xy->property("profileSampleCount").toUInt() == 9, "Rotation recovery loses original data");
  auto labels = qSharedPointerDynamicCast<QCPAxisTickerText>(xy->xAxis2->ticker());
  check(labels && labels->ticks().size() == 3 && labels->ticks().values().contains("B") &&
            xy->xAxis2->tickLabelFont().pointSize() == 7,
        "Station annotations missing");
  const double bOffset = table->item(1, 1)->text().toDouble();
  click(xy, bOffset, 10.);
  check(table->currentRow() == 1 && table->selectionModel()->selectedRows().size() == 1 &&
            map->graph(3)->dataCount() == 1,
        "Clicking a column does not select its map/table station");
  const auto bPosition = projection.transform({survey->get_station_data("B").position()})[0];
  check(near(map->graph(3)->data()->at(0)->key, bPosition[1] / 1000.) &&
            near(map->graph(3)->data()->at(0)->value, bPosition[0] / 1000.),
        "Highlighted map station is wrong");
  const auto dragRange = xy->xAxis->range();
  const QPointF dragStart(xy->xAxis->coordToPixel(bOffset), xy->yAxis->coordToPixel(10.));
  const QPointF dragEnd = dragStart + QPointF(25, 0);
  QMouseEvent down(QEvent::MouseButtonPress, dragStart, dragStart, Qt::LeftButton, Qt::LeftButton,
                   Qt::NoModifier);
  QMouseEvent move(QEvent::MouseMove, dragEnd, dragEnd, Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
  QMouseEvent up(QEvent::MouseButtonRelease, dragEnd, dragEnd, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(xy, &down);
  QApplication::sendEvent(xy, &move);
  QApplication::sendEvent(xy, &up);
  check(table->currentRow() == 1, "Panning changes station selection");
  xy->xAxis->setRange(dragRange);
  auto *tabs = control<QTabWidget>(window, "profileTabs");
  check(tabs->count() == 2, "Tensor and strike views are not combined");
  tabs->setCurrentIndex(1);
  QApplication::processEvents();
  check(pt->isVisible() && control<QCustomPlot>(window, "profileRosePlot")->isVisible(),
        "Combined tab hides a plot");
  auto *tensorPalette = control<QComboBox>(window, "profileTensorPalette");
  const auto beforeColor = pt->findChild<QCPColorScale *>()->gradient().color(.5, QCPRange(0, 1));
  tensorPalette->setCurrentIndex((tensorPalette->currentIndex() + 1) % tensorPalette->count());
  check(control<QComboBox>(window, "profilePalette")->currentIndex() == tensorPalette->currentIndex() &&
            beforeColor != pt->findChild<QCPColorScale *>()->gradient().color(.5, QCPRange(0, 1)),
        "Combined-tab colormap does not update plots");
  click(pt, table->item(2, 1)->text().toDouble(), 10.);
  check(table->currentRow() == 2 && table->item(2, 0)->text() == "C", "Tensor column selection failed");
  table->setCurrentCell(0, 0);
  const auto aPosition = projection.transform({survey->get_station_data("A").position()})[0];
  check(near(map->graph(3)->data()->at(0)->value, aPosition[0] / 1000.),
        "Table selection does not update map");
  tabs->setCurrentIndex(0);
  QApplication::processEvents();
  const auto yRange = xy->yAxis->range();
  check(xy->yAxis->rangeReversed(), "Period direction wrong");
  xy->xAxis->setRange(.2, 1.5);
  check(near(xy->xAxis2->range().lower, .2) && near(yx->xAxis->range().lower, .2) &&
            near(pt->xAxis->range().upper, 1.5),
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
  check(table->currentRow() == 2 && table->item(2, 0)->text() == "A", "Reverse changes highlighted station");
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

    QFile::copy(directory + "/profile.pdf", artifacts + "/profile.pdf");
  }
  control<QPushButton>(window, "profileClear")->click();
  check(table->rowCount() == 0 && map->graph(3)->dataCount() == 0 &&
            !control<QPushButton>(window, "profileExtract")->isEnabled(),
        "Clear leaves selection");
}
} // namespace
int main(int argc, char **argv) {
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  try {
    calculations();
    strikeEstimation();
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
