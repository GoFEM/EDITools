#include "ui/dialogs/SurveyCoordinatesDialog.h"
#include "core/MTSurveyData.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QTemporaryDir>
#include <cmath>
#include <fstream>
#include <iostream>

namespace {
void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
template <class T> T *control(QWidget &parent, const char *name) {
  auto *value = parent.findChild<T *>(name);
  check(value != nullptr, "Missing UTM control");
  return value;
}
void calculations() {
  std::vector<std::array<double, 3>> points{{40., -119., 0.}, {42., -103., 0.}};
  const auto suggested = SurveyCoordinates::suggested(points);
  check(suggested.zone == 12 && suggested.north, "Default is not geographic extent centre");
  for (int i = 0; i < 20; ++i)
    points.push_back(points.front());
  check(SurveyCoordinates::suggested(points).zone == 12, "Station density biases default zone");
  points = {{10., 179., 0.}, {12., -179., 0.}};
  const auto bounds = UTM::extent(points);
  check(std::abs(bounds.east - bounds.west - 2.) < 1e-10 &&
            std::abs(bounds.unwrap(-179.) - bounds.unwrap(179.) - 2.) < 1e-10,
        "Date-line extent discontinuity");
  const auto dateZone = SurveyCoordinates::suggested(points);
  check(dateZone.zone == 60 || dateZone.zone == 1, "Date-line default zone wrong");
  auto analysis = UTM::analyze(points, dateZone);
  check(analysis.valid() && analysis.stationZones.size() == 2 && analysis.outsideZone == 1 &&
            analysis.extentMaxError < .01,
        "Date-line analysis wrong");
  const auto overlays = UTM::regions(bounds);
  check(overlays.size() == 2 && overlays[0].zone == 60 && overlays[1].zone == 1, "Date-line overlays wrong");
  check(UTM::zoneAt(60., 4.) == 32 && UTM::zoneAt(75., 8.) == 31 && UTM::zoneAt(75., 10.) == 33 &&
            UTM::zoneAt(75., 22.) == 35 && UTM::zoneAt(75., 34.) == 37,
        "Special zone selection wrong");
  for (const auto &r : UTM::regions({0., 42., 55., 84.})) {
    for (double x : {r.west + 1e-6, r.east - 1e-6})
      for (double y : {r.south + 1e-6, r.top - 1e-6})
        check(UTM::zoneAt(y, x) == r.zone, "Overlay crosses nominal zone boundary");
  }
  auto c = SurveyCoordinates::calculate({{0., 9., 0.}}, 32, true, false);
  analysis = UTM::analyze({{0., 9., 0.}}, c);
  check(analysis.valid() && std::abs(analysis.stations[0].scale - .9996) < 1e-8 &&
            std::abs(analysis.stations[0].convergence) < 1e-8,
        "Central-meridian scale/convergence wrong");
  // Compare factor scale with a small projected segment and exact WGS84 equatorial arc.
  const double delta = 1e-4;
  auto segment = c.transform({{0., 12. - delta, 0.}, {0., 12. + delta, 0.}});
  const double measured = std::hypot(segment[1][0] - segment[0][0], segment[1][1] - segment[0][1]) /
                          (6378137. * 2. * delta * std::acos(-1.) / 180.);
  analysis = UTM::analyze({{0., 12., 0.}}, c);
  check(std::abs(measured - analysis.stations[0].scale) < 1e-8, "Scale does not match projected distance");
  points = {{-1., 8., 0.}, {1., 14., 0.}};
  analysis = UTM::analyze(points, c);
  check(analysis.valid() && analysis.oppositeHemisphere == 1 && analysis.outsideZone == 1 &&
            analysis.stationZones.size() == 2 && analysis.extentMaxError >= analysis.stationMaxError,
        "Zone/hemisphere crossing analysis wrong");
  const auto centered = UTM::analyze(points, SurveyCoordinates::calculate(points, 32, true, true));
  check(centered.extentMaxError == analysis.extentMaxError, "Centering changed distortion");
  check(!UTM::analyze({{0., 99., 0.}}, c).valid(), "Singular projection accepted");
}
void dialog(const QString &directory) {
  const QString input = directory + "/survey.gofem";
  {
    std::ofstream file(input.toStdString());
    for (const auto *name : {"A", "B", "C"})
      file << "RealZxy 1 Plane_wave " << name << " .002 .0001\nImagZxy 1 Plane_wave " << name
           << " .002 .0001\n";
  }
  MTSurveyData survey;
  survey.load_from_gofem(input.toStdString());
  survey.get_station_data("A").set_position({40., -116., 100.});
  survey.get_station_data("B").set_position({41., -111., 200.});
  survey.get_station_data("C").set_position({42., -106., 300.});
  const auto original = survey.geographic_locations();
  SurveyCoordinatesDialog preview(survey);
  preview.show();
  QApplication::processEvents();
  auto *zone = control<QSpinBox>(preview, "utmZone");
  auto *table = control<QTableWidget>(preview, "utmStations");
  auto *candidates = control<QTableWidget>(preview, "utmCandidates");
  auto *map = control<QCustomPlot>(preview, "utmZoneMap");
  check(zone->value() == 12 && table->rowCount() == 3 && candidates->rowCount() >= 3 && map->itemCount() > 3,
        "Zone preview incomplete");
  check(control<QLabel>(preview, "utmCrossing")->text().contains("2 outside"), "Crossing summary wrong");
  zone->setValue(11);
  check(map->graph(0)->dataCount() == 1 && map->graph(1)->dataCount() == 2, "Zone map not updated");
  auto *ok = preview.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok);
  check(ok->isEnabled(), "High but finite distortion prevents zone selection");
  check(control<QLabel>(preview, "utmStatus")->text().startsWith("Above"), "High distortion not reported");
  control<QPushButton>(preview, "utmSuggested")->click();
  check(zone->value() == 12, "Survey-centre reset failed");
  control<QCheckBox>(preview, "utmCenter")->setChecked(true);
  check(std::abs(control<QLineEdit>(preview, "utmOriginEasting")->text().toDouble()) > 100000.,
        "Preview origin missing");
  check(!survey.coordinates().utm, "Preview changed survey before acceptance");
  const auto artifact = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if (!artifact.isEmpty())
    preview.grab().save(artifact + "/utm-dialog.png");
  ok->click();
  check(survey.coordinates().zone == 12 && survey.coordinates().centered &&
            survey.geographic_locations() == original,
        "Applying coordinates changed original positions");
  const auto applied = survey.coordinates();
  SurveyCoordinatesDialog cancelled(survey);
  check(control<QSpinBox>(cancelled, "utmZone")->value() == 12 &&
            control<QCheckBox>(cancelled, "utmCenter")->isChecked(),
        "Dialog lost saved settings");
  control<QSpinBox>(cancelled, "utmZone")->setValue(13);
  cancelled.reject();
  check(survey.coordinates().zone == applied.zone &&
            survey.coordinates().origin_easting == applied.origin_easting,
        "Cancel changed projection");
}
} // namespace
int main(int argc, char **argv) {
  QApplication app(argc, argv);
  try {
    calculations();
    QTemporaryDir dir;
    check(dir.isValid(), "Temporary directory failed");
    dialog(dir.path());
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cout << "UTM diagnostics and dialog checks passed.\n";
}
