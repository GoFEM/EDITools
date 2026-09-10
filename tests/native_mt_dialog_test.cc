#include "ExportNativeMTDialog.h"
#include "SurveyCoordinatesDialog.h"
#include "include/MTSurveyData.h"
#include <QApplication>
#include <QClipboard>
#include <QSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QTemporaryDir>
#include <QTimer>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdint>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace {
void check(bool ok, const char *message)
{
  if(!ok) throw std::runtime_error(message);
}
template<class T> T *widget(QWidget &parent, const char *name)
{
  auto *result = parent.findChild<T *>(name);
  check(result != nullptr, "Missing dialog control");
  return result;
}

void check_receiver_origin(const std::string &path, const SurveyCoordinates &coordinates)
{
  std::ifstream input(path);
  std::string line;
  const std::string prefix = "# UTM origin (easting northing, m): ";
  unsigned count = 0;
  bool has_zone = false;
  while(std::getline(input, line)) {
    has_zone |= line.find("# WGS84 / UTM " + std::to_string(coordinates.zone) +
                         (coordinates.north ? "N" : "S")) == 0;
    if(line.find(prefix) == 0) {
      std::istringstream values(line.substr(prefix.size()));
      double easting, northing;
      std::string extra;
      check(bool(values >> easting >> northing) && !(values >> extra) &&
            easting == coordinates.origin_easting && northing == coordinates.origin_northing,
            "Receiver comment lost UTM origin precision or axis order");
      ++count;
    }
  }
  check(has_zone && count == 1, "Receiver file needs a UTM zone and exactly one commented origin line");
}

void coordinate_checks(const QString &directory)
{
  auto projection = SurveyCoordinates::calculate({{{0.,9.,10.}}}, 32, true, false);
  auto p = projection.transform({{{0.,9.,10.}}})[0];
  check(std::abs(p[1]-500000.) < 1e-6 && std::abs(p[0]) < 1e-6 && p[2] == 10., "UTM central meridian");
  auto south = SurveyCoordinates::suggested({{{-33.8688,151.2093,20.}}});
  check(south.zone == 56 && !south.north, "Southern zone selection");
  p = south.transform({{{-33.8688,151.2093,20.}}})[0];
  check(std::abs(p[1]-334368.633648097) < 1e-6 && std::abs(p[0]-6250948.345385009) < 1e-6, "Southern UTM false northing");
  const std::vector<std::array<double,3>> points{{{55.,12.,100.}}, {{55.1,12.1,200.}}, {{55.02,12.02,300.}}};
  const auto centered = SurveyCoordinates::calculate(points, 32, true, true);
  const auto local = centered.transform(points);
  check(std::abs(local[0][0]+local[1][0]) < 1e-6 && std::abs(local[0][1]+local[1][1]) < 1e-6,
        "Center is not the projected bounding-box midpoint");
  bool rejected = false;
  try {projection.transform({{{90.,9.,0.}}});} catch(const std::exception &) {rejected = true;}
  check(rejected, "Polar latitude must be rejected");

  // Sidecars reflect emitted data, masks, period selection and disabled stations.
  const auto input = (directory + "/frequencies.gofem").toStdString();
  {
    std::ofstream data(input);
    data << "RealZxy 1.2345678901234567 Plane_wave S01 1 0.1\n"
         << "ImagZxy 1.2345678901234567 Plane_wave S01 2 0.1\n"
         << "RealZxy 2 Plane_wave S01 1 0.1\nImagZxy 2 Plane_wave S01 2 0.1\n"
         << "RealZxy 3 Plane_wave S01 1 0.1\nImagZxy 3 Plane_wave S01 2 0.1\n"
         << "RealZxy 4 Plane_wave S02 1 0.1\nImagZxy 4 Plane_wave S02 2 0.1\n";
  }
  MTSurveyData survey;
  survey.load_from_gofem(input);
  survey.get_station_data("S01").set_position(points[0]);
  survey.get_station_data("S02").set_position(points[1]);
  survey.get_station_data("S01").set_data_mask(RealZxy, 3., false);
  survey.set_active_flag("S02", false);
  const auto output = (directory + "/legacy").toStdString();
  const double precise = 1.2345678901234567;
  survey.write_gofem(output, {RealZxy}, {1./precise, .5, 1./3., .25});
  check_receiver_origin(output + ".recvs", SurveyCoordinates::suggested(survey.geographic_locations()));
  {
    std::ifstream coords(output + ".recvs");
    const auto position = NativeMT::read_receivers(coords).at("S01");
    check(std::abs(position[0]-6098907.825129169) < 1e-6 &&
          std::abs(position[1]-308124.367862458) < 1e-6,
          "GoFEM receiver file must use UTM even when the map is geographic");
  }
  {
    std::ifstream freq(output + ".freqs");
    double a, b;
    check(bool(freq >> a >> b) && a == precise && b == 2., "Frequency file must contain exact sorted values without a header");
    std::string extra;
    check(!(freq >> extra), "Frequency file contains data that were not exported");
  }
  survey.set_coordinates(SurveyCoordinates::calculate(survey.geographic_locations(), 32, true, true));
  const auto before = survey.get_stations_locations();
  check(survey.closest_station_name(before[0][0],before[0][1]) == "S01", "Map picking ignored projected coordinates");
  survey.write_gofem(output, {RealZxy}, {.5});
  check_receiver_origin(output + ".recvs", survey.coordinates());
  {
    std::ifstream coords(output + ".recvs");
    auto recvs = NativeMT::read_receivers(coords);
    check(recvs.at("S01") == before[0], "GoFEM export ignored survey coordinates");
    std::ifstream freq(output + ".freqs");
    double f; std::string extra;
    check(bool(freq >> f) && f == 2. && !(freq >> extra), "Period selection or frequency file replacement failed");
  }
  std::stringstream archive;
  {boost::archive::binary_oarchive out(archive); out << survey;}
  MTSurveyData restored;
  {boost::archive::binary_iarchive in(archive); in >> restored;}
  check(restored.coordinates().zone == 32 && restored.coordinates().centered &&
        restored.coordinates().origin_easting == survey.coordinates().origin_easting &&
        restored.coordinates().origin_northing == survey.coordinates().origin_northing,
        "Project save/load lost UTM settings");
  check(restored.get_stations_locations() == before && restored.geographic_locations() == survey.geographic_locations(),
        "Project save/load changed coordinates");
  restored.set_coordinates(SurveyCoordinates{});
  check(restored.get_stations_locations() == survey.geographic_locations(), "Restoring original coordinates failed");
  restored.ensure_utm_coordinates();
  check(restored.coordinates().utm && !restored.coordinates().centered, "Automatic native coordinates failed");
}
}
int main(int argc, char **argv)
{
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  try {
    QTemporaryDir dir;
    check(dir.isValid(), "Cannot create temporary directory");
    coordinate_checks(dir.path());
    // Optional existing project for checking the pre-coordinate archive layout.
    if(argc == 2) {
      std::shared_ptr<MTSurveyData> old_survey;
      const auto read_project = [&](bool has_header) {
        std::ifstream file(argv[1], std::ios::binary);
        boost::archive::binary_iarchive archive(file);
        if(has_header) {
          std::uint32_t magic = 0, version = 0;
          archive >> magic >> version;
          check(magic == 0x45444954, "Legacy unmarked project");
        }
        archive >> old_survey;
      };
      try {read_project(true);} catch(const std::exception &) {read_project(false);}
      check(old_survey != nullptr && old_survey->n_stations() > 0 && !old_survey->coordinates().utm,
            "Existing project failed to load its original geographic coordinates");
      std::cout << "Existing project loaded with " << old_survey->n_stations() << " stations.\n";
    }
    const auto source = dir.path() + "/input.gofem";
    {
      std::ofstream input(source.toStdString());
      input << "RealZxy 1 Plane_wave S01 0.002 0.0001\nImagZxy 1 Plane_wave S01 0.002 0.0001\n";
    }
    MTSurveyData survey;
    survey.load_from_gofem(source.toStdString());
    survey.get_station_data("S01").set_error_floor(0.);
    survey.get_station_data("S01").set_position({{55., 12., 320.}});
    SurveyCoordinatesDialog coordinates_dialog(survey);
    widget<QSpinBox>(coordinates_dialog, "utmZone")->setValue(32);
    widget<QCheckBox>(coordinates_dialog, "utmCenter")->setChecked(true);
    auto *origin_e = widget<QLineEdit>(coordinates_dialog, "utmOriginEasting");
    auto *origin_n = widget<QLineEdit>(coordinates_dialog, "utmOriginNorthing");
    check(origin_e->isReadOnly() && std::abs(origin_e->text().toDouble()-691875.632137542) < 1e-6, "Copyable easting wrong");
    check(std::abs(origin_n->text().toDouble()-6098907.825129169) < 1e-6, "Copyable northing wrong");
    widget<QPushButton>(coordinates_dialog, "copyUtmOrigin")->click();
    check(QApplication::clipboard()->text().contains(origin_e->text()) &&
          QApplication::clipboard()->text().contains(origin_n->text()), "Copy origin omitted values");
    coordinates_dialog.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->click();
    check(survey.coordinates().utm && survey.coordinates().centered, "Coordinate action did not apply");
    check(survey.geographic_locations()[0] == std::array<double,3>{{55.,12.,320.}}, "Lost original latitude/longitude");
    check(survey.get_stations_locations()[0] == std::array<double,3>{{0.,0.,320.}}, "Survey not centered");
    ExportNativeMTDialog dialog(survey, dir.path());
    auto *types = widget<QListWidget>(dialog, "nativeTypes");
    check(types->count() == 24, "Dialog must offer all 24 scalar types");
    for(int i = 0; i < types->count(); ++i) {
      auto *item = types->item(i);
      check((item->checkState() == Qt::Checked) == (i < 8), "Unexpected default component selection");
      const int type = item->data(Qt::UserRole).toInt();
      item->setCheckState(type == RealZxy || type == ImagZxy || type == PhsZxy || type == RhoZxy ? Qt::Checked : Qt::Unchecked);
    }
    check(widget<QListWidget>(dialog, "nativePeriods")->selectedItems().size() == 1, "Default period selection");
    check(dialog.findChild<QLineEdit *>("nativeReceiverFile") == nullptr, "Export still asks for coordinates");
    check(!widget<QCheckBox>(dialog, "nativeConjugate")->isChecked(), "Responses changed by default");
    const auto destination = dir.path() + "/survey.data";
    QString failure;
    bool saved = false, file_selected = false;
    QString last_modal;
    QTimer automation;
    QObject::connect(&automation, &QTimer::timeout, [&] {
      auto *active = QApplication::activeModalWidget();
      const QString modal = active ? QString(active->metaObject()->className()) + ": " + active->windowTitle() : "none";
      if(modal != last_modal) {
        std::cerr << "Export workflow: " << modal.toStdString() << '\n';
        last_modal = modal;
      }
      if(auto *message = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
        if(message->windowTitle() != "MT export complete") failure = message->text();
        message->button(QMessageBox::Ok)->click();
        if(!failure.isEmpty()) dialog.reject();
      } else if(auto *file = qobject_cast<QFileDialog *>(QApplication::activeModalWidget())) {
        if(!file_selected) {
          file->selectFile(destination);
          file_selected = true;
        } else {
          // Let QFileSystemModel finish loading the selected directory before
          // accepting. Re-selecting on every tick can restart that async load.
          QMetaObject::invokeMethod(file, "accept", Qt::DirectConnection);
          saved = true;
        }
      }
    });
    automation.start(20);
    QTimer::singleShot(0, [&] {
      dialog.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Save)->click();
    });
    QTimer::singleShot(10000, [&] {
      failure = "Dialog export timed out";
      if(auto *modal = qobject_cast<QDialog *>(QApplication::activeModalWidget())) modal->reject();
      dialog.reject();
    });
    const int code = dialog.exec();
    check(failure.isEmpty(), failure.toStdString().c_str());
    check(code == QDialog::Accepted && saved, "Export dialog did not complete");
    check(dialog.exportedFile() == destination, "Export destination lost");
    std::ifstream data(destination.toStdString());
    auto rows = NativeMT::read_observations(data);
    check(rows.size() == 4, "Wrong dialog data selection");
    check(rows[0].type == RealZxy && rows[0].value == .002 && rows[0].error == .0001, "SI impedance changed");
    check(rows[1].type == ImagZxy && rows[1].value == .002, "Stored convention changed");
    check(rows[3].type == PhsZxy && rows[3].value == 45., "Stored phase convention changed");
    std::ifstream coordinates((dir.path() + "/survey.recvs").toStdString());
    check_receiver_origin((dir.path() + "/survey.recvs").toStdString(), survey.coordinates());
    auto receivers = NativeMT::read_receivers(coordinates);
    check(receivers.at("S01") == std::array<double,3>{{0.,0.,-320.}}, "GUI coordinate conversion failed");
    std::ifstream frequencies((dir.path() + "/survey.freqs").toStdString());
    std::string line;
    std::getline(frequencies, line);
    check(line == "1", "Frequency companion missing or wrong");
    std::cout << "Native MT dialog export workflow passed.\n";
  } catch(const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}
