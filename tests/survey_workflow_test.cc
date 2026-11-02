#include "include/SurveyHistory.h"
#include "include/SurveyExport.h"
#include "ExportPreview.h"
#include "ExportNativeMTDialog.h"
#include "ExportGOFEMDialog.h"
#include "StationComparisonWindow.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
template <class T> T *control(QWidget &parent, const char *name) {
  auto *result = parent.findChild<T *>(name);
  check(result, "Missing workflow control");
  return result;
}
std::shared_ptr<MTSurveyData> fixture(const QString &directory) {
  const auto path = directory + "/survey.gofem";
  {
    std::ofstream file(path.toStdString());
    for (const auto *name : {"A", "B", "C"})
      for (double f : {1., 10.}) {
        file << "RealZxy " << f << " Plane_wave " << name << " .002 .0001\nImagZxy " << f << " Plane_wave "
             << name << " .003 .0002\n"
             << "RealTzx " << f << " Plane_wave " << name << " 0 .01\nImagTzx " << f << " Plane_wave " << name
             << " 0 .01\n";
      }
  }
  auto survey = std::make_shared<MTSurveyData>();
  survey->load_from_gofem(path.toStdString());
  survey->get_station_data("A").set_position({10., 179.9, 100.});
  survey->get_station_data("B").set_position({10., -179.9, 200.});
  survey->get_station_data("C").set_position({10., 179., 300.});
  return survey;
}
void historyChecks(const std::shared_ptr<MTSurveyData> &survey) {
  const auto *stationAddress = &survey->get_station_data("A");
  int restored = 0;
  SurveyHistory history([&] { ++restored; });
  history.reset(survey);
  const auto masks = survey->get_station_data("A").impedance_mask();
  check(!history.record("No change") && history.stack().count() == 0, "No-op pollutes history");
  survey->get_station_data("A").set_data_mask(RealZxy, 1., false, true);
  history.record("Masks");
  const auto masked = survey->get_station_data("A").impedance_mask();
  double v, e;
  survey->get_station_data("A").scalar_value(RealZxy, 1, v, e);
  const double originalError = e;
  survey->set_error_floor(.2);
  history.record("Errors");
  survey->set_coordinates(SurveyCoordinates::calculate(survey->geographic_locations(), 60, true, true));
  history.record("Coordinates");
  survey->rename_station("C", "RENAMED");
  history.record("Rename");
  survey->remove_station("B");
  history.record("Delete");
  survey->decimate();
  history.record("Decimate");
  check(history.stack().count() == 6, "Edit transaction count wrong");
  history.stack().undo();
  check(survey->get_station_data("A").frequencies().size() == 2, "Undo decimation lost data");
  check(&survey->get_station_data("A") == stationAddress, "Undo invalidated a surviving station reference");
  history.stack().undo();
  check(survey->is_station_present("B"), "Undo deletion failed");
  history.stack().undo();
  check(survey->is_station_present("C") && !survey->is_station_present("RENAMED"), "Undo rename failed");
  history.stack().undo();
  check(!survey->coordinates().utm, "Undo coordinates failed");
  history.stack().undo();
  survey->get_station_data("A").scalar_value(RealZxy, 1, v, e);
  check(e == originalError && survey->get_station_data("A").impedance_mask() == masked,
        "Undo error floor changed masks/errors");
  history.stack().undo();
  check(survey->get_station_data("A").impedance_mask() == masks && restored == 6, "Undo masks failed");
  history.stack().redo();
  check(survey->get_station_data("A").impedance_mask() == masked, "Redo masks failed");
  survey->set_active_flag("B", false);
  history.record("Disable");
  check(!history.stack().canRedo() && history.stack().count() == 2, "New edit did not discard redo branch");
  const auto original = history.original();
  check(original->is_active("B") && original->get_station_data("A").impedance_mask() == masks,
        "Original survey mutated");
  for (int i = 0; i < 25; ++i) {
    survey->set_active_flag("B", i % 2 == 0);
    history.record("Activity");
  }
  check(history.stack().count() == 20, "History limit ignored");
  *survey = *original;
  history.reset(survey);
  check(!history.stack().canUndo(), "Reload kept previous history");
  auto grid = std::make_shared<Topography::Grid>();
  survey->topography().grid = grid;
  survey->topography().style.opacity = .4;
  history.record("Terrain");
  history.stack().undo();
  check(!survey->topography().grid, "Undo terrain failed");
  history.stack().redo();
  check(survey->topography().grid == grid && survey->topography().style.opacity == .4,
        "Redo copied or changed immutable terrain");
  survey->restore(*original);
  history.reset(survey);
}
void exportChecks(const std::shared_ptr<MTSurveyData> &survey) {
  MTStationData phaseOnly;
  phaseOnly.set_name("PHASE");
  phaseOnly.set_size(1, true);
  phaseOnly.set_frequencies({1.});
  phaseOnly.set_data(1., {PhsZxy}, {42.}, {2.});
  NativeMT::Options phaseOptions;
  phaseOptions.types = {PhsZxy};
  phaseOptions.periods = {1.};
  phaseOptions.horizontal_axes_confirmed = phaseOptions.conventions_confirmed = true;
  phaseOptions.coordinate_description = "Synthetic Cartesian coordinates";
  phaseOptions.receivers["PHASE"] = {{0., 0., 0.}};
  check(NativeMT::prepare({&phaseOnly}, phaseOptions).observations.front().value == 42.,
        "Explicit phase-only data lost on export");
  survey->get_station_data("A").set_data_mask(RealZxy, 1., false);
  survey->set_active_flag("C", false);
  // Disabled stations with invalid coordinates must not block an otherwise valid export.
  survey->get_station_data("C").set_position({NAN, NAN, NAN});
  NativeMT::Options options;
  options.types = {RealZxy, ImagZxy, RealTzx, ImagTzx};
  options.periods = {.1, 1.};
  const auto native = SurveyExport::prepare(*survey, options, SurveyExport::Format::Native);
  const auto gofem = SurveyExport::prepare(*survey, options, SurveyExport::Format::GoFEM);
  check(!survey->coordinates().utm, "Preview changed saved projection");
  check(native.observations.size() == 14 && native.masked == 2 && native.receivers.size() == 2,
        "Export selection or counts wrong");
  check(native.receivers.at("A")[2] == -100. && gofem.receivers.at("A")[2] == 100.,
        "Vertical export conventions wrong");
  std::istringstream input(gofem.data_text);
  const auto rows = MTResponseData::read(input, MTResponseData::Format::GoFEM);
  check(rows.size() == gofem.observations.size(), "GoFEM preview differs from written rows");
  for (size_t i = 0; i < rows.size(); ++i)
    check(rows[i].value == gofem.observations[i].value && rows[i].error == gofem.observations[i].error &&
              rows[i].type == gofem.observations[i].type,
          "GoFEM writer changed prepared values");
  ExportPreview preview;
  check(preview.update(*survey, options, SurveyExport::Format::Native), "Valid preview rejected");
  check(control<QTableWidget>(preview, "exportStations")->rowCount() == 2 &&
            control<QTableWidget>(preview, "exportPeriods")->rowCount() == 2,
        "Preview tables wrong");
  check(control<QLabel>(preview, "exportDetails")->text().contains("UTM origin"), "Preview omitted origin");
  auto bad = options;
  bad.periods.clear();
  check(!preview.update(*survey, bad, SurveyExport::Format::Native) && preview.result().observations.empty(),
        "Invalid selection left stale prepared output");
  const auto artifact = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  ExportNativeMTDialog dialog(*survey, {});
  dialog.show();
  QApplication::processEvents();
  if (!artifact.isEmpty())
    dialog.grab().save(artifact + "/export-preview.png");
  ExportGOFEMDialog legacy(*survey);
  legacy.show();
  QApplication::processEvents();
  check(legacy.exportPreview().result().observations.size() == 6, "GoFEM dialog preview incorrect");
  if (!artifact.isEmpty())
    legacy.grab().save(artifact + "/gofem-preview.png");
}
void comparisonChecks(const std::shared_ptr<MTSurveyData> &survey) {
  StationComparisonWindow window;
  window.setData(survey, {"A"});
  window.show();
  QApplication::processEvents();
  auto *count = control<QSpinBox>(window, "compareCount");
  count->setValue(1);
  control<QPushButton>(window, "compareNearest")->click();
  auto *list = control<QListWidget>(window, "compareStations");
  check(list->item(0)->data(Qt::UserRole).toString() == "A" &&
            list->item(1)->data(Qt::UserRole).toString() == "B" && list->item(1)->checkState() == Qt::Checked,
        "Date-line nearest neighbour wrong");
  auto *plot = control<QCustomPlot>(window, "compareFirstPlot");
  check(plot->graphCount() == 2, "Checked station curves missing");
  check(plot->graph(0)->dataCount() == 2 && std::isnan(plot->graph(0)->data()->at(1)->value),
        "Masked period is not a curve gap");
  list->item(1)->setCheckState(Qt::Unchecked);
  check(plot->graphCount() == 1, "Station checkbox did not hide curve");
  auto *map = control<QCustomPlot>(window, "compareMap");
  check(map->xAxis->range().size() < 2., "Date-line map discontinuity");
  QCPGraph *neighbour = nullptr;
  for (int i = 0; i < map->graphCount(); ++i)
    if (map->graph(i)->name() == "B")
      neighbour = map->graph(i);
  check(neighbour, "Neighbour is missing from map");
  const auto point = neighbour->data()->constBegin();
  const QPointF click(map->xAxis->coordToPixel(point->key), map->yAxis->coordToPixel(point->value));
  QMouseEvent press(QEvent::MouseButtonPress, click, Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
  QMouseEvent release(QEvent::MouseButtonRelease, click, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(map, &press);
  QApplication::sendEvent(map, &release);
  check(list->item(1)->checkState() == Qt::Checked && plot->graphCount() == 2,
        "Map click did not toggle station curves");
  control<QComboBox>(window, "compareFamily")->setCurrentIndex(1);
  check(plot->graph(0)->dataCount() == 2 && plot->graph(0)->data()->at(0)->value == 0.,
        "Valid zero tipper omitted");
  control<QCheckBox>(window, "compareErrors")->setChecked(true);
  check(plot->plottableCount() == 4, "Comparison errors missing");
  control<QComboBox>(window, "compareFamily")->setCurrentIndex(0);
  const auto artifact = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if (!artifact.isEmpty())
    window.grab().save(artifact + "/station-comparison.png");
  window.setData({});
  check(plot->graphCount() == 0, "Unloading leaves comparison curves");
}
void exportNameChecks(const QString &directory) {
  auto survey = fixture(directory);
  survey->rename_station("A", "001");
  survey->rename_station("B", "2B");
  survey->rename_station("C", "S_001");
  NativeMT::Options options;
  options.types = {RealZxy, ImagZxy}; options.periods = {.1, 1.};
  for(auto format: {SurveyExport::Format::Native, SurveyExport::Format::GoFEM}) {
    const auto original = SurveyExport::prepare(*survey, options, format);
    check(original.receivers.count("001") && original.receivers.count("2B") && original.name_mapping_text.empty(),
          "Digit-leading names must export unchanged by default");
    options.compatible_names = true;
    const auto aliased = SurveyExport::prepare(*survey, options, format);
    check(aliased.original_names.at("S_001_2") == "001" && aliased.original_names.at("S_2B") == "2B" &&
          aliased.receivers.count("S_001"), "Aliases collided with an existing station name");
    check(aliased.name_mapping_text.find("\"001\",\"S_001_2\"") != std::string::npos,
          "Name mapping companion missing");
    std::istringstream data(aliased.data_text), receivers(aliased.receiver_text);
    const auto rows = MTResponseData::read(data, format == SurveyExport::Format::Native ?
                                         MTResponseData::Format::Native : MTResponseData::Format::GoFEM);
    const auto positions = NativeMT::read_receivers(receivers);
    check(rows.size() == original.observations.size(), "Aliasing lost observations");
    for(const auto &row: rows) {
      const auto it = aliased.original_names.find(row.receiver);
      const auto name = it == aliased.original_names.end() ? row.receiver : it->second;
      check(positions.at(row.receiver) == original.receivers.at(name), "Receiver alias changed coordinates");
      bool found = false;
      for(const auto &before: original.observations)
        if(before.receiver == name && before.frequency == row.frequency && before.type == row.type)
          found = before.value == row.value && before.error == row.error;
      check(found, "Aliasing changed data or uncertainties");
    }
    options.compatible_names = false;
  }
  survey->rename_station("2B", "2 B,\"x\"");
  ExportNativeMTDialog dialog(*survey, {});
  auto *preview = control<ExportPreview>(dialog, "exportPreview");
  auto *save = dialog.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Save);
  check(!save->isEnabled() && control<QLabel>(dialog, "exportSummary")->text().contains("Use compatible names"),
        "Invalid names should offer a remedy");
  auto *names = control<QCheckBox>(dialog, "exportCompatibleNames");
  names->setChecked(true); QApplication::processEvents();
  check(save->isEnabled() && preview->result().original_names.at("S_2_B__x_") == "2 B,\"x\"",
        "Compatible-name option did not unblock export");
  check(preview->result().name_mapping_text.find("\"2 B,\"\"x\"\"\"") != std::string::npos,
        "Mapping CSV did not escape the original name");
  check(survey->is_station_present("001") && survey->is_station_present("2 B,\"x\""),
        "Export preview renamed project stations");
  const auto path = directory + "/aliased.data";
  QTimer::singleShot(0, [&] {
    auto *chooser = qobject_cast<QFileDialog *>(QApplication::activeModalWidget());
    check(chooser, "Export file chooser missing");
    chooser->selectFile(path);
    QMetaObject::invokeMethod(chooser, "accept", Qt::DirectConnection);
  });
  check(saveSurveyExport(&dialog, preview->result(), SurveyExport::Format::Native, directory) == path,
        "Aliased export was not saved");
  for(const auto &entry: std::map<QString, std::string>{{path, preview->result().data_text},
        {directory + "/aliased.recvs", preview->result().receiver_text},
        {directory + "/aliased.freqs", preview->result().frequency_text},
        {directory + "/aliased.names.csv", preview->result().name_mapping_text}}) {
    QFile file(entry.first);
    check(file.open(QIODevice::ReadOnly) && file.readAll().toStdString() == entry.second,
          "Saved export companions differ from preview");
  }
  names->setChecked(false); QApplication::processEvents();
  check(!save->isEnabled() && preview->result().observations.empty(), "Unchecked conversion left stale export");
  ExportGOFEMDialog gofem(*survey);
  control<QCheckBox>(gofem, "exportCompatibleNames")->setChecked(true); QApplication::processEvents();
  check(gofem.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->isEnabled() &&
        gofem.exportPreview().result().original_names.size() == 2, "GoFEM alias control did not refresh preview");
}
} // namespace
int main(int argc, char **argv) {
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  try {
    QTemporaryDir directory;
    check(directory.isValid(), "Temporary directory failed");
    auto survey = fixture(directory.path());
    historyChecks(survey);
    exportChecks(survey);
    comparisonChecks(survey);
    exportNameChecks(directory.path());
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cout << "Survey history, comparison and export preview checks passed.\n";
}
