#include "ui/dialogs/EDIImportDialog.h"
#include "core/PeriodResampling.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QTimer>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
void check(bool ok, const char *message) { if(!ok) throw std::runtime_error(message); }
template<class T> T *control(QWidget &parent, const char *name) {
  auto *result = parent.findChild<T *>(name); check(result, name); return result;
}
template<class T> std::string bytes(const T &value) {
  std::ostringstream out; { boost::archive::binary_oarchive archive(out); archive << value; } return out.str();
}
double scalar(const MTStationData &station, double period, double *error = nullptr) {
  const auto &frequencies = station.frequencies();
  const auto f = std::find_if(frequencies.begin(), frequencies.end(), [period](double f) { return std::abs(f * period - 1.) < 1e-12; });
  double value, e; check(f != frequencies.end() && station.scalar_value(RealZxy, unsigned(f - frequencies.begin()), value, e), "Missing interpolated impedance");
  if(error) *error = e;
  return value;
}
PeriodResampling::Status status(const MTSurveyData &survey, unsigned component, double period) {
  for(const auto &record: survey.resampling_info().records)
    if(record.component == component && std::abs(record.period / period - 1.) < 1e-12) return record.status;
  throw std::runtime_error("Missing import interpolation audit");
}
}

void edi_import_interpolation_tests(const QString &directory) {
  const auto path = directory + "/import-layout.gofem";
  {
    std::ofstream out(path.toStdString()); out << std::setprecision(17);
    const std::map<std::string, std::vector<double>> layouts{{"OLD_A", {.5, 2., 3., 10.}}, {"OLD_B", {4., 200.}}, {"NEW", {1., 2., 4., 100.}}};
    for(const auto &entry: layouts) for(double period: entry.second) {
      for(auto type: {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy, RealTzx, ImagTzx, RealTzy, ImagTzy}) {
        const double value = type == RealZxy || type == ImagZxy ? period :
          type == RealZyx || type == ImagZyx ? -period : type == RealZxx || type == RealZyy || type == ImagZxx || type == ImagZyy ? 0. : .1;
        out << Datum::convert_type_to_string(type) << ' ' << 1. / period << " PW " << entry.first << ' ' << value << " .01\n";
      }
    }
  }
  MTSurveyData original; original.load_from_gofem(path.toStdString());
  original.get_station_data("NEW").set_error_floor(.2);
  original.get_station_data("NEW").set_data_mask(RealZxy, .5, false);
  original.get_station_data("OLD_A").set_data_mask(RealZxx, .5, false);
  original.set_active_flag("OLD_B", false); // Disabled stations still define the existing layout.
  const auto before = bytes(original);
  auto candidate = original;
  for(bool accept: {false, true}) {
    QTimer::singleShot(0, [&] {
      auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget()); check(dialog, "Import dialog missing");
      auto *mode = control<QComboBox>(*dialog, "ediPeriodMode"); mode->setCurrentIndex(2);
      check(control<QDoubleSpinBox>(*dialog, "ediMaximumRatio")->value() == 2., "Import gap limit differs from resampling default");
      check(control<QCheckBox>(*dialog, "ediMaskBarriers")->isChecked(), "Import crosses masks by default");
      auto *table = control<QTableWidget>(*dialog, "ediPeriodGroups");
      check(table->rowCount() == 6 && table->item(0, 0)->text() == "0.5" && table->item(5, 0)->text() == "200", "Existing-period union is wrong");
      check(table->item(2, 4)->toolTip().contains("mask barrier"), "Masked-component exclusion not explained");
      check(table->item(4, 4)->toolTip().contains("Gap too large"), "Gap exclusion not explained");
      control<QDoubleSpinBox>(*dialog, "ediMaximumRatio")->setValue(4.);
      control<QCheckBox>(*dialog, "ediMaskBarriers")->setChecked(false);
      check(table->item(2, 1)->text() == "0" && table->item(2, 2)->text() == "6" && table->item(2, 3)->text() == "4",
            "Preview does not report interpolated Z/T and derived PT separately");
      check(table->item(0, 4)->text() == "10", "Out-of-range values not excluded");
      check(bytes(candidate) == before, "Preview mutated the import candidate");
      const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
      if(accept && !artifacts.isEmpty()) {
        QApplication::processEvents();
        dialog->grab().save(artifacts + "/edi-interpolation.png");
      }
      if(accept) dialog->accept(); else dialog->reject();
    });
    check(reviewEDIImport(nullptr, candidate, {"NEW"}, {}) == accept, "Import accept/cancel result wrong");
    if(!accept) check(bytes(candidate) == before, "Cancel changed import data");
  }
  for(const auto *name: {"OLD_A", "OLD_B"}) check(bytes(candidate.get_station_data(name)) == bytes(original.get_station_data(name)), "Interpolation altered an existing station");
  check(candidate.get_survey_name() == original.get_survey_name(), "Interpolated import renamed project");
  check(candidate.get_station_data("NEW").frequencies().size() == 6, "Incoming station not on target layout");
  double e1, e4, actualError;
  const double a = scalar(original.get_station_data("NEW"), 1., &e1), b = scalar(original.get_station_data("NEW"), 4., &e4);
  const double weight = std::log(3.) / std::log(4.);
  const double result = scalar(candidate.get_station_data("NEW"), 3., &actualError);
  check(std::abs(result - ((1. - weight) * a + weight * b)) < 1e-12 && std::abs(actualError - ((1. - weight) * e1 + weight * e4)) < 1e-12, "Import interpolation/error propagation differs from shared engine");
  using S = PeriodResampling::Status;
  check(status(candidate, 1, .5) == S::Outside && status(candidate, 1, 200.) == S::Outside, "Import extrapolated");
  check(status(candidate, 1, 10.) == S::Gap && status(candidate, 1, 2.) == S::Masked, "Gap or exact mask lost");
  check(status(candidate, 4, 2.) == S::Existing && status(candidate, 4, 3.) == S::Interpolated, "Z masks affected tippers");
  check(bytes(*candidate.resampling_source()) == before, "Original import not retained");
  std::istringstream input(bytes(candidate)); MTSurveyData restored;
  { boost::archive::binary_iarchive archive(input); archive >> restored; }
  check(bytes(*restored.resampling_source()) == before && restored.resampling_info().records.size() == 60, "Saved project lost import source/audit");

  // No existing layout, no overlapping values and oversized grids all leave
  // Import disabled, and can be recovered by switching to Keep original.
  for(int scenario = 0; scenario < 3; ++scenario) {
    candidate = original;
    std::set<std::string> incoming{"NEW"};
    if(scenario == 0) incoming = {"NEW", "OLD_A", "OLD_B"};
    else {
      std::set<double> frequencies;
      for(int i = 0; i < (scenario == 1 ? 1 : 2001); ++i) frequencies.insert(10000. + i);
      for(const auto *name: {"OLD_A", "OLD_B"}) {
        auto &station = candidate.get_station_data(name); station.set_size(frequencies.size(), true); station.set_frequencies(frequencies);
      }
    }
    const auto unchanged = bytes(candidate);
    QTimer::singleShot(0, [&] {
      auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget()); check(dialog, "Import dialog missing");
      auto *mode = control<QComboBox>(*dialog, "ediPeriodMode"); mode->setCurrentIndex(2);
      auto *buttons = dialog->findChild<QDialogButtonBox *>(); check(buttons && !buttons->button(QDialogButtonBox::Ok)->isEnabled(), "Invalid interpolation was applicable");
      mode->setCurrentIndex(0); check(buttons->button(QDialogButtonBox::Ok)->isEnabled(), "Cannot recover to original periods");
      dialog->accept();
    });
    check(reviewEDIImport(nullptr, candidate, incoming, {}) && bytes(candidate) == unchanged, "Keep-original mode applied a stale interpolation preview");
  }
}
