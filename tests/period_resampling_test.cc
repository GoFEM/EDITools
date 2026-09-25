#include "PeriodLayoutWindow.h"
#include "include/MTMapData.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QTabWidget>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <cmath>
#include <fstream>
#include <sstream>

namespace {
using namespace PeriodResampling;
void check(bool value, const char *message) { if(!value) throw std::runtime_error(message); }
void near(double a, double b) { check(std::abs(a - b) < 1e-10 * std::max({1., std::abs(a), std::abs(b)}), "Resampling numeric mismatch"); }
template<class T> T *widget(QWidget &parent, const char *name) {
  auto *result = parent.findChild<T *>(name); check(result, name); return result;
}
unsigned index(const MTStationData &station, double p) {
  const int result = MTMapData::nearest_period(station.frequencies(), p, 1e-13);
  check(result >= 0, "Target period absent from common layout"); return result;
}
double scalar(const MTStationData &station, RealDataType type, double p, double expectedError = -1.) {
  double value, error;
  check(station.scalar_value(type, index(station, p), value, error), "Expected valid resampled scalar");
  if(expectedError >= 0.) near(error, expectedError);
  return value;
}
Status status(const MTSurveyData &survey, const std::string &name, unsigned component, double p) {
  for(const auto &r: survey.resampling_info().records)
    if(r.station == name && r.component == component && std::abs(r.period - p) < p * 1e-12) return r.status;
  throw std::runtime_error("Missing audit record");
}
template<class F> void rejects(F action) {
  try { action(); } catch(const std::invalid_argument &) { return; }
  throw std::runtime_error("Invalid resampling settings accepted");
}
}

void period_resampling_tests(const QString &directory)
{
  using namespace PeriodResampling;
  const auto path = directory + "/resampling-input.gofem";
  std::ofstream rows(path.toStdString());
  auto write = [&](const char *station, double frequency, double z, double imaginary, double error) {
    for(auto type: {RealZxx, ImagZxx, RealZyy, ImagZyy}) rows << Datum::convert_type_to_string(type) << ' ' << frequency << " PW " << station << " 0 " << error << '\n';
    rows << "RealZxy " << frequency << " PW " << station << ' ' << z << ' ' << error << '\n';
    rows << "ImagZxy " << frequency << " PW " << station << ' ' << imaginary << ' ' << error << '\n';
    rows << "RealZyx " << frequency << " PW " << station << ' ' << -z << ' ' << error << '\n';
    rows << "ImagZyx " << frequency << " PW " << station << ' ' << -2*z << ' ' << error << '\n';
    for(auto type: {RealTzx, ImagTzx, RealTzy, ImagTzy}) rows << Datum::convert_type_to_string(type) << ' ' << frequency << " PW " << station << " 0 " << error << '\n';
  };
  write("A", 1., 1., 1., .1); write("A", .25, 3., 5., .3); write("A", .01, 8., 9., .5);
  write("B", .25, 2., 3., .2);
  write("D", 1., 1., 1., .1); write("D", .25, 3., 5., .3);
  for(auto type: {RealTzx, ImagTzx, RealTzy, ImagTzy}) rows << Datum::convert_type_to_string(type) << " .5 PW D 0 .7\n";
  for(auto type: {PTxx, PTxy, PTyx, PTyy}) {
    rows << Datum::convert_type_to_string(type) << " 1 PW C 2 .1\n";
    rows << Datum::convert_type_to_string(type) << " .25 PW C 4 .3\n";
  }
  rows.close();
  auto survey = std::make_shared<MTSurveyData>(); survey->load_from_gofem(path.toStdString());
  survey->get_station_data("A").set_position({{55., 9., 100.}});
  survey->get_station_data("B").set_position({{55.01, 9., 100.}});
  survey->get_station_data("C").set_position({{55.02, 9., 100.}});
  survey->get_station_data("D").set_position({{55.03, 9., 100.}});
  const std::vector<double> periods{.5, 1., 2., 4., 10., 100., 200.};
  auto output = resample(*survey, periods, {4., true});
  const auto &a = output->get_station_data("A");
  near(scalar(a, RealZxy, 2., .2), 2.); near(scalar(a, ImagZxy, 2., .2), 3.);
  near(scalar(a, PTxx, 2.), 2.); near(scalar(a, PTyy, 2.), 1.5);
  near(scalar(a, PhsZxy, 2.), std::atan2(3., 2.) * 180. / std::acos(-1.));
  near(scalar(a, RhoZxy, 2.), 13. / (4e-7 * std::acos(-1.) * std::acos(-1.)));
  near(scalar(a, RealTzx, 2.), 0.);
  near(scalar(output->get_station_data("D"), RealZxy, 2., .2), 2.);
  near(scalar(output->get_station_data("D"), RealTzx, 2., .7), 0.);
  check(status(*output, "D", 1, 2.) == Status::Interpolated && status(*output, "D", 4, 2.) == Status::Existing,
        "Existing component at a period prevented filling another component's small gap");
  near(scalar(output->get_station_data("C"), PTxx, 2., .2), 3.);
  check(status(*output, "A", 1, 2.) == Status::Interpolated && status(*output, "A", 6, 2.) == Status::Derived,
        "Wrong interpolation/derivation provenance");
  check(status(*output, "A", 1, .5) == Status::Outside && status(*output, "A", 1, 200.) == Status::Outside,
        "Resampling extrapolated");
  check(status(*output, "A", 1, 10.) == Status::Gap, "Resampling filled a large gap");
  check(status(*output, "B", 1, 4.) == Status::Existing && status(*output, "B", 1, 2.) == Status::Outside,
        "Single-point component was lost or extrapolated");
  for(const auto &name: output->get_stations_names()) check(output->get_station_data(name).frequencies() == a.frequencies(), "Station grids are not common");
  check(!a.impedance_mask()[1][index(a, 10.)], "Unsupported target entered fitting/export masks");
  check(survey->get_station_data("A").frequencies().size() == 3 && output->response_observations().empty(), "Source changed or stale response rows survived");
  check(output->resampling_source() && output->resampling_source()->get_station_data("A").position() == a.position(), "Source snapshot or locations lost");
  auto &originalA = survey->get_station_data("A");
  originalA.set_data_mask(RealZxy, .25, false);
  auto blocked = resample(*survey, {1., 2., 4., 8., 100.}, {100., true});
  check(status(*blocked, "A", 1, 2.) == Status::Masked && status(*blocked, "A", 1, 4.) == Status::Masked,
        "Masked sample was replaced or bridged");
  auto copied = blocked->get_station_data("A"); copied.set_data_mask(RealZxy, .25, true);
  near(scalar(copied, RealZxy, 4., .3), 3.); // Exact masked samples retain their original values/errors.
  auto bridge = resample(*survey, {2., 4.}, {100., false});
  check(status(*bridge, "A", 1, 2.) == Status::Interpolated && status(*bridge, "A", 1, 4.) == Status::Masked,
        "Optional mask bridging changed exact masked samples");
  check(status(*resample(*survey, {1.01}, {4., false}), "A", 1, 1.01) == Status::Gap,
        "Near-endpoint target bypassed the bounding-period gap limit");
  originalA.set_data_mask(RealZxy, .25, true);
  originalA.set_data_mask(PTxx, .25, false);
  check(status(*resample(*survey, {2.}, {4., true}), "A", 6, 2.) == Status::Masked, "Tensor derivation bypassed its masks");
  originalA.set_data_mask(PTxx, .25, true);
  originalA.set_data(.25, {RealTzy}, {NAN}, {.3});
  check(!available(status(*resample(*survey, {2.}, {4., true}), "A", 5, 2.)), "Incomplete complex endpoint produced a tipper");
  originalA.set_data(.25, {RealTzy}, {0.}, {.3});
  originalA.set_data(.25, {RealTzx}, {0.}, {NAN});
  check(status(*resample(*survey, {2.}, {4., true}), "A", 4, 2.) == Status::InvalidError, "Invalid endpoint uncertainty accepted");
  originalA.set_data(.25, {RealTzx}, {0.}, {.3});
  survey->set_active_flag("B", false);
  auto disabled = resample(*survey, {4.}, {4., true});
  check(!disabled->is_active("B"), "Disabled station re-enabled");
  check(PeriodResampling::layout(disabled->get_station_data("B"), 1).front().status == Status::Masked, "Disabled station counted as coverage");
  survey->set_active_flag("B", true);
  std::stringstream saved;
  { boost::archive::binary_oarchive archive(saved); archive << *output; }
  MTSurveyData restored;
  { boost::archive::binary_iarchive archive(saved); archive >> restored; }
  near(scalar(restored.get_station_data("A"), RealZxy, 2.), 2.);
  near(scalar(restored.resampling_source()->get_station_data("A"), RealZxy, 4.), 3.);
  check(restored.resampling_info().records.size() == output->resampling_info().records.size(), "Resampling audit lost on save/load");
  near(logarithmic_grid(1., 100., 2)[1], std::sqrt(10.));
  check(validate_periods({4., 1., 4.}) == std::vector<double>({1., 4.}), "Custom sorting/deduplication failed");
  check(resample(*survey, {1e-8, 1.1e-8}, {})->get_unique_periods().size() == 2,
        "Common layout collapsed distinct short-period samples");
  rejects([&] { resample(*survey, {}, {}); }); rejects([&] { resample(*survey, {NAN}, {}); });
  rejects([&] { resample(*survey, {0.}, {}); }); rejects([&] { resample(*survey, {1.}, {1., true}); });
  rejects([&] { logarithmic_grid(10., 1., 6); }); rejects([&] { logarithmic_grid(1e-20, 1e20, 100); });

  std::shared_ptr<MTSurveyData> applied;
  PeriodLayoutWindow window(nullptr, [&](std::shared_ptr<MTSurveyData> data) { applied = std::move(data); });
  window.setData(survey); window.show(); QApplication::processEvents();
  widget<QComboBox>(window, "layoutGridMode")->setCurrentIndex(2);
  widget<QLineEdit>(window, "layoutCustomPeriods")->setText("4, 2; 1 2 .5 10 100 200");
  widget<QDoubleSpinBox>(window, "layoutMaximumRatio")->setValue(4.);
  auto *apply = widget<QPushButton>(window, "layoutApply");
  check(!apply->isEnabled(), "Unpreviewed settings can be applied");
  widget<QPushButton>(window, "layoutPreview")->click();
  check(apply->isEnabled() && widget<QTableWidget>(window, "layoutCounts")->item(0, 4)->text() == "1", "Preview count disagrees with backend");
  apply->click(); check(applied && applied->resampling_info().periods == periods, "Apply lost the preview grid");
  near(scalar(applied->get_station_data("A"), RealZxy, 2.), 2.);
  widget<QComboBox>(window, "layoutGroup")->setCurrentIndex(1);
  check(apply->isEnabled(), "Display filter invalidated a valid preview");
  widget<QDoubleSpinBox>(window, "layoutMaximumRatio")->setValue(2.);
  check(!apply->isEnabled(), "Changing gap limit left a stale preview applicable");
  widget<QLineEdit>(window, "layoutCustomPeriods")->setText("1 invalid 3");
  widget<QPushButton>(window, "layoutPreview")->click(); check(!apply->isEnabled(), "Invalid custom input applied");
  widget<QComboBox>(window, "layoutGridMode")->setCurrentIndex(1);
  widget<QComboBox>(window, "layoutReference")->setCurrentText("B");
  widget<QPushButton>(window, "layoutPreview")->click(); apply->click();
  check(applied->resampling_info().periods == std::vector<double>{4.}, "Reference station grid not used");
  window.setData(output); widget<QPushButton>(window, "layoutOpenSource")->click();
  check(applied->get_station_data("A").frequencies().size() == 3, "Cannot reopen retained source");
  window.setData({}); widget<QPushButton>(window, "layoutPreview")->click(); check(!apply->isEnabled(), "Empty survey allowed resampling");
}
