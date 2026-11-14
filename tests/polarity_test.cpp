#include "core/PolarityAnalysis.h"
#include "core/SurveyExport.h"
#include "ui/workspace/SurveyHistory.h"
#include "ui/analysis/PolarityWindow.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QPushButton>
#include <QTableWidget>
#include <QTemporaryDir>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool value, const char *message) { if(!value) throw std::runtime_error(message); }
void checkNear(double a, double b) { check(std::isfinite(a) && std::abs(a - b) < 1e-9 * std::max(1., std::abs(b)), "Numerical mismatch"); }
double value(const MTStationData &s, RealDataType type, unsigned i = 0) {
  double v, e; check(s.scalar_value(type, i, v, e), "Missing scalar"); return v;
}
std::string bytes(const MTStationData &s) {
  std::ostringstream out; { boost::archive::binary_oarchive archive(out); archive << s; } return out.str();
}
MTStationData station() {
  MTStationData s; s.set_size(3, true); s.set_frequencies({1., 10., 100.}); s.set_name("TEST");
  for(double f : s.frequencies()) s.set_data(f,
    {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy, RealTzx, ImagTzx, RealTzy, ImagTzy},
    {.3, .2, 2., 3., -4., -2., .4, -.3, .1, .2, -.3, .4}, std::vector<double>(12, .01));
  s.calculate_apparent_resistivity(); s.calculate_phase(); s.calculate_phase_tensor(); s.set_error_floor(.1);
  return s;
}
void transformations() {
  const auto original = station();
  const RealDataType zr[] = {RealZxx, RealZxy, RealZyx, RealZyy}, zi[] = {ImagZxx, ImagZxy, ImagZyx, ImagZyy};
  const RealDataType rho[] = {RhoZxx, RhoZxy, RhoZyx, RhoZyy}, phs[] = {PhsZxx, PhsZxy, PhsZyx, PhsZyy};
  const RealDataType pt[] = {PTxx, PTxy, PTyx, PTyy};
  for(unsigned bits = 0; bits < 32; ++bits) {
    std::array<bool, 5> channels; for(unsigned i = 0; i < 5; ++i) channels[i] = bits & (1 << i);
    auto s = original;
    s.set_data_mask(RealZxx, 10., false); s.set_data_mask(RealTzy, 100., false); s.set_data_mask(PTxy, 1., false);
    const auto zm = s.impedance_mask(), tm = s.tipper_mask(), pm = s.phase_tensor_mask();
    s.reverse_polarity(channels);
    check(s.impedance_mask() == zm && s.tipper_mask() == tm && s.phase_tensor_mask() == pm, "Correction changed masks");
    s.set_data_mask(RealZxx, 10., true); s.set_data_mask(RealTzy, 100., true); s.set_data_mask(PTxy, 1., true);
    auto derived = s; derived.calculate_phase_tensor(); derived.calculate_phase();
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned c = 0; c < 4; ++c) {
        const int sign = channels[c / 2] != channels[2 + c % 2] ? -1 : 1;
        checkNear(value(s, zr[c], i), sign * value(original, zr[c], i));
        checkNear(value(s, zi[c], i), sign * value(original, zi[c], i));
        checkNear(value(s, rho[c], i), value(original, rho[c], i));
        checkNear(value(s, pt[c], i), value(derived, pt[c], i));
        checkNear(std::remainder(value(s, phs[c], i) - value(derived, phs[c], i), 360.), 0.);
        for(auto type : {zr[c], zi[c], rho[c], phs[c], pt[c]}) {
          double v, before, after; original.scalar_value(type, i, v, before); s.scalar_value(type, i, v, after); checkNear(after, before);
        }
      }
      for(unsigned c = 0; c < 2; ++c) for(auto type : (c == 0 ? std::vector<RealDataType>{RealTzx, ImagTzx} : std::vector<RealDataType>{RealTzy, ImagTzy})) {
        const int sign = channels[4] != channels[2 + c] ? -1 : 1;
        checkNear(value(s, type, i), sign * value(original, type, i));
        double v, before, after; original.scalar_value(type, i, v, before); s.scalar_value(type, i, v, after); checkNear(after, before);
      }
    }
    s.reverse_polarity(channels);
    for(auto type : {RealZxy, ImagZyx, RealTzx, ImagTzy, PTxy, PTyx, PhsZxy, PhsZyx}) checkNear(value(s, type), value(original, type));
  }
  MTStationData sparse; sparse.set_size(1, true); sparse.set_frequencies({1.});
  sparse.set_data(1., {PhsZxy, RhoZxy, PTxy, RealTzx}, {45., 100., .3, .2}, {2., 5., .01, .01});
  sparse.reverse_polarity({{true, false, false, false, false}});
  checkNear(value(sparse, PhsZxy), -135.); checkNear(value(sparse, RhoZxy), 100.); checkNear(value(sparse, PTxy), .3); checkNear(value(sparse, RealTzx), .2);
  sparse.reverse_polarity({{false, false, true, false, false}});
  checkNear(value(sparse, PTxy), -.3); checkNear(value(sparse, RealTzx), -.2);
  MTStationData restored; std::istringstream input(bytes(sparse)); boost::archive::binary_iarchive archive(input); archive >> restored;
  check(bytes(sparse) == bytes(restored), "Project serialization lost correction");
}
std::shared_ptr<MTSurveyData> fixture(const QString &directory, bool phaseOnly = false, double frequencyRatio = 1., bool seam = false) {
  const auto path = directory + "/survey.data";
  {
    std::ofstream out(path.toStdString());
    for(unsigned n = 0; n < 7; ++n) for(double f : {.1, 1., 10., 100.}) {
      const auto name = std::string(1, char('A' + n));
      const double freq = f * (n > 1 ? frequencyRatio : 1.);
      if(phaseOnly) {
        const double xy = seam ? (n % 2 ? -179. : 179.) : 45.;
        out << "PhsZxy " << freq << " Plane_wave " << name << ' ' << xy << " 2\nPhsZyx " << freq << " Plane_wave " << name << " -135 2\n";
      } else {
        const auto s = station();
        for(auto type : {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy, RealTzx, ImagTzx, RealTzy, ImagTzy})
          out << Datum::convert_type_to_string(type) << ' ' << freq << " Plane_wave " << name << ' ' << value(s, type) << " .01\n";
      }
    }
  }
  auto survey = std::make_shared<MTSurveyData>(); survey->load_from_gofem(path.toStdString());
  for(unsigned n = 0; n < 7; ++n) survey->get_station_data(std::string(1, char('A' + n))).set_position({{40. + .01 * n, -115., 100.}});
  return survey;
}
void detection(const QString &directory) {
  using PolarityAnalysis::analyse;
  auto survey = fixture(directory);
  auto results = analyse(*survey);
  for(const auto &s : results) for(const auto &c : s.components) check(c.supported() && !c.candidate(), "Healthy station falsely flagged");
  auto &target = survey->get_station_data("A");
  target.reverse_polarity({{true, false, false, false, false}});
  results = analyse(*survey);
  check(results[0].components[0].candidate() && !results[0].components[1].candidate(), "Single electric reversal missed");
  for(unsigned i = 1; i < results.size(); ++i) check(!results[i].components[0].candidate(), "Reversed neighbour contaminated healthy stations");
  target.reverse_polarity({{false, true, false, false, false}});
  results = analyse(*survey);
  check(results[0].components[0].candidate() && results[0].components[1].candidate(), "Both electric reversals missed");
  check(!analyse(*survey, {5, .01, 30.})[0].components[0].supported(), "Radius ignored");
  target.set_position({{NAN, NAN, NAN}});
  check(!analyse(*survey)[0].components[0].supported(), "Missing location used");
  target.set_position({{40., -115., 100.}});
  target.set_data_mask(RealZxy, .1, false); target.set_data_mask(RealZxy, 1., false);
  check(!analyse(*survey)[0].components[0].candidate(), "Too few periods flagged");
  target.set_active(false);
  check(analyse(*survey)[0].components[0].usable == 0, "Disabled target used");
  survey = fixture(directory, false, 1.01);
  check(!analyse(*survey)[0].components[0].supported(), "Compared unmatched periods");
  survey = fixture(directory, false, 1.0005);
  check(analyse(*survey)[0].components[0].supported(), "Nearby periods not matched");
  survey->set_error_floor(1.);
  check(!analyse(*survey)[0].components[0].supported(), "Low-SNR samples used");
  survey = fixture(directory, true, 1., true);
  check(analyse(*survey)[0].components[0].supported() && !analyse(*survey)[0].components[0].candidate(), "Circular phase seam mishandled");
  survey->get_station_data("A").reverse_polarity({{true, false, false, false, false}});
  check(analyse(*survey)[0].components[0].candidate(), "Phase-only reversal missed");
  survey = fixture(directory, true);
  for(const auto &name : survey->get_stations_names()) {
    auto &s = survey->get_station_data(name);
    for(double f : s.frequencies()) s.set_data(f, {PhsZxy}, {120.}, {2.});
  }
  check(!analyse(*survey)[0].components[0].candidate(), "Genuine shared out-of-quadrant phases flagged");
  for(const auto &name : survey->get_stations_names()) survey->get_station_data(name).reverse_polarity({{true, true, false, false, false}});
  check(!analyse(*survey)[0].components[0].candidate(), "Claimed detection of a common survey-wide reversal");
  for(const auto *name : {"B", "C", "D", "E"}) survey->set_active_flag(name, false);
  check(!analyse(*survey)[0].components[0].supported(), "Disabled neighbours used");
  survey = fixture(directory);
  survey->get_station_data("A").reverse_polarity({{false, false, true, false, false}});
  results = analyse(*survey);
  check(!results[0].components[0].candidate() && results[0].components[1].candidate(), "Magnetic-column reversal missed");
}
template<class T> T *control(QWidget &window, const char *name) {
  auto *result = window.findChild<T *>(name); check(result, "Missing polarity control"); return result;
}
void dialog(const QString &directory) {
  auto survey = fixture(directory);
  auto &target = survey->get_station_data("A"); target.reverse_polarity({{true, false, false, false, false}});
  target.set_error_floor(.05); target.set_data_mask(RealZxx, 1., false);
  const auto before = bytes(target), neighbour = bytes(survey->get_station_data("B"));
  int changed = 0; SurveyHistory history; history.reset(survey);
  PolarityWindow window(nullptr, history, [&] { ++changed; }); window.setData(survey, "A");
  window.show(); QApplication::processEvents();
  auto *apply = control<QPushButton>(window, "polarityApply");
  check(!apply->isEnabled(), "Correction preselected");
  control<QPushButton>(window, "polaritySuggest")->click();
  check(control<QCheckBox>(window, "polarityEx")->isChecked() && !control<QCheckBox>(window, "polarityEy")->isChecked(), "Suggestion wrong");
  check(bytes(target) == before && history.stack().count() == 0, "Preview edited survey");
  auto *plot = control<QCustomPlot>(window, "polarityXY");
  check(plot->graphCount() == 3 && plot->graph(2)->dataCount() == 4, "Preview missing");
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if(!artifacts.isEmpty()) window.grab().save(artifacts + "/polarity-review.png");
  apply->click();
  check(changed == 1 && history.stack().count() == 1 && !apply->isEnabled(), "Apply/history state wrong");
  check(value(target, RealZxy) > 0. && bytes(survey->get_station_data("B")) == neighbour, "Wrong station corrected");
  NativeMT::Options options; options.types = {RealZxy, ImagZxy}; options.periods = {.01, .1, 1., 10.};
  for(auto format : {SurveyExport::Format::Native, SurveyExport::Format::GoFEM}) {
    const auto prepared = SurveyExport::prepare(*survey, options, format);
    check(!prepared.observations.empty(), "Corrected export empty");
    for(const auto &row : prepared.observations) if(row.receiver == "A" && row.type == RealZxy) checkNear(row.value, 2.);
  }
  history.stack().undo(); check(bytes(target) == before, "Undo did not restore exact station");
  history.stack().redo(); check(value(target, RealZxy) > 0., "Redo failed");
  window.setData(survey, "B"); check(!apply->isEnabled(), "Stale channel selection after refresh");
}
}
int main(int argc, char **argv) {
  QApplication app(argc, argv); QTemporaryDir directory;
  try { transformations(); detection(directory.path()); dialog(directory.path()); }
  catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
  return 0;
}
