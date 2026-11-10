#include "AutoMaskDialog.h"
#include "mainwindow.h"
#include "include/MTSurveyData.h"
#include "include/SurveyHistory.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QAction>
#include <QElapsedTimer>
#include <QListWidget>
#include <QProgressDialog>
#include <QSettings>
#include <QTabBar>
#include <QThread>
#include <QCheckBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool ok, const char *text) { if(!ok) throw std::runtime_error(text); }
std::string bytes(const MTSurveyData &survey) {
  std::ostringstream out; boost::archive::binary_oarchive archive(out); archive << survey; return out.str();
}
MTSurveyData fixture(const QString &directory, bool spike = true) {
  const auto path = directory + "/auto-mask.data";
  {
    std::ofstream out(path.toStdString()); out.precision(17);
    for(unsigned n = 0; n < 6; ++n) for(unsigned p = 0; p < 9; ++p) {
      const double f = std::pow(2., -double(p)); const char name = 'A' + n;
      const bool bad = spike && n == 0 && p == 4;
      const std::vector<RealDataType> types{RhoZxy, PhsZxy, RhoZyx, PhsZyx, RhoZxx, PhsZxx, RhoZyy, PhsZyy, RealTzx, ImagTzx, RealTzy, ImagTzy, PTxx, PTxy, PTyx, PTyy};
      const std::vector<double> values{100.*std::pow(2., .1*p)*(bad ? 10. : 1.), 45.+p+(bad ? 60. : 0.), 150., -135., 1., 45., 1., 45., .1, .2, .3, .4, 1., .1, .1, 1.};
      for(unsigned i = 0; i < types.size(); ++i) out << Datum::convert_type_to_string(types[i]) << ' ' << f << " Plane_wave " << name << ' ' << values[i] << ' ' << (i >= 8 && i < 12 ? .01 : 1.) << '\n';
    }
  }
  MTSurveyData survey; survey.load_from_gofem(path.toStdString());
  for(unsigned n = 0; n < 6; ++n) survey.get_station_data(std::string(1, 'A'+n)).set_position({{40.+.01*n, -115., 0.}});
  return survey;
}
const AutoMask::Point &point(const AutoMask::Result &result, double frequency = 1./16., RealDataType component = RealZxy) {
  for(const auto &p : result.points) if(p.station == "A" && p.frequency == frequency && p.component == component) return p;
  throw std::runtime_error("Expected assessment missing");
}
template<class T> T *control(QWidget &d, const char *name) {
  auto *result = d.findChild<T *>(name); check(result, "Missing auto-mask control"); return result;
}
void engine(const QString &directory) {
  auto survey = fixture(directory); const auto original = bytes(survey);
  AutoMask::Options options; options.criteria[0] = false; options.minimumFlags = 2;
  auto result = AutoMask::analyse(survey, {"A"}, options);
  check(point(result).candidate && point(result).votes() == 2, "Isolated spike not nominated by period and spatial checks");
  for(const auto &p : result.points) if(p.candidate) check(p.component == RealZxy && p.frequency == 1./16., "Healthy point nominated");
  check(!std::isfinite(point(result, 1.).scores[1]), "Endpoint extrapolated");
  for(unsigned k = 0; k < AutoMask::criterionCount; ++k) {
    options.criteria.fill(false); options.criteria[k] = true; options.minimumFlags = 1;
    result = AutoMask::analyse(survey, {"A"}, options);
    check(point(result).candidate == (k == 1 || k == 2), "Single criterion switch ignored or missing impedance used");
    for(unsigned j = 0; j < AutoMask::criterionCount; ++j) if(j != k) check(!std::isfinite(point(result).scores[j]) && !point(result).flagged[j], "Disabled criterion computed or voted");
  }
  options.criteria = {{false, false, false}};
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(!p.candidate, "All-disabled checks nominated a mask");
  options.criteria = {{false, true, true}}; options.minimumFlags = 3;
  check(!point(AutoMask::analyse(survey, {"A"}, options)).candidate, "Minimum agreement ignored");
  options.minimumFlags = 2; options.threshold.fill(100.);
  check(!point(AutoMask::analyse(survey, {"A"}, options)).candidate, "Thresholds ignored");
  options = {}; options.criteria[0] = false; options.components = {{false, true, false, false}};
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(p.component == RealZyx && !p.candidate, "Component toggle ignored");
  options.components[0] = true; options.maximumGapRatio = 1.5;
  check(!std::isfinite(point(AutoMask::analyse(survey, {"A"}, options)).scores[1]), "Large gap bridged");
  options.maximumGapRatio = 3.; options.radiusKm = .01;
  check(!std::isfinite(point(AutoMask::analyse(survey, {"A"}, options)).scores[2]), "Neighbour radius ignored");
  options.radiusKm = 100.; options.minimumPhaseDifference = 90.;
  check(!point(AutoMask::analyse(survey, {"A"}, options)).flagged[2], "Minimum angular departure ignored");
  options.minimumPhaseDifference = 15.;
  check(bytes(survey) == original, "Analysis changed values, masks or errors");
  auto cancelled = AutoMask::analyse(survey, {"A", "B"}, options, [](unsigned done, unsigned) { return done == 0; });
  check(cancelled.cancelled && cancelled.points.empty() && bytes(survey) == original, "Cancellation retained candidates or mutated data");
  survey.get_station_data("A").set_data_mask(RealZxy, 1./8., false);
  check(!std::isfinite(point(AutoMask::analyse(survey, {"A"}, options)).scores[1]), "Masked reference gap bridged");
  survey.get_station_data("A").set_data_mask(RealZxy, 1./16., false);
  result = AutoMask::analyse(survey, {"A"}, options);
  for(const auto &p : result.points) check(p.component != RealZxy || p.frequency != 1./16., "Masked target reassessed");
  survey = fixture(directory, false);
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(!p.candidate && !p.flagged[1] && !p.flagged[2], "Smooth trend flagged");
  for(const auto &name : survey.get_stations_names()) {
    auto &s = survey.get_station_data(name);
    for(double f : s.frequencies()) s.set_data(f, {PhsZxy, RhoZxy}, {name == "A" ? -179. : 179., name == "A" ? 10000. : 100.}, {1., 1.});
  }
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(!p.flagged[1] && !p.flagged[2], "Circular phase seam or static amplitude shift flagged");
  for(const auto *name : {"B", "C", "D"}) survey.set_active_flag(name, false);
  check(!std::isfinite(point(AutoMask::analyse(survey, {"A"}, options)).scores[2]), "Disabled neighbours used");
  survey.set_active_flag("A", false); result = AutoMask::analyse(survey, {"A"}, options);
  check(result.disabled == 1 && result.points.empty(), "Disabled target analysed");
}
void impedance(const QString &directory) {
  auto survey = fixture(directory); auto &station = survey.get_station_data("A");
  const double pi = std::acos(-1.), mu = 4e-7*pi;
  station.set_size(24, true); std::set<double> frequencies;
  for(unsigned i = 0; i < 24; ++i) frequencies.insert(1./std::pow(10., -3.+i*7./23.));
  station.set_frequencies(frequencies);
  for(double f : frequencies) {
    const auto top = std::sqrt(std::complex<double>(0., 2*pi*f*mu*10.)), bottom = std::sqrt(std::complex<double>(0., 2*pi*f*mu*1000.));
    const auto k = std::tanh(std::sqrt(std::complex<double>(0., 2*pi*f*mu/10.))*1500.);
    const auto z = top*(bottom+top*k)/(top+bottom*k); const double e = .001*std::abs(z);
    station.set_data(f, {RealZxy, ImagZxy}, {z.real(), z.imag()}, {e, e});
  }
  station.calculate_apparent_resistivity(); station.calculate_phase(); station.set_error_floor(0.);
  const auto original = bytes(survey);
  AutoMask::Options options; options.criteria = {{true, false, false}}; options.components = {{true, false, false, false}}; options.minimumFlags = 1;
  options.rhoPlus.intervalsPerDecade = 4;
  auto result = AutoMask::analyse(survey, {"A"}, options);
  for(const auto &p : result.points) check(!std::isfinite(p.scores[0]) && !p.candidate, "Provisional rho+ fit contributed a vote");
  options.rhoPlus.intervalsPerDecade = 12;
  const auto direct = ImpedanceAnalysis::rhoPlus(ImpedanceAnalysis::curve(station, ImpedanceAnalysis::Quantity::XY), options.rhoPlus);
  check(direct.converged && direct.gridStable, "Stable rho+ fixture failed");
  result = AutoMask::analyse(survey, {"A"}, options);
  check(result.points.size() == direct.observed.size(), "Rho+ sample count differs");
  for(unsigned i = 0; i < direct.observed.size(); ++i) {
    const double expected = std::max(std::abs(direct.rhoResidual[i]), std::abs(direct.phaseResidual[i]));
    check(std::abs(result.points[i].scores[0]-expected) < 1e-10, "Auto mask did not use main-data propagated errors");
  }
  check(bytes(survey) == original, "Rho+ analysis changed errors or masks");
}
void noisyTail(const QString &directory) {
  auto survey = fixture(directory, false); auto &station = survey.get_station_data("A");
  const std::array<double, 3> frequencies{{1./64., 1./128., 1./256.}}, rho{{.25, 30., 250.}}, phase{{5., 30., 90.}};
  for(unsigned i = 0; i < frequencies.size(); ++i)
    station.set_data(frequencies[i], {RhoZxy, PhsZxy}, {rho[i], phase[i]}, {2.*rho[i], 80.});
  const auto original = bytes(survey);
  AutoMask::Options options;
  auto result = AutoMask::analyse(survey, {"A"}, options);
  for(const auto &p : result.points) {
    const bool noisy = p.component == RealZxy && p.frequency <= frequencies.front();
    check(p.flagged[3] == noisy && (!noisy || p.candidate), "Low-SNR tail or endpoints escaped uncertainty screening");
    if(noisy) check(p.votes() == 1, "Fixture no longer reproduces noise hidden by normalized residuals");
  }
  std::exception_ptr dialogError;
  QTimer::singleShot(0, [&] {
    auto *d = qobject_cast<QDialog *>(QApplication::activeModalWidget());
    try {
      check(d, "Uncertainty review did not open");
      for(unsigned k = 0; k < 3; ++k)
        control<QCheckBox>(*d, ("autoMaskCriterion" + QString::number(k)).toUtf8().constData())->setChecked(false);
      auto *analyse = control<QPushButton>(*d, "autoMaskAnalyze"), *apply = control<QPushButton>(*d, "autoMaskApply");
      auto *table = control<QTableWidget>(*d, "autoMaskCandidates");
      analyse->click(); check(table->rowCount() == 3 && apply->isEnabled(), "Uncertainty-only review omitted noisy endpoints");
      for(int i = 0; i < table->rowCount(); ++i)
        check(table->item(i, 0)->checkState() == Qt::Checked && table->item(i, 7)->text().toDouble() > 1., "Uncertainty suggestions or ratio missing");
      auto *criterion = control<QCheckBox>(*d, "autoMaskCriterion3"); criterion->setChecked(false);
      check(!analyse->isEnabled() && !apply->isEnabled() && table->rowCount() == 0 && !control<QDoubleSpinBox>(*d, "autoMaskRhoError")->isEnabled(), "Disabled uncertainty left stale suggestions or active settings");
      criterion->setChecked(true);
      control<QDoubleSpinBox>(*d, "autoMaskRhoError")->setValue(300.);
      control<QDoubleSpinBox>(*d, "autoMaskPhaseError")->setValue(90.);
      analyse->click(); check(!apply->isEnabled() && table->rowCount() == 0, "UI ignored custom uncertainty limits");
      control<QDoubleSpinBox>(*d, "autoMaskRhoError")->setValue(50.);
      control<QDoubleSpinBox>(*d, "autoMaskPhaseError")->setValue(15.);
      analyse->click(); check(table->rowCount() == 3 && apply->isEnabled(), "Restored uncertainty limits failed");
      const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
      if(!artifacts.isEmpty()) { QApplication::processEvents(); d->grab().save(artifacts + "/auto-mask-uncertainty.png"); }
      check(bytes(survey) == original, "Uncertainty preview changed errors or masks");
    } catch(...) { dialogError = std::current_exception(); }
    if(d) d->reject();
  });
  check(!reviewAutoMasks(nullptr, survey, {"A"}, false), "Uncertainty preview applied on cancel");
  if(dialogError) std::rethrow_exception(dialogError);
  options.criteria[3] = false;
  for(double f : frequencies) check(!point(AutoMask::analyse(survey, {"A"}, options), f).candidate, "Disabled uncertainty check still nominated tail");
  options.criteria = {{false, false, false, true}}; options.minimumFlags = 2;
  for(double f : frequencies) check(!point(AutoMask::analyse(survey, {"A"}, options), f).candidate, "Uncertainty bypassed requested agreement count");
  options.minimumFlags = 1; options.threshold[3] = 10.;
  for(double f : frequencies) check(!point(AutoMask::analyse(survey, {"A"}, options), f).candidate, "Uncertainty multiplier ignored");
  options.threshold[3] = 1.; options.maximumRelativeRhoError = 3.; options.maximumPhaseError = 90.;
  for(double f : frequencies) check(!point(AutoMask::analyse(survey, {"A"}, options), f).candidate, "Custom uncertainty limits ignored");
  options.maximumRelativeRhoError = .5;
  check(point(AutoMask::analyse(survey, {"A"}, options), frequencies.back()).candidate, "Rho uncertainty alone did not flag");
  options.maximumRelativeRhoError = 3.; options.maximumPhaseError = 15.;
  check(point(AutoMask::analyse(survey, {"A"}, options), frequencies.back()).candidate, "Phase uncertainty alone did not flag");
  check(bytes(survey) == original, "Uncertainty screening rewrote input errors");
  // A missing error is unavailable, not zero and not evidence of poor precision.
  station.set_data(frequencies.back(), {RhoZxy, PhsZxy}, {100., 45.}, {NAN, NAN});
  check(!std::isfinite(point(AutoMask::analyse(survey, {"A"}, options), frequencies.back()).scores[3]), "Unknown errors treated as measured precision");
  station.set_data_mask(RealZxy, frequencies.front(), false);
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points)
    check(p.component != RealZxy || p.frequency != frequencies.front(), "Already masked noisy point reselected");
}
int pending(QDialog &d, unsigned pane) {
  auto *plot = control<QCustomPlot>(d, ("groupMaskPlot" + QString::number(pane)).toUtf8().constData()); int count = 0;
  for(int i = 0; i < plot->graphCount(); ++i) if(plot->graph(i)->property("maskState").toString() == "pending")
    for(auto it = plot->graph(i)->data()->constBegin(); it != plot->graph(i)->data()->constEnd(); ++it) count += std::isfinite(it->value);
  return count;
}
void tippers(const QString &directory) {
  auto survey = fixture(directory, false); auto &station = survey.get_station_data("A");
  station.set_data(1./16., {RealTzx, ImagTzx}, {-.8, .7}, {.01, .01});
  AutoMask::Options options; options.components = {{false, false, true, true}};
  auto result = AutoMask::analyse(survey, {"A"}, options);
  const auto &spike = point(result, 1./16., RealTzx);
  check(spike.candidate && spike.flagged[1] && spike.flagged[2] && !spike.flagged[0], "Tipper trend or complex neighbour check missed spike");
  for(const auto &p : result.points) check(!std::isfinite(p.scores[0]), "Rho+ evaluated tipper data");
  station.set_data(1./256., {RealTzy, ImagTzy}, {0., 0.}, {.2, .2});
  result = AutoMask::analyse(survey, {"A"}, options);
  check(point(result, 1./256., RealTzy).flagged[3], "Noisy zero tipper at endpoint missed");
  const auto original = survey;
  std::vector<bool> checked; for(const auto &p : result.points) checked.push_back(p.candidate);
  auto edits = AutoMask::maskPlan(survey, result.points, checked, true, true);
  AutoMask::apply(survey, edits);
  check(station.impedance_mask() == original.get_station_data("A").impedance_mask() &&
        station.phase_tensor_mask() == original.get_station_data("A").phase_tensor_mask(), "Tipper masking propagated into Z/PT");
  double v, e;
  const auto &f = station.frequencies(); const auto i = std::find(f.begin(), f.end(), 1./16.)-f.begin();
  check(!station.scalar_value(RealTzx, i, v, e) && !station.scalar_value(ImagTzx, i, v, e) && station.scalar_value(RealTzy, i, v, e), "Tipper real/imaginary masks or component independence wrong");
  survey = original;
  for(const auto &name : survey.get_stations_names()) for(double frequency : survey.get_station_data(name).frequencies())
    survey.get_station_data(name).set_data(frequency, {RealTzx, ImagTzx, RealTzy, ImagTzy}, {0., 0., 0., 0.}, {.01, .01, .01, .01});
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(!p.candidate, "Well-resolved zero tipper flagged");
  options.criteria = {{false, false, false, true}};
  survey.get_station_data("A").set_data(1./16., {RealTzx, ImagTzx}, {0., 0.}, {.2, .2});
  options.components[2] = false;
  for(const auto &p : AutoMask::analyse(survey, {"A"}, options).points) check(p.component != RealTzx, "Disabled tipper component analysed");
}
void settle() {
  QElapsedTimer clock; clock.start();
  while(clock.elapsed() < 230) { QApplication::processEvents(); QThread::msleep(5); }
}
void workflow(const QString &directory) {
  auto survey = std::make_shared<MTSurveyData>(fixture(directory, false));
  for(const auto *name : {"A", "C"}) survey->get_station_data(name).set_data(1./16., {RealTzx, ImagTzx}, {.1, .2}, {.2, .2});
  survey->get_station_data("A").set_data(1./32., {PhsZxy}, {50.}, {40.});
  survey->get_station_data("B").set_data(1./16., {RealTzy, ImagTzy}, {.3, .4}, {.2, .2});
  const auto path = directory + "/workflow.mtd";
  { std::ofstream file(path.toStdString(), std::ios::binary); boost::archive::binary_oarchive archive(file); archive << survey << std::map<std::string, MTSurveyData>{}; }
  MainWindow window;
  auto action = [&](const char *name) { return control<QAction>(window, name); };
  QTimer::singleShot(0, [&] {
    auto *dialog = qobject_cast<QFileDialog *>(QApplication::activeModalWidget());
    if(dialog) { dialog->selectFile(path); QMetaObject::invokeMethod(dialog, "accept", Qt::DirectConnection); }
  });
  action("actionLoad_project")->trigger();
  auto *list = control<QListWidget>(window, "stationList"); check(list->count() == 6, "Workflow project did not load");
  list->setCurrentRow(0); list->clearSelection(); list->item(0)->setSelected(true); list->item(2)->setSelected(true);
  auto *rho = control<QCustomPlot>(window, "plot11"), *phase = control<QCustomPlot>(window, "plot12"), *tip = control<QCustomPlot>(window, "plot21"), *pt = control<QCustomPlot>(window, "plot22");
  auto markers = [](QCustomPlot *plot) { int n = 0; for(int i = 0; i < plot->itemCount(); ++i) n += plot->item(i)->objectName() == "autoMaskMarker"; return n; };
  for(const auto *name : {"actionAuto_mask_check0", "actionAuto_mask_check1", "actionAuto_mask_check2"}) action(name)->setChecked(false);
  action("actionAuto_mask_preview")->setChecked(true); settle();
  check(markers(rho) == 2 && markers(phase) == 2 && markers(tip) == 2 && markers(pt) == 4, "Main plots omitted candidates or linked markers");
  check(!action("actionUndo")->isEnabled(), "Preview wrote history");
  rho->xAxis->setRange(2., 100.); rho->yAxis->setRange(1., 1000.);
  action("actionAuto_mask_Tzx")->setChecked(false); settle(); check(markers(tip) == 0 && markers(rho) == 2, "Component switch left tipper preview");
  check(rho->xAxis->range() == QCPRange(2., 100.) && rho->yAxis->range() == QCPRange(1., 1000.), "Auto-mask refresh reset zoom");
  action("actionAuto_mask_Tzx")->setChecked(true); settle();
  action("actionTipper_arrows")->setChecked(true); settle(); check(markers(tip) == 2, "Tipper arrow candidates missing");
  action("actionTipper_arrows")->setChecked(false); settle();
  action("actionAuto_mask_preview")->setChecked(false); check(markers(rho) == 0 && markers(tip) == 0, "Disabling preview left markers");
  action("actionAuto_mask_preview")->setChecked(true); settle();
  action("actionAuto_mask_apply_next")->trigger();
  check(list->currentItem()->text() == "B" && list->selectedItems().size() == 2 && list->item(0)->isSelected() && list->item(2)->isSelected(), "Apply/next changed batch selection or failed navigation");
  check(!action("actionAuto_mask_apply")->isEnabled(), "Navigation retained stale Apply"); settle();
  check(control<QLabel>(window, "autoMaskStatus")->text().contains("B:"), "Current station not identified while preserving batch selection");
  check(tip->graph(5)->dataCount() == 0, "Apply current also masked next station");
  action("actionAuto_mask_previous")->trigger(); settle();
  check(rho->graph(5)->dataCount() == 1 && rho->graph(4)->dataCount() == 1 && tip->graph(4)->dataCount() == 1 && tip->graph(6)->dataCount() == 1, "Apply missed previewed Z or tipper masks");
  action("actionUndo")->trigger(); settle();
  check(!action("actionUndo")->isEnabled() && rho->graph(5)->dataCount() == 0 && tip->graph(4)->dataCount() == 0 && markers(tip) == 2, "Undo did not restore masks and preview");
  QTimer::singleShot(0, [&] { if(auto *p = window.findChild<QProgressDialog *>("autoMaskBatchProgress")) p->cancel(); });
  action("actionAuto_mask_selected")->trigger(); settle();
  check(!action("actionUndo")->isEnabled() && tip->graph(4)->dataCount() == 0, "Cancelled batch applied changes");
  action("actionAuto_mask_selected")->trigger(); settle();
  check(tip->graph(4)->dataCount() == 1 && rho->graph(5)->dataCount() == 1, "Batch missed current selected station");
  action("actionAuto_mask_next")->trigger(); settle(); check(tip->graph(5)->dataCount() == 0, "Batch changed unselected station");
  action("actionAuto_mask_next")->trigger(); settle(); check(tip->graph(4)->dataCount() == 1 && rho->graph(5)->dataCount() == 0, "Batch missed second selection or propagated tipper to Z");
  check(list->item(0)->isSelected() && list->item(2)->isSelected(), "Batch navigation lost selection");
  action("actionUndo")->trigger(); settle(); check(!action("actionUndo")->isEnabled() && tip->graph(4)->dataCount() == 0, "Selected-site batch was not one Undo step");
  action("actionRedo")->trigger(); settle(); check(tip->graph(4)->dataCount() == 1, "Batch redo failed");
  list->clearSelection(); check(!action("actionAuto_mask_selected")->isEnabled(), "Batch enabled for empty selection");
  list->setCurrentRow(5, QItemSelectionModel::NoUpdate); check(!action("actionAuto_mask_next")->isEnabled(), "Next enabled beyond last station");
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if(!artifacts.isEmpty()) {
    list->setCurrentRow(1, QItemSelectionModel::NoUpdate); settle();
    action("actionClassic_menus")->setChecked(false); control<QTabBar>(window, "ribbonTabs")->setCurrentIndex(2);
    window.showNormal(); window.resize(1500, 900); settle(); window.grab().save(artifacts + "/auto-mask-workflow.png");
  }
}
void dialog(const QString &directory) {
  auto survey = std::make_shared<MTSurveyData>(fixture(directory)); const auto original = bytes(*survey);
  SurveyHistory history; history.reset(survey);
  for(unsigned mode = 0; mode < 3; ++mode) {
    std::exception_ptr error;
    QTimer::singleShot(0, [&] {
      auto *d = qobject_cast<QDialog *>(QApplication::activeModalWidget());
      try {
        check(d && d->objectName() == "autoMaskDialog", "Auto mask dialog missing");
        auto *apply = control<QPushButton>(*d, "autoMaskApply"); auto *analyse = control<QPushButton>(*d, "autoMaskAnalyze");
        check(!apply->isEnabled(), "Unanalysed survey can be masked");
        check(control<QSpinBox>(*d, "autoMaskMinimumFlags")->value() == 1, "Default still requires multiple checks for noisy points");
        control<QSpinBox>(*d, "autoMaskMinimumFlags")->setValue(2);
        control<QCheckBox>(*d, "autoMaskCriterion3")->setChecked(false);
        auto *rho = control<QCheckBox>(*d, "autoMaskCriterion0"); rho->setChecked(false);
        check(!control<QDoubleSpinBox>(*d, "autoMaskThreshold0")->isEnabled(), "Disabled check retained editable threshold");
        analyse->click(); auto *table = control<QTableWidget>(*d, "autoMaskCandidates");
        check(apply->isEnabled() && table->rowCount() >= 1, "No candidates in dialog");
        check(pending(*d, 0) == 2 && pending(*d, 2) == 0 && pending(*d, 3) == 4, "Plot preview or linkage incorrect");
        check(bytes(*survey) == original, "Preview altered source data");
        table->item(0, 0)->setCheckState(Qt::Unchecked);
        check(!apply->isEnabled() && pending(*d, 0) == 0, "Unchecking suggestion left mask plan");
        table->item(0, 0)->setCheckState(Qt::Checked);
        control<QCheckBox>(*d, "autoMaskPair")->setChecked(false);
        control<QCheckBox>(*d, "autoMaskLink")->setChecked(false);
        check(pending(*d, 0) == 1 && pending(*d, 3) == 0, "Optional associated edits ignored");
        auto *period = control<QCheckBox>(*d, "autoMaskCriterion1"); period->setChecked(false);
        check(!apply->isEnabled() && table->rowCount() == 0, "Changed criteria left stale masks");
        analyse->click(); check(!apply->isEnabled(), "Minimum agreement silently lowered");
        control<QSpinBox>(*d, "autoMaskMinimumFlags")->setValue(1); analyse->click();
        check(apply->isEnabled(), "Single enabled check cannot nominate candidates");
        control<QCheckBox>(*d, "autoMaskFlaggedOnly")->setChecked(false);
        check(table->rowCount() == 36, "All assessed points including tippers missing");
        if(mode == 0) {
          const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
          if(!artifacts.isEmpty()) { QApplication::processEvents(); d->grab().save(artifacts + "/auto-mask.png"); }
          const QString csv = directory + "/audit.csv";
          QTimer::singleShot(0, [&] {
            auto *file = qobject_cast<QFileDialog *>(QApplication::activeModalWidget());
            if(file) { file->selectFile(csv); QMetaObject::invokeMethod(file, "accept", Qt::DirectConnection); }
          });
          control<QPushButton>(*d, "autoMaskExport")->click();
          std::ifstream in(csv.toStdString()); const std::string text((std::istreambuf_iterator<char>(in)), {});
          check(text.find("criterion_0_enabled=0") != std::string::npos && text.find("criterion_2_enabled=1") != std::string::npos && text.find(directory.toStdString()) == std::string::npos,
                "Audit missing criterion choices or leaked local path");
          d->reject();
        } else if(mode == 1) {
          survey->get_station_data("B").set_data_mask(RealZxy, 1., false);
          apply->click(); check(!apply->isEnabled() && d->isVisible(), "Stale survey accepted mask plan"); d->reject();
          survey->get_station_data("B").set_data_mask(RealZxy, 1., true);
        } else apply->click();
      } catch(...) { error = std::current_exception(); if(d) d->reject(); }
    });
    const bool changed = reviewAutoMasks(nullptr, *survey, {"A"}, true);
    if(error) std::rethrow_exception(error);
    check(changed == (mode == 2), "Apply/cancel result wrong");
    if(mode != 2) check(bytes(*survey) == original, "Cancel changed source");
  }
  history.record("Auto masking");
  auto expected = fixture(directory); expected.get_station_data("A").set_data_mask(RealZxy, 1./16., false);
  check(bytes(*survey) == bytes(expected), "Applied unselected masks or changed measured errors");
  history.stack().undo(); check(bytes(*survey) == original, "Undo failed");
  history.stack().redo(); check(bytes(*survey) == bytes(expected), "Redo failed");
}
}
int main(int argc, char **argv) {
  QApplication app(argc, argv); QTemporaryDir directory;
  QCoreApplication::setOrganizationName("EDIToolsTests"); QCoreApplication::setApplicationName("AutoMask");
  QSettings::setDefaultFormat(QSettings::IniFormat); QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
  try { engine(directory.path()); impedance(directory.path()); noisyTail(directory.path()); tippers(directory.path()); dialog(directory.path()); workflow(directory.path()); }
  catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
