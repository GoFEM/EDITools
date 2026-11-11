#include "GroupMaskDialog.h"
#include "include/MTSurveyData.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QTableWidget>
#include <QTimer>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
#include <functional>
#include <sstream>

namespace {
void check(bool ok, const char *message) { if(!ok) throw std::runtime_error(message); }
template<class T> T *control(QDialog &dialog, const char *name) {
  auto *result = dialog.findChild<T *>(name); check(result, "Missing group-mask control"); return result;
}
std::string bytes(const MTSurveyData &survey) {
  std::ostringstream out; boost::archive::binary_oarchive archive(out); archive << survey; return out.str();
}
bool review(MTSurveyData &survey, bool linked, bool accept, const std::function<void(QDialog &)> &edit) {
  const auto before = bytes(survey);
  std::exception_ptr error;
  QTimer::singleShot(0, [&] {
    auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
    try {
      check(dialog && dialog->objectName() == "groupMaskDialog", "Group-mask dialog did not open");
      edit(*dialog);
      check(bytes(survey) == before, "Preview changed the source survey");
    } catch(...) { error = std::current_exception(); }
    if(dialog) { if(accept && !error) dialog->accept(); else dialog->reject(); }
  });
  const bool result = editGroupMasks(nullptr, survey, {"A", "B"}, linked);
  if(error) std::rethrow_exception(error);
  return result;
}
unsigned index(const MTSurveyData &survey, const char *name, double frequency) {
  const auto &f = survey.get_station_data(name).frequencies();
  return std::find(f.begin(), f.end(), frequency) - f.begin();
}
int plotted(QDialog &dialog, unsigned pane, const char *state) {
  auto *plot = control<QCustomPlot>(dialog, ("groupMaskPlot" + QString::number(pane)).toUtf8().constData());
  int count = 0;
  for(int i = 0; i < plot->graphCount(); ++i) if(plot->graph(i)->property("maskState").toString() == state)
    for(auto point = plot->graph(i)->data()->constBegin(); point != plot->graph(i)->data()->constEnd(); ++point)
      count += std::isfinite(point->value);
  return count;
}
}

void group_mask_tests(const QString &directory) {
  const auto path = directory + "/group-mask.gofem";
  const std::vector<RealDataType> types{RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy,
                                      RealTzx, ImagTzx, RealTzy, ImagTzy, PTxx, PTxy, PTyx, PTyy};
  {
    std::ofstream file(path.toStdString());
    for(const auto *name : {"A", "B", "C"})
      for(double f : std::string(name) == "B" ? std::vector<double>{1., .2} : std::vector<double>{1., 1.0005, .1})
        for(auto type : types) file << Datum::convert_type_to_string(type) << ' ' << f << " Plane_wave " << name << " .002 .0001\n";
  }
  MTSurveyData survey; survey.load_from_gofem(path.toStdString()); survey.set_error_floor(.1);
  survey.set_active_flag("B", false);
  const auto original = survey;
  auto zxy = [](QDialog &d) { control<QCheckBox>(d, "groupMaskZxy")->setChecked(true); };
  check(!review(survey, true, false, [&](QDialog &d) {
    check(!control<QPushButton>(d, "groupMaskApply")->isEnabled(), "No component selection enabled Apply");
    for(unsigned pane = 0; pane < 4; ++pane)
      check(plotted(d, pane, "observed") == 20 && plotted(d, pane, "pending") == 0, "Preview omitted station curves or included unselected stations");
    zxy(d);
    check(plotted(d, 0, "pending") == 10 && plotted(d, 1, "pending") == 10 && plotted(d, 2, "pending") == 0 && plotted(d, 3, "pending") == 20,
          "Plotted changes disagree with Z/PT masking plan");
    auto *plot = control<QCustomPlot>(d, "groupMaskPlot0");
    plot->xAxis->setRange(1., 5.); plot->yAxis->setRange(.01, 1.);
    control<QCheckBox>(d, "groupMaskVisible0_1")->setChecked(false);
    control<QCheckBox>(d, "groupMaskLink")->setChecked(false);
    check(plotted(d, 0, "pending") == 10 && plotted(d, 3, "pending") == 0, "Visibility changed masking or stale linked markers remained");
    check(plot->xAxis->range() == QCPRange(1., 5.) && plot->yAxis->range() == QCPRange(.01, 1.) &&
          control<QCustomPlot>(d, "groupMaskPlot2")->xAxis->range() == plot->xAxis->range(), "Preview lost zoom or period-axis synchronization");
    for(int i = 0; i < plot->graphCount(); ++i) if(plot->graph(i)->property("component").toUInt() == 1)
      check(!plot->graph(i)->visible(), "Component visibility omitted some station curves or highlights");
    control<QCheckBox>(d, "groupMaskLink")->setChecked(true);
    auto *table = control<QTableWidget>(d, "groupMaskPreview");
    check(table->rowCount() == 2 && table->item(0, 2)->text() == "6" && table->item(0, 3)->text() == "0" &&
          table->item(0, 4)->text() == "12" && table->item(1, 2)->text() == "4", "Linked preview counts or disabled-station scope wrong");
    const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
    if(!artifacts.isEmpty()) {
      control<QCheckBox>(d, "groupMaskVisible0_1")->setChecked(true);
      control<QPushButton>(d, "groupMaskFitView")->click();
      QApplication::processEvents(); d.grab().save(artifacts + "/group-mask-curves.png");
    }
  }) && bytes(survey) == bytes(original), "Cancel changed masks");
  check(review(survey, true, true, [&](QDialog &d) {
    zxy(d); control<QComboBox>(d, "groupMaskPeriods")->setCurrentIndex(1);
    control<QLineEdit>(d, "groupMaskLower")->setText("2"); control<QLineEdit>(d, "groupMaskUpper")->setText("1");
    check(!control<QPushButton>(d, "groupMaskApply")->isEnabled(), "Inverted range allowed");
    for(unsigned pane = 0; pane < 4; ++pane) check(plotted(d, pane, "pending") == 0, "Invalid range left stale highlight markers");
    control<QLineEdit>(d, "groupMaskLower")->setText("nan");
    check(!control<QPushButton>(d, "groupMaskApply")->isEnabled(), "Nonfinite range allowed");
    control<QLineEdit>(d, "groupMaskLower")->setText("1e0");
    check(plotted(d, 0, "pending") == 4 && plotted(d, 3, "pending") == 8, "Range markers included nearby periods");
  }), "Inclusive range did not apply");
  for(const auto *name : {"A", "B"}) {
    const auto &s = survey.get_station_data(name);
    for(unsigned i = 0; i < s.frequencies().size(); ++i) {
      const bool enabled = s.frequencies()[i] != 1.;
      check(s.impedance_mask(true)[1][i] == enabled && s.impedance_mask(true)[0][i] == enabled && s.impedance_mask(true)[2][i] && s.impedance_mask(true)[3][i], "Z selection missed its diagonal or changed the other pair");
      for(unsigned c = 0; c < 4; ++c) check(s.phase_tensor_mask(true)[c][i] == enabled, "Z/PT link missed a component or period");
      check(s.tipper_mask(true)[0][i] && s.tipper_mask(true)[1][i], "Z masking changed tippers");
    }
  }
  check(review(survey, true, true, [](QDialog &d) {
    control<QComboBox>(d, "groupMaskOperation")->setCurrentIndex(1);
    control<QCheckBox>(d, "groupMaskPTxy")->setChecked(true);
    check(plotted(d, 0, "masked") == 4 && plotted(d, 0, "observed") == 16 && plotted(d, 0, "pending") == 4 &&
          plotted(d, 3, "pending") == 2, "Unmask preview lost masked data or included unrelated PT entries");
    auto *plot = control<QCustomPlot>(d, "groupMaskPlot0");
    for(int i = 0; i < plot->graphCount(); ++i) if(plot->graph(i)->property("maskState") == "pending" && plot->graph(i)->dataCount())
      check(plot->graph(i)->scatterStyle().shape() == QCPScatterStyle::ssCircle && plot->graph(i)->pen().color() == QColor("#14804a"), "Unmask marker style indistinguishable from masking");
  }), "Linked unmask did not apply");
  const auto a1 = index(survey, "A", 1.);
  check(survey.get_station_data("A").impedance_mask(true)[1][a1] && survey.get_station_data("A").phase_tensor_mask(true)[1][a1] &&
        !survey.get_station_data("A").phase_tensor_mask(true)[0][a1], "PT unmask recursively changed other PT entries");
  check(review(survey, false, true, [](QDialog &d) {
    control<QCheckBox>(d, "groupMaskTzx")->setChecked(true);
    control<QComboBox>(d, "groupMaskPeriods")->setCurrentIndex(2);
    control<QLineEdit>(d, "groupMaskLower")->setText("1"); control<QLineEdit>(d, "groupMaskUpper")->setText("5");
  }), "Outside range did not apply");
  check(survey.get_station_data("A").tipper_mask(true)[0][a1] && !survey.get_station_data("A").tipper_mask(true)[0][index(survey, "A", .1)] &&
        survey.get_station_data("B").tipper_mask(true)[0][index(survey, "B", .2)], "Outside range included endpoints");
  check(review(survey, false, true, [](QDialog &d) {
    control<QCheckBox>(d, "groupMaskTzy")->setChecked(true);
    control<QComboBox>(d, "groupMaskPeriods")->setCurrentIndex(3);
    check(!control<QPushButton>(d, "groupMaskApply")->isEnabled(), "Empty period selection enabled Apply");
    auto *list = control<QListWidget>(d, "groupMaskPeriodList");
    for(int i = 0; i < list->count(); ++i) if(list->item(i)->data(Qt::UserRole).toDouble() == 1. / 1.0005) list->item(i)->setSelected(true);
    check(plotted(d, 2, "pending") == 2 && plotted(d, 0, "pending") == 0, "Tipper preview did not highlight both complex parts independently");
  }), "Individual period did not apply");
  check(!survey.get_station_data("A").tipper_mask(true)[1][index(survey, "A", 1.0005)] && survey.get_station_data("A").tipper_mask(true)[1][a1],
        "Exact period selection masked a nearby period");
  auto restoredMasks = survey;
  for(const auto *name : {"A", "B"}) for(auto type : types) restoredMasks.get_station_data(name).mask_type(type, true);
  check(bytes(restoredMasks) == bytes(original), "Group masks changed values, errors, activity or unselected stations");
  survey = original;
  auto all = [](QDialog &d) { control<QPushButton>(d, "groupMaskAll")->click(); };
  check(review(survey, true, true, all), "All components did not apply");
  check(!review(survey, true, true, [&](QDialog &d) {
    all(d); check(!control<QPushButton>(d, "groupMaskApply")->isEnabled(), "No-op masks enabled Apply");
  }), "Repeated mask reported changes");
  check(review(survey, true, true, [&](QDialog &d) { all(d); control<QComboBox>(d, "groupMaskOperation")->setCurrentIndex(1); }), "Unmask all failed");
  check(bytes(survey) == bytes(original), "Unmask all did not restore original flags");
  check(!review(survey, false, true, [&](QDialog &d) { all(d); control<QPushButton>(d, "groupMaskClear")->click(); }), "Clear components still applied changes");
  for(bool enabled: {false, true}) {
    check(review(survey, false, true, [&](QDialog &d) {
      control<QCheckBox>(d, "groupMaskZyx")->setChecked(true);
      control<QComboBox>(d, "groupMaskOperation")->setCurrentIndex(enabled ? 1 : 0);
      check(plotted(d, 0, "pending") == 10 && plotted(d, 3, "pending") == 0,
            "Unlinked group preview missed YX/YY pairing or changed PT");
    }), "Unlinked YX/YY group edit did not apply");
    for(const auto *name: {"A", "B"}) {
      const auto &s = survey.get_station_data(name);
      for(unsigned i = 0; i < s.frequencies().size(); ++i)
        check(s.impedance_mask(true)[0][i] && s.impedance_mask(true)[1][i] &&
              s.impedance_mask(true)[2][i] == enabled && s.impedance_mask(true)[3][i] == enabled,
              "Unlinked group operation edited the wrong pair");
      check(s.phase_tensor_mask(true) == original.get_station_data(name).phase_tensor_mask(true), "Unlinked group operation changed PT");
    }
  }
  check(bytes(survey) == bytes(original), "Paired group unmask did not restore the survey");
}
