#include "mainwindow.h"
#include "include/NativeMT.h"
#include "include/FitStatistics.h"
#include "PeriodMapWindow.h"
#include "PeriodLayoutWindow.h"
#include "include/SurveyReport.h"
#include <QApplication>
#include <QFileDialog>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QTreeWidget>
#include <QTabWidget>
#include <QSettings>
#include <QTemporaryDir>
#include <QTimer>
#include <QRegularExpression>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>

void period_resampling_tests(const QString &directory);

namespace {
void check(bool ok, const std::string &message)
{
  if(!ok) throw std::runtime_error(message);
}
template<class T> T *widget(QWidget &parent, const char *name)
{
  auto *result = parent.findChild<T *>(name);
  check(result != nullptr, std::string("Missing control: ") + name);
  return result;
}
void map_gesture(QCustomPlot *plot, QPointF start, QPointF end, bool drag = false,
                 Qt::KeyboardModifiers modifiers = Qt::NoModifier)
{
  QMouseEvent press(QEvent::MouseButtonPress, start, Qt::LeftButton, Qt::LeftButton, modifiers);
  QMouseEvent move(QEvent::MouseMove, end, Qt::NoButton, Qt::LeftButton, modifiers);
  QMouseEvent release(QEvent::MouseButtonRelease, end, Qt::LeftButton, Qt::NoButton, modifiers);
  QApplication::sendEvent(plot, &press);
  if(drag) QApplication::sendEvent(plot, &move);
  QApplication::sendEvent(plot, &release);
  QApplication::processEvents();
}
void file_action(QWidget &window, const char *name, const QStringList &paths)
{
  bool handled = false;
  QTimer timer;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    for(auto *top: QApplication::topLevelWidgets()) {
      auto *dialog = qobject_cast<QFileDialog *>(top);
      if(!dialog || !dialog->isVisible()) continue;
      dialog->setDirectory(QFileInfo(paths.front()).absolutePath());
      QStringList quoted;
      for(const auto &path: paths) quoted.push_back('"' + QFileInfo(path).fileName() + '"');
      auto *nameEdit = dialog->findChild<QLineEdit *>("fileNameEdit");
      if(!nameEdit) continue;
      nameEdit->setText(quoted.join(' '));
      handled = true;
      timer.stop();
      QMetaObject::invokeMethod(dialog, "accept", Qt::QueuedConnection);
      return;
    }
  });
  timer.start(10);
  if(auto *action = window.findChild<QAction *>(name)) action->trigger();
  else widget<QPushButton>(window, name)->click();
  timer.stop();
  check(handled, std::string("File dialog did not open: ") + name);
}
void write(const QString &path, const std::string &contents)
{
  std::ofstream output(path.toStdString());
  output << contents;
  check(bool(output), "Failed writing fixture");
}
QByteArray read_report(const QString &path)
{
  QFile file(path); check(file.open(QIODevice::ReadOnly), "Report missing"); return file.readAll();
}
void report_pages(const QString &path, unsigned expected)
{
  const auto pdf = read_report(path);
  check(pdf.startsWith("%PDF-"), "Report is not PDF");
  auto pages = QRegularExpression("/Type\\s*/Page\\b").globalMatch(QString::fromLatin1(pdf.constData(), pdf.size()));
  unsigned count = 0; while(pages.hasNext()) { pages.next(); ++count; }
  check(count == expected, "Report must have an overview and exactly one page per included station");
}
void check_report_axes(const std::array<SurveyReport::Options::Axis, 4> &axes)
{
  unsigned found = 0;
  for(auto *top: QApplication::topLevelWidgets()) {
    auto *plot = qobject_cast<QCustomPlot *>(top);
    if(!plot) continue;
    for(unsigned i = 0; i < axes.size(); ++i) if(plot->objectName() == QString("surveyReportPlot%1").arg(i)) {
      ++found;
      if(!axes[i].autoscale) check(plot->yAxis->range() == QCPRange(axes[i].lower, axes[i].upper), "PDF changed a fixed Y range");
      else check(plot->yAxis->range() != QCPRange(axes[i].lower, axes[i].upper), "PDF used stale fixed limits for an automatic Y axis");
    }
  }
  check(found == 4, "Could not inspect all rendered report panels");
}
void export_report(MainWindow &window, const QString &path,
                   const std::array<SurveyReport::Options::Axis, 4> *expectedAxes = nullptr)
{
  auto *stations = widget<QListWidget>(window, "stationList");
  auto *responses = widget<QListWidget>(window, "responsesList");
  const int stationRow = stations->currentRow(), responseRow = responses->currentRow();
  QTimer inspectTimer; unsigned inspected = 0;
  if(expectedAxes) {
    QObject::connect(&inspectTimer, &QTimer::timeout, [&] {
      for(auto *top: QApplication::topLevelWidgets()) if(top->objectName() == "surveyReportPlot0") {
        check_report_axes(*expectedAxes); ++inspected; break;
      }
    });
    inspectTimer.start(0);
  }
  QTimer timer; bool handled = false;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    auto *dialog = window.findChild<QDialog *>("surveyReportDialog");
    if(!dialog || !dialog->isVisible()) return;
    check(widget<QComboBox>(*dialog, "reportResponse")->currentData().toString() ==
          (responses->currentItem() ? responses->currentItem()->toolTip() : QString()), "Report should select current response");
    check(widget<QToolButton>(*dialog, "reportHelp")->toolTip().contains("overview"), "Report help missing");
    handled = true; timer.stop(); dialog->accept();
  });
  // Exercise the default suffix when the remaining filename has no other dots.
  const auto enteredPath = QFileInfo(path).completeBaseName().contains('.') ? path : path.chopped(4);
  timer.start(10); file_action(window, "actionExport_survey_report", {enteredPath});
  check(handled, "Report options did not open");
  check(!expectedAxes || inspected > 0, "GUI fixed ranges were not checked during export");
  report_pages(path, stations->count() + 1);
  check(stations->currentRow() == stationRow && responses->currentRow() == responseRow, "Report changed current selection");
}
void main_axis_ranges(MainWindow &window, bool autoscale, const std::array<SurveyReport::Options::Axis, 4> &axes)
{
  QTimer timer; bool handled = false;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    for(auto *top: QApplication::topLevelWidgets()) {
      auto *dialog = qobject_cast<QDialog *>(top);
      if(!dialog || !dialog->isVisible() || dialog->windowTitle() != "Axis ranges") continue;
      const auto edits = dialog->findChildren<QLineEdit *>();
      check(edits.size() == 8, "Unexpected main axis controls");
      dialog->findChild<QCheckBox *>()->setChecked(autoscale);
      for(unsigned i = 0; i < axes.size(); ++i) {
        edits[2 * i]->setText(QString::number(axes[i].lower)); edits[2 * i + 1]->setText(QString::number(axes[i].upper));
      }
      handled = true; timer.stop(); dialog->accept(); return;
    }
  });
  timer.start(0); widget<QAction>(window, "actionPlot_axis_ranges")->trigger();
  check(handled, "Main axis dialog did not open");
}
void survey_reports(const QString &directory)
{
  // Period counts distinguish complete tensors, partial components, masks, and
  // missing values without counting real/imaginary parts as separate periods.
  const auto layoutInput = directory + "/report-periods.gofem";
  std::ostringstream layout;
  for(int f: {1, 16}) for(const auto *type: {"RealZxx", "ImagZxx", "RealZxy", "ImagZxy", "RealZyx", "ImagZyx", "RealZyy", "ImagZyy"})
    layout << type << ' ' << f << " Plane_wave LAYOUT 0 0\n";
  layout << "RealZxy 2 Plane_wave LAYOUT 1 .1\nImagZxy 2 Plane_wave LAYOUT 1 .1\n"
            "RealZyx 2 Plane_wave LAYOUT 1 .1\nImagZyx 2 Plane_wave LAYOUT 1 .1\n"
            "RealZxx 4 Plane_wave LAYOUT 1 .1\nRealTzx 8 Plane_wave LAYOUT 0 0\n";
  write(layoutInput, layout.str());
  MTSurveyData layoutSurvey; layoutSurvey.load_from_gofem(layoutInput.toStdString());
  auto &layoutStation = layoutSurvey.get_station_data("LAYOUT");
  layoutStation.set_data_mask(RealZxx, 4., false); layoutStation.set_data_mask(RealZxx, 16., false);
  const auto summary = SurveyReport::summarize_periods(layoutStation);
  check(summary.counts[0] == std::array<unsigned, 4>{{1, 2, 1, 1}}, "Wrong full/partial/masked/missing impedance period counts");
  check(summary.counts[1] == std::array<unsigned, 4>{{0, 1, 0, 4}}, "A lone finite tipper scalar must count as one partial period");
  check(summary.periods.size() == 5 && summary.periods.at(1.)[0] == SurveyReport::PeriodState::Full &&
        summary.periods.at(.25)[0] == SurveyReport::PeriodState::Masked, "Coverage must retain exact periods and all-masked entries");
  for(const auto &counts: summary.counts) check(counts[0] + counts[1] + counts[2] + counts[3] == 5, "Period categories must be exhaustive and disjoint");
  layoutStation.set_active(false);
  const auto disabledSummary = SurveyReport::summarize_periods(layoutStation);
  check(disabledSummary.counts[0] == std::array<unsigned, 4>{{0, 0, 4, 1}}, "Disabled periods must remain distinct from absent data");
  const auto input = directory + "/report-input.gofem";
  write(input, "RealZxy 1 Plane_wave MAP .002 .0001\nImagZxy 1 Plane_wave MAP .002 .0001\n"
               "RealZxy 2 Plane_wave MAP .003 .0001\nImagZxy 2 Plane_wave MAP .003 .0001\n"
               "RealTzx 1 Plane_wave NO_LOCATION .1 .01\nImagTzx 1 Plane_wave NO_LOCATION .2 .01\n");
  MTSurveyData survey; survey.load_from_gofem(input.toStdString());
  survey.get_station_data("MAP").set_position({{55., 9., 100.}});
  survey.get_station_data("MAP").set_data_mask(RealZxy, 2., false);
  survey.get_station_data("NO_LOCATION").set_position({{NAN, NAN, NAN}});
  survey.set_active_flag("NO_LOCATION", false);
  const auto masks = survey.get_station_data("MAP").impedance_mask();
  SurveyReport::Options options;
  std::vector<unsigned> progress;
  const auto all = directory + "/report-all.pdf";
  check(SurveyReport::write_pdf(all, survey, nullptr, options, [&](unsigned done, unsigned total) {
    check(total == 3, "Report progress has wrong page total"); progress.push_back(done); return true;
  }), "Report export failed");
  report_pages(all, 3); check(progress == std::vector<unsigned>({0, 1, 2, 3}), "Report progress missing pages");
  options.includeDisabled = false; options.responseName = "GoFEM response";
  options.errorBars = false; options.components[0].fill(false);
  MTSurveyData response; response.load_from_gofem(input.toStdString());
  const auto enabled = directory + "/report-enabled.pdf";
  check(SurveyReport::write_pdf(enabled, survey, &response, options), "Response report export failed"); report_pages(enabled, 2);
  const auto original = read_report(enabled);
  auto fixedOptions = options; fixedOptions.includeDisabled = true;
  fixedOptions.axes = {{{false, 10., 100.}, {true, 1234., 5678.}, {false, -.3, .3}, {false, -5., 5.}}};
  unsigned checkedPages = 0;
  check(SurveyReport::write_pdf(directory + "/report-ranges.pdf", survey, &response, fixedOptions, [&](unsigned done, unsigned) {
    if(done >= 2) { check_report_axes(fixedOptions.axes); ++checkedPages; }
    return true;
  }), "Fixed range export failed");
  check(checkedPages == 2, "Fixed ranges must persist across all station pages, including empty panels");
  check(!SurveyReport::write_pdf(enabled, survey, nullptr, options, [](unsigned done, unsigned) { return done < 1; }), "Report cancellation ignored");
  check(read_report(enabled) == original, "Cancellation replaced existing report");
  check(!survey.is_active("NO_LOCATION") && survey.get_station_data("MAP").impedance_mask() == masks, "Report changed source masks");
  survey.set_active_flag("MAP", false);
  bool rejected = false;
  try { SurveyReport::write_pdf(enabled, survey, nullptr, options); } catch(const std::runtime_error &) { rejected = true; }
  check(rejected && read_report(enabled) == original, "Empty selection must preserve existing destination");
}
std::string gofem_text(const std::vector<NativeMT::Observation> &rows)
{
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(17) << "  # GoFEM response fixture\n\n";
  for(const auto &row: rows)
    text << Datum::convert_type_to_string(row.type) << ' ' << row.frequency << " Plane_wave "
         << row.receiver << ' ' << row.value << ' ' << row.error << '\n';
  return text.str();
}
double value(QCustomPlot &plot, int graph, double period)
{
  const auto data = plot.graph(graph)->data();
  const auto it = data->findBegin(period, false);
  check(it != data->end() && it->key == period, "Missing response period");
  return it->value;
}
void near(double actual, double expected)
{
  check(std::isfinite(actual) && std::abs(actual - expected) < 1e-10, "Unexpected plotted value");
}
void select_station(MainWindow &window, const QString &name)
{
  auto *list = widget<QListWidget>(window, "stationList");
  const auto items = list->findItems(name, Qt::MatchExactly);
  check(items.size() == 1, "Station name did not match observations");
  list->setCurrentItem(items.front());
}

QCheckBox *component_box(QCustomPlot &plot, int component)
{
  return widget<QCheckBox>(plot, ("componentVisibility" + QString::number(component)).toUtf8().constData());
}

void set_visibility(QCustomPlot &plot, int mask)
{
  for(int i = 0; i < 4; ++i) {
    auto *box = component_box(plot, i);
    if(box->isChecked() != bool(mask & (1 << i))) box->click();
  }
}

void check_visibility(QCustomPlot &plot, int mask, bool errorBars = true)
{
  for(int i = 0; i < plot.graphCount(); ++i)
    check(plot.graph(i)->visible() == bool(mask & (1 << (i % 4))), "Component visibility is inconsistent");
  int errorCount = 0;
  for(int i = 0; i < plot.plottableCount(); ++i)
    if(auto *bars = dynamic_cast<QCPErrorBars *>(plot.plottable(i))) {
      ++errorCount;
      check(bars->visible() == (errorBars && bars->dataPlottable()->visible()), "Hidden component has visible error bars");
    }
  check(errorCount == 8, "Error bars were lost");
  check(plot.legend->itemCount() == 4 && plot.findChildren<QCheckBox *>().size() == 4,
        "Each legend item must keep exactly one visibility checkbox");
  for(int i = 0; i < 4; ++i) {
    auto *box = component_box(plot, i);
    check(box->isChecked() == bool(mask & (1 << i)), "Legend checkbox is out of sync");
    check(!box->isHidden() && box->isEnabled() && plot.legend->outerRect().contains(box->geometry()),
          "Legend checkbox cannot be used to restore a hidden component");
  }
}

void range_dialog(QDialog &analysis, const std::function<void(QDialog &)> &edit)
{
  bool handled = false;
  std::exception_ptr failure;
  QTimer timer;
  QObject::connect(&timer, &QTimer::timeout, [&] {
    auto *dialog = analysis.findChild<QDialog *>("fitAxisRangesDialog");
    if(!dialog || !dialog->isVisible()) return;
    timer.stop();
    handled = true;
    try { edit(*dialog); }
    catch(...) { failure = std::current_exception(); dialog->reject(); }
  });
  timer.start(10);
  widget<QPushButton>(analysis, "fitPlotRanges")->click();
  timer.stop();
  if(failure) std::rethrow_exception(failure);
  check(handled, "Plot range dialog did not open");
}

void statistics_controls(QDialog &analysis)
{
  auto *periods = widget<QCustomPlot>(analysis, "fitPeriods");
  auto *checks = widget<QTreeWidget>(analysis, "fitResponses");
  const double originalRms = checks->topLevelItem(1)->data(1, Qt::UserRole).toDouble();
  const auto originalCount = checks->topLevelItem(1)->text(2);
  auto *map = widget<QCustomPlot>(analysis, "fitMapA");
  const int labels = map->itemCount(), points = map->graphCount();
  check(labels > 0, "RMS station labels are missing");
  widget<QCheckBox>(analysis, "fitStationNames")->setChecked(false);
  check(map->itemCount() == 0 && map->graphCount() == points, "Hiding RMS labels changed its station points");
  widget<QCheckBox>(analysis, "fitStationNames")->setChecked(true);
  check(map->itemCount() == labels, "RMS labels cannot be restored");
  auto set_axis = [](QDialog &dialog, const QString &prefix, double lower, double upper) {
    widget<QCheckBox>(dialog, (prefix + "Auto").toUtf8().constData())->setChecked(false);
    widget<QLineEdit>(dialog, (prefix + "Min").toUtf8().constData())->setText(QString::number(lower));
    widget<QLineEdit>(dialog, (prefix + "Max").toUtf8().constData())->setText(QString::number(upper));
  };
  range_dialog(analysis, [&](QDialog &dialog) {
    set_axis(dialog, "fitOverallY", 0., 4.);
    set_axis(dialog, "fitPeriodsX", 0., 10.);
    set_axis(dialog, "fitPeriodsY", .5, 3.);
    set_axis(dialog, "fitHistogramX", -2., 2.);
    set_axis(dialog, "fitComponentsY", 0., 4.);
    set_axis(dialog, "fitStationsY", 0., 5.);
    auto *buttons = dialog.findChild<QDialogButtonBox *>();
    buttons->button(QDialogButtonBox::Ok)->click();
    check(dialog.isVisible() && !widget<QLabel>(dialog, "fitAxisRangeError")->text().isEmpty(),
          "Invalid logarithmic limits were accepted");
    set_axis(dialog, "fitPeriodsX", .1, 10.);
    set_axis(dialog, "fitPeriodsY", 3., .5);
    buttons->button(QDialogButtonBox::Ok)->click();
    check(dialog.isVisible(), "Reversed axis limits were accepted");
    set_axis(dialog, "fitPeriodsY", .5, 3.);
    buttons->button(QDialogButtonBox::Ok)->click();
  });
  auto check_ranges = [&] {
    near(periods->xAxis->range().lower, .1); near(periods->xAxis->range().upper, 10.);
    near(periods->yAxis->range().lower, .5); near(periods->yAxis->range().upper, 3.);
    auto *histogram = widget<QCustomPlot>(analysis, "fitHistogram");
    near(histogram->xAxis->range().lower, -2.); near(histogram->xAxis->range().upper, 2.);
    near(widget<QCustomPlot>(analysis, "fitOverall")->yAxis->range().upper, 4.);
    near(widget<QCustomPlot>(analysis, "fitComponents")->yAxis->range().upper, 4.);
    near(widget<QCustomPlot>(analysis, "fitStations")->yAxis->range().upper, 5.);
  };
  check_ranges();
  widget<QComboBox>(analysis, "fitCurveColors")->setCurrentText("Thermal");
  QCPColorGradient curveGradient(QCPColorGradient::gpThermal);
  check(periods->graph(0)->pen().color() == QColor::fromRgb(curveGradient.color(0., QCPRange(0., 1.))),
        "Curve palette did not change");
  widget<QComboBox>(analysis, "fitColorMap")->setCurrentText("Grayscale");
  widget<QCheckBox>(analysis, "fitReverseColors")->setChecked(true);
  widget<QCheckBox>(analysis, "fitMapRangeAuto")->setChecked(false);
  widget<QDoubleSpinBox>(analysis, "fitMapRangeMin")->setValue(.5);
  widget<QDoubleSpinBox>(analysis, "fitMapRangeMax")->setValue(3.);
  widget<QCheckBox>(analysis, "fitHeatRangeAuto")->setChecked(false);
  widget<QDoubleSpinBox>(analysis, "fitHeatRangeMin")->setValue(.25);
  widget<QDoubleSpinBox>(analysis, "fitHeatRangeMax")->setValue(5.);
  auto check_colors = [&] {
    const auto expectedGradient = QCPColorGradient(QCPColorGradient::gpGrayscale).inverted();
    for(const auto *name: {"fitMapA", "fitMapB", "fitHeatmapA", "fitHeatmapB"}) {
      const bool heatmap = QString(name).contains("Heatmap");
      const auto expectedRange = heatmap ? QCPRange(.25, 5.) : QCPRange(.5, 3.);
      auto *plot = widget<QCustomPlot>(analysis, name);
      auto *scale = dynamic_cast<QCPColorScale *>(plot->plotLayout()->element(1, 1));
      check(scale && scale->dataRange() == expectedRange, "Color limits changed or differ between panels");
      check(scale->gradient() == expectedGradient, "Colormap or reversal was lost");
      if(heatmap) {
        auto *map = dynamic_cast<QCPColorMap *>(plot->plottable(0));
        check(map && map->gradient() == expectedGradient && map->dataRange() == expectedRange,
              "Heatmap did not follow its color scale");
      } else {
        check(plot->graph(0)->scatterStyle().brush().color() == QColor(Qt::white),
              "Station value below the color limit did not use the endpoint color");
      }
    }
  };
  // Use the exactly fitting response in both panels to check endpoint colors.
  widget<QComboBox>(analysis, "fitDetailB")->setCurrentIndex(0);
  check_colors();
  checks->topLevelItem(1)->setCheckState(0, Qt::Unchecked);
  check_ranges(); check_colors();
  checks->topLevelItem(1)->setCheckState(0, Qt::Checked);
  auto *normalization = widget<QComboBox>(analysis, "fitErrorSource");
  normalization->setCurrentIndex(1);
  check_ranges(); check_colors();
  normalization->setCurrentIndex(0);
  near(checks->topLevelItem(1)->data(1, Qt::UserRole).toDouble(), originalRms);
  check(checks->topLevelItem(1)->text(2) == originalCount, "Display controls changed the comparison counts");
  auto *histogram = widget<QCustomPlot>(analysis, "fitHistogram");
  double retained = 0.;
  for(const auto &point: *histogram->graph(1)->data()) retained += point.value;
  near(retained, originalCount.section(" / ", 0, 0).toDouble());
  range_dialog(analysis, [](QDialog &dialog) {
    widget<QPushButton>(dialog, "fitAxesAutoscale")->click();
    dialog.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->click();
  });
  near(periods->yAxis->range().lower, 0.);
  check(periods->yAxis->range().upper != 3., "Autoscale failed to clear fixed limits");
  widget<QCheckBox>(analysis, "fitMapRangeAuto")->setChecked(true);
  widget<QCheckBox>(analysis, "fitHeatRangeAuto")->setChecked(true);
  auto *scale = dynamic_cast<QCPColorScale *>(widget<QCustomPlot>(analysis, "fitMapA")->plotLayout()->element(1, 1));
  near(scale->dataRange().lower, 0.);
  check(!widget<QDoubleSpinBox>(analysis, "fitMapRangeMin")->isEnabled(), "Auto color limits remain editable");
}

void statistics(const QString &directory)
{
  const auto observedPath = directory + "/statistics-observed.txt";
  const auto responsePath = directory + "/statistics-response.txt";
  write(observedPath, "1 S01 impedance_xy real 1 .1\n1 S01 impedance_xy imag 2 .1\n");
  write(responsePath, "1 S01 impedance_xy real 2 .2\n1 S01 impedance_xy imag 4 .4\n"
                      "1 S99 impedance_xy real 5 .2\n2 S01 impedance_xy real 5 .2\n");
  MTSurveyData observed, response;
  observed.load_from_native_responses(observedPath.toStdString());
  response.load_from_native_responses(responsePath.toStdString());
  auto result = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  check(result.total.count == 2 && result.available == 4, "Fit comparison lost sparse coverage counts");
  near(result.total.rms(), 5.);
  check(result.residuals == std::vector<double>({5., 5.}), "Scalar response errors or residual sign changed");
  near(result.periods.at(1.).rms(), 5.);
  near(result.components.at("Zxy").rms(), 5.);
  near(result.stations.at("S01").rms(), 5.);
  near(result.station_periods.at("S01").at(1.).rms(), 5.);
  result = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Observed);
  near(result.total.rms(), std::sqrt(250.));
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., false);
  result = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  check(result.total.count == 0 && std::isnan(result.total.rms()), "Masked comparisons must be empty");
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., true);
  observed.set_active_flag("S01", false);
  check(FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response).total.count == 0,
        "Disabled stations entered the statistics");
  check(FitStatistics::component_name(RhoZxy) == "Rho xy" &&
        FitStatistics::component_name(PTyx) == "PTyx" &&
        FitStatistics::component_name(ImagTzy) == "Tzy", "Component grouping is incorrect");

  observed.set_active_flag("S01", true);
  write(responsePath, "1.0000001 S01 impedance_xy real .9 .2\n");
  response.load_from_native_responses(responsePath.toStdString());
  result = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  near(result.residuals.at(0), -.5);
  check(result.periods.count(1.) == 1, "Rounded frequency did not align with observed periods");
  write(observedPath, "1 S01 phase_xy value 45 2\n");
  write(responsePath, "1 S01 phase_xy value 405 2\n");
  observed.load_from_native_responses(observedPath.toStdString());
  response.load_from_native_responses(responsePath.toStdString());
  near(FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response).total.rms(), 0.);
}

void gofem_responses(const QString &directory)
{
  const auto path = directory + "/sparse.gofem", observedPath = directory + "/gofem-observed.txt";
  write(path, " \t#! response header\r\n #\r\n \t\r\n"
              " ImagZxy 1 Plane_wave S01 4 .4\r\n"
              "RealTzx 1 Plane_wave S01 0 .02 # valid zero\r\n"
              "RealZxy 4 Plane_wave S01 3 0\r\n"
              "RealZxy 1 Plane_wave S01 2 .2\r\n");
  write(observedPath, "1 S01 impedance_xy real 1 .1\n1 S01 impedance_xy imag 2 .1\n"
                      "1 S01 induction_x real .1 .01\n4 S01 impedance_xy real 1 .2\n");
  MTSurveyData response, observed;
  response.load_responses(path.toStdString());
  observed.load_from_native_responses(observedPath.toStdString());
  const auto &station = response.get_station_data("S01");
  check(station.frequencies() == dvector({1., 4.}) && response.response_observations().size() == 4,
        "GoFEM sparse keys or source rows were lost");
  near(response.response_observations()[0].error, .4);
  near(response.response_observations()[3].error, .2);
  double scalar, error;
  check(station.scalar_value(RealTzx, 0, scalar, error) && scalar == 0., "GoFEM discarded a valid zero");
  for(auto type: {ImagTzx, RealTzy, ImagTzy, RealZxx, ImagZxx, PTxx})
    check(!station.scalar_value(type, 0, scalar, error), "GoFEM fabricated a missing component");
  check(!station.scalar_value(ImagZxy, 1, scalar, error) && !station.scalar_value(RhoZxy, 1, scalar, error),
        "Incomplete GoFEM impedance produced a derived curve");
  MTMapData::PhaseTensor tensor;
  std::array<double, 2> vector;
  check(!MTMapData::phase_tensor(station, 0, tensor) && !MTMapData::induction_vector(station, 0, false, vector),
        "Incomplete GoFEM components produced map glyphs");
  auto fit = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  check(fit.available == 4 && fit.total.count == 3, "GoFEM statistics invented scalars or counted zero errors");
  near(fit.total.rms(), 5.);
  check(fit.residuals == std::vector<double>({5., -5., 5.}), "GoFEM merged real/imaginary response errors");
  fit = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Observed);
  check(fit.total.count == 4, "Zero-error GoFEM response cannot use observed errors");
  near(fit.total.rms(), std::sqrt(175.));
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., false);
  fit = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  check(fit.total.count == 1, "GoFEM statistics ignored observed masks");
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., true);

  // Both missing values and individual scalar errors must survive a project archive.
  std::stringstream saved;
  { boost::archive::binary_oarchive archive(saved); archive << response; }
  MTSurveyData restored;
  { boost::archive::binary_iarchive archive(saved); archive >> restored; }
  near(FitStatistics::compare(observed, restored, FitStatistics::ErrorSource::Response).total.rms(), 5.);
  check(restored.response_observations().size() == 4 &&
        !restored.get_station_data("S01").scalar_value(ImagTzx, 0, scalar, error), "Project lost sparse GoFEM data");

  // Reloads replace all prior rows; malformed input leaves the current data intact.
  write(path, "RealTzy 1e-8 Plane_wave station:03 0 .01\nRealTzy 1.005e-8 Plane_wave station:03 .1 .02\n");
  response.load_from_gofem(path.toStdString());
  check(response.n_stations() == 1 && !response.is_station_present("S01") && response.response_observations().size() == 2,
        "GoFEM reload retained old response data");
  check(response.get_station_data("station:03").scalar_value(RealTzy, 1, scalar, error), "GoFEM collapsed distinct small frequencies");
  near(scalar, .1);
  for(const auto *invalid: {"RealZxy 0 Plane_wave S01 1 .1", "RealZxy 1 Plane_wave S01 1 -.1",
                           "RealZxy 1 Plane_wave S01 nan .1", "RealZxy 1 Plane_wave S01 1 inf",
                           "RealZxy 1 Plane_wave S01 1", "RealZxy 1 Plane_wave S01 1 .1extra",
                           "Unknown 1 Plane_wave S01 1 .1", "log10RhoZxy 1 Plane_wave S01 1 .1",
                           "RealZxy 1 Plane_wave S01 1 .1\nRealZxy 1 Other_source S01 2 .2"}) {
    write(path, std::string("# invalid input\n") + invalid + '\n');
    bool rejected = false;
    try { response.load_from_gofem(path.toStdString()); }
    catch(const std::exception &e) { rejected = std::string(e.what()).find("GoFEM line ") == 0; }
    check(rejected && response.is_station_present("station:03") && response.response_observations().size() == 2,
          "Invalid GoFEM input was accepted or changed existing data");
  }
  write(path, " # empty file\n\n");
  bool rejected = false;
  try { response.load_from_gofem(path.toStdString()); }
  catch(const std::exception &e) { rejected = std::string(e.what()).find("Empty GoFEM") == 0; }
  check(rejected, "Empty GoFEM file was accepted");

  // Explicit rho/phase/tensor values across frequencies must survive derivation.
  // Every supported scalar must have identical native and GoFEM representations.
  std::vector<NativeMT::Observation> rows;
  for(double frequency: {1., 2.})
    for(const auto &mapping: NativeMT::mappings())
      rows.push_back({frequency, "S01", mapping.type, double(rows.size() + 1), .01 * (rows.size() + 1)});
  std::reverse(rows.begin(), rows.end());
  write(path, gofem_text(rows));
  response.load_from_gofem(path.toStdString());
  check(response.response_observations().size() == rows.size(), "GoFEM lost a scalar type");
  const auto expected = MTResponseData::stations(rows);
  for(const auto &row: rows) {
    const auto index = row.frequency == 1. ? 0 : 1;
    double expectedValue, expectedError;
    check(response.get_station_data("S01").scalar_value(row.type, index, scalar, error) &&
          expected.at("S01").scalar_value(row.type, index, expectedValue, expectedError), "Missing scalar after GoFEM import");
    near(scalar, row.value); near(scalar, expectedValue); near(error, expectedError);
  }
  const auto &preserved = response.response_observations();
  for(unsigned i = 0; i < rows.size(); ++i) {
    near(preserved[i].value, rows[i].value); near(preserved[i].error, rows[i].error);
  }
}

void period_maps(const QString &directory)
{
  MTMapData::PhaseTensor tensor;
  const double degrees = 180. / std::acos(-1.);
  check(MTMapData::phase_tensor({{2., 0., 0., 1.}}, tensor), "Diagonal tensor rejected");
  near(tensor.phiMax, std::atan(2.) * degrees); near(tensor.phiMin, 45.);
  near(tensor.axisRatio, .5); near(tensor.azimuth, 0.); near(tensor.skew, 0.);
  check(MTMapData::phase_tensor({{1.75, std::sqrt(3.) / 4., std::sqrt(3.) / 4., 1.25}}, tensor), "Rotated tensor rejected");
  near(tensor.azimuth, 30.); near(tensor.axisRatio, .5);
  check(MTMapData::phase_tensor({{2., 1., 0., 1.}}, tensor), "Non-symmetric tensor rejected");
  // The ellipse's major direction is the eigenvector of Phi * transpose(Phi).
  near(tensor.azimuth, .5 * std::atan2(2., 4.) * degrees);
  check(MTMapData::phase_tensor({{2., 0., 0., -1.}}, tensor), "Negative determinant tensor rejected");
  near(tensor.phiMin, -45.); near(tensor.axisRatio, .5);
  check(MTMapData::phase_tensor({{1., 0., 0., 1.}}, tensor), "Isotropic tensor rejected");
  near(tensor.axisRatio, 1.);
  check(!MTMapData::phase_tensor({{0., 0., 0., 0.}}, tensor) &&
        !MTMapData::phase_tensor({{1., NAN, 0., 1.}}, tensor), "Undefined tensor produced an ellipse");
  check(MTMapData::nearest_period({0., NAN, 1.0001, 1., .1}, 1., .001) == 3 &&
        MTMapData::nearest_period({1.0001}, 1., .001) == 0 &&
        MTMapData::nearest_period({1.0001}, 1., 0.) == -1 &&
        MTMapData::nearest_period({1.}, 10., .25) == -1, "Period matching exceeded the tolerance");

  const auto observedPath = directory + "/map-observed.txt", responsePath = directory + "/map-response.gofem";
  auto rows = [](const std::string &station, const std::string &frequency) {
    std::ostringstream text;
    for(const auto *row: {"impedance_xx real 0", "impedance_xx imag 0", "impedance_xy real 1", "impedance_xy imag 1",
                         "impedance_yx real -1", "impedance_yx imag -2", "impedance_yy real 0", "impedance_yy imag 0",
                         "induction_x real .2", "induction_x imag .05", "induction_y real .1", "induction_y imag -.1"})
      text << frequency << ' ' << station << ' ' << row << " .01\n";
    return text.str();
  };
  write(observedPath, rows("A", "1") + rows("B", "1") + rows("A", ".1") + "1 C induction_x real .2 .01\n");
  std::istringstream responseRows(rows("A", "1"));
  write(responsePath, gofem_text(NativeMT::read_observations(responseRows)));
  auto survey = std::make_shared<MTSurveyData>();
  survey->load_from_native_responses(observedPath.toStdString());
  survey->get_station_data("A").set_position({{55., 9., 100.}});
  survey->get_station_data("B").set_position({{55.01, 9., 100.}});
  survey->get_station_data("C").set_position({{55., 9.02, 100.}});
  std::map<std::string, MTSurveyData> responses;
  responses[responsePath.toStdString()].load_from_gofem(responsePath.toStdString());

  auto station = survey->get_station_data("A");
  const auto index = MTMapData::nearest_period(station.frequencies(), 1., 0.);
  std::array<double, 2> vector;
  check(MTMapData::phase_tensor(station, index, tensor), "Impedance-derived tensor rejected");
  near(tensor.axisRatio, .5); near(tensor.azimuth, 0.);
  station.set_data_mask(RealZxx, 1., false);
  check(!MTMapData::phase_tensor(station, index, tensor), "Masked impedance entered the tensor map");
  station.set_data_mask(RealZxx, 1., true);
  station.set_data_mask(PTxy, 1., false);
  check(!MTMapData::phase_tensor(station, index, tensor), "Masked tensor entered the map");
  station.set_data_mask(PTxy, 1., true);
  station.set_data(1., {RealZxx, RealZxy, RealZyx, RealZyy}, {1., 2., 2., 4.}, {.01, .01, .01, .01});
  check(!MTMapData::phase_tensor(station, index, tensor), "Singular real impedance produced an ellipse");
  station.set_data(1., {RealTzx, RealTzy}, {0., 0.}, {.01, .01});
  check(MTMapData::induction_vector(station, index, false, vector), "Valid zero vector rejected");
  near(vector[0], 0.); near(vector[1], 0.);
  station.set_data_mask(RealTzy, 1., false);
  check(!MTMapData::induction_vector(station, index, false, vector), "Incomplete vector was plotted");
  check(!MTMapData::induction_vector(survey->get_station_data("C"), 0, false, vector), "Missing east component treated as zero");
  MTStationData explicitTensor;
  explicitTensor.set_size(1, true); explicitTensor.set_frequencies({1.});
  explicitTensor.set_data(1., {PTxx, PTxy, PTyx, PTyy}, {2., 0., 0., 1.}, {.01, .01, .01, .01});
  check(MTMapData::phase_tensor(explicitTensor, 0, tensor), "Explicit phase tensor without impedance was rejected");

  int maskChanges = 0;
  PeriodMapWindow window(nullptr, [&] { ++maskChanges; });
  window.setData(survey, responses); window.show(); QApplication::processEvents();
  auto *plot = widget<QCustomPlot>(window, "periodMapPlot");
  auto *period = widget<QComboBox>(window, "mapPeriod");
  auto *dataset = widget<QComboBox>(window, "mapDataset");
  near(period->currentData().toDouble(), 1.);
  check(period->count() == 2 && dataset->count() == 2 && plot->property("phaseTensorCount").toUInt() == 2,
        "Default period map coverage is incorrect");
  auto *ellipse = widget<QCPCurve>(*plot, "phaseTensor_A");
  near(ellipse->property("axisRatio").toDouble(), .5);
  auto *arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  const QPointF center = arrow->start->coords(), delta = arrow->end->coords() - center;
  const double size = widget<QDoubleSpinBox>(window, "mapArrowSize")->value();
  near(delta.x(), -.1 * size); near(delta.y(), -.2 * size);
  const auto majorPoint = ellipse->data()->constBegin();
  near(majorPoint->key, center.x());
  near(majorPoint->value - center.y(), .5 * widget<QDoubleSpinBox>(window, "mapEllipseSize")->value());
  check(!plot->findChild<QCPItemLine *>("inductionReal_C"), "Incomplete station vector displayed");
  widget<QComboBox>(window, "mapArrowConvention")->setCurrentIndex(1);
  arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  near((arrow->end->coords() + delta - center).manhattanLength(), 0.);
  widget<QDoubleSpinBox>(window, "mapArrowSize")->setValue(size * 2.);
  arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  near((arrow->end->coords() + 2. * delta - center).manhattanLength(), 0.);
  widget<QCheckBox>(window, "mapImagArrows")->setChecked(true);
  arrow = widget<QCPItemLine>(*plot, "inductionImag_A");
  near(arrow->end->coords().x() - center.x(), -.1 * size * 2.);
  near(arrow->end->coords().y() - center.y(), .05 * size * 2.);
  widget<QCheckBox>(window, "mapRealArrows")->setChecked(false);
  check(!plot->findChild<QCPItemLine *>("inductionReal_A") && plot->findChild<QCPCurve *>("phaseTensor_A"),
        "Arrow visibility affected phase tensors");
  widget<QCheckBox>(window, "mapImagArrows")->setChecked(false);
  check(!plot->legend->visible(), "Tensor-only map retained an empty arrow legend");
  widget<QCheckBox>(window, "mapImagArrows")->setChecked(true);
  widget<QCheckBox>(window, "mapPhaseTensors")->setChecked(false);
  check(!plot->findChild<QCPCurve *>("phaseTensor_A") && plot->findChild<QCPItemLine *>("inductionImag_A"),
        "Tensor visibility affected imaginary arrows");
  widget<QCheckBox>(window, "mapPhaseTensors")->setChecked(true);
  check(!plot->findChild<QCPItemText *>("periodMapLabel_A"), "Station names should start hidden");
  widget<QCheckBox>(window, "mapStationNames")->setChecked(true);
  check(widget<QCPItemText>(*plot, "periodMapLabel_A")->text() == "A", "Station label toggle failed");
  widget<QComboBox>(window, "mapColorMap")->setCurrentText("Grayscale");
  widget<QCheckBox>(window, "mapReverseColors")->setChecked(true);
  widget<QDoubleSpinBox>(window, "mapColorMin")->setValue(10.);
  widget<QDoubleSpinBox>(window, "mapColorMax")->setValue(60.);
  auto *scale = dynamic_cast<QCPColorScale *>(plot->plotLayout()->element(1, 1));
  check(scale && scale->dataRange() == QCPRange(10., 60.) &&
        scale->gradient() == QCPColorGradient(QCPColorGradient::gpGrayscale).inverted(), "Map colors did not follow controls");
  check(widget<QCPCurve>(*plot, "phaseTensor_A")->brush().color() == QColor(Qt::black), "Color limits do not clip to palette endpoints");
  widget<QPushButton>(window, "mapNextPeriod")->click();
  near(period->currentData().toDouble(), 10.);
  check(plot->property("phaseTensorCount").toUInt() == 1, "Unavailable period was extrapolated");
  widget<QPushButton>(window, "mapPreviousPeriod")->click();

  // Selection works on glyphs and on center dots after masking. Only this period changes.
  auto &observedA = survey->get_station_data("A");
  const int otherIndex = MTMapData::nearest_period(observedA.frequencies(), 10., 0.);
  auto *mask = widget<QPushButton>(window, "mapMask");
  auto *unmask = widget<QPushButton>(window, "mapUnmask");
  auto *group = widget<QComboBox>(window, "mapMaskGroup");
  check(!mask->isEnabled(), "Mask enabled without selected stations");
  widget<QCheckBox>(window, "mapSelectStations")->setChecked(true);
  auto selectedCount = [&] { return widget<QCPGraph>(*plot, "mapSelectedStations")->data()->size(); };
  auto centerPixel = [&] { return QPointF(plot->xAxis->coordToPixel(center.x()), plot->yAxis->coordToPixel(center.y())); };
  auto arrowTip = widget<QCPItemLine>(*plot, "inductionImag_A")->end->pixelPosition();
  map_gesture(plot, arrowTip, arrowTip);
  check(selectedCount() == 1 && mask->isEnabled(), "Clicking an arrow did not select its station");
  mask->click();
  check(maskChanges == 1 && !observedA.tipper_mask()[0][index] && !observedA.tipper_mask()[3][index] &&
        observedA.tipper_mask()[0][otherIndex] && observedA.impedance_mask()[0][index] &&
        !plot->findChild<QCPItemLine *>("inductionImag_A") && plot->findChild<QCPCurve *>("phaseTensor_A"),
        "Tipper masking changed other periods/components or left masked arrows");
  check(FitStatistics::compare(*survey, responses.begin()->second, FitStatistics::ErrorSource::Response).total.count == 8,
        "Tipper masking did not exclude all four scalars from statistics");
  group->setCurrentIndex(1); mask->click();
  check(!observedA.phase_tensor_mask()[0][index] && observedA.impedance_mask()[0][index] &&
        !plot->findChild<QCPCurve *>("phaseTensor_A"), "Tensor-only mask changed impedance or retained the ellipse");
  // Round-trip the same station masks that projects serialize.
  std::stringstream savedMasks;
  { boost::archive::binary_oarchive archive(savedMasks); archive << *survey; }
  MTSurveyData restoredMasks;
  { boost::archive::binary_iarchive archive(savedMasks); archive >> restoredMasks; }
  check(!restoredMasks.get_station_data("A").phase_tensor_mask()[0][index] &&
        !restoredMasks.get_station_data("A").tipper_mask()[0][index], "Map masks were lost on serialization");
  widget<QPushButton>(window, "mapClearSelection")->click();
  map_gesture(plot, centerPixel(), centerPixel());
  check(selectedCount() == 1, "Masked station center could not be selected");
  unmask->click();
  group->setCurrentIndex(0); unmask->click();
  group->setCurrentIndex(2); mask->click();
  check(!observedA.impedance_mask()[3][index] && !observedA.phase_tensor_mask()[3][index] &&
        observedA.impedance_mask()[3][otherIndex] && observedA.tipper_mask()[0][index],
        "Impedance + tensor masking changed the wrong components or period");
  check(FitStatistics::compare(*survey, responses.begin()->second, FitStatistics::ErrorSource::Response).total.count == 4,
        "Impedance masking did not exclude eight response scalars");
  unmask->click();
  // Ctrl-click toggles; dragging selects several centers, and disabled stations are skipped.
  map_gesture(plot, centerPixel(), centerPixel(), false, Qt::ControlModifier);
  check(selectedCount() == 0, "Ctrl-click did not remove the selected station");
  auto *bArrow = widget<QCPItemLine>(*plot, "inductionImag_B");
  const auto bPixel = bArrow->start->pixelPosition();
  map_gesture(plot, centerPixel() - QPointF(8, 8), centerPixel() + QPointF(8, 8), true);
  check(selectedCount() == 1, "Box selection failed");
  map_gesture(plot, bPixel, bPixel, false, Qt::ControlModifier);
  check(selectedCount() == 2, "Ctrl-click did not add a station");
  survey->set_active_flag("B", false); window.setData(survey, responses);
  mask->click(); survey->set_active_flag("B", true);
  check(survey->get_station_data("B").impedance_mask()[0][0], "Map masking modified a disabled station");
  unmask->click();
  widget<QPushButton>(window, "mapNextPeriod")->click();
  group->setCurrentIndex(0); mask->click();
  check(!observedA.tipper_mask()[0][otherIndex] && observedA.tipper_mask()[0][index] &&
        survey->get_station_data("B").tipper_mask()[0][0], "Masking extrapolated into an unavailable period");
  unmask->click(); widget<QPushButton>(window, "mapPreviousPeriod")->click();
  widget<QPushButton>(window, "mapClearSelection")->click();
  // An ellipse boundary is also a selection target.
  ellipse = widget<QCPCurve>(*plot, "phaseTensor_A");
  const auto edge = ellipse->data()->constBegin();
  const QPointF edgePixel(plot->xAxis->coordToPixel(edge->key), plot->yAxis->coordToPixel(edge->value));
  map_gesture(plot, edgePixel, edgePixel);
  check(selectedCount() == 1, "Clicking an ellipse did not select its station");
  dataset->setCurrentIndex(1);
  check(period->count() == 1 && plot->property("phaseTensorCount").toUInt() == 1, "Response map did not update");
  group->setCurrentIndex(2); mask->click();
  check(!observedA.impedance_mask()[0][index] && plot->property("phaseTensorCount").toUInt() == 1 &&
        responses.begin()->second.get_station_data("A").impedance_mask()[0][0],
        "Masking from a computed map modified the response instead of the observations");
  unmask->click();
  widget<QCheckBox>(window, "mapSelectStations")->setChecked(false);
  check(plot->interactions().testFlag(QCP::iRangeDrag), "Turning selection off did not restore panning");
  survey->set_active_flag("A", false);
  window.setData(survey, responses);
  check(dataset->currentIndex() == 1 && plot->property("phaseTensorCount").toUInt() == 0,
        "Refresh ignored disabled station or changed dataset selection");
  survey->set_active_flag("A", true);
  // A station away from the central meridian uses geographic north, not grid north.
  survey->get_station_data("A").set_position({{55., 12., 100.}});
  window.setData(survey, responses);
  widget<QCheckBox>(window, "mapRealArrows")->setChecked(true);
  arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  const auto projection = SurveyCoordinates::suggested(survey->geographic_locations());
  const auto northPoints = projection.transform({{{55., 12., 100.}}, {{55.00001, 12., 100.}}});
  const QPointF north(northPoints[1][1] - northPoints[0][1], northPoints[1][0] - northPoints[0][0]);
  const auto unitNorth = north / std::hypot(north.x(), north.y());
  const auto expected = widget<QDoubleSpinBox>(window, "mapArrowSize")->value() *
    (.2 * unitNorth + .1 * QPointF(unitNorth.y(), -unitNorth.x()));
  near((arrow->end->coords() - arrow->start->coords() - expected).manhattanLength(), 0.);
  check(scale->dataRange() == QCPRange(10., 60.), "Refreshing maps reset color limits");
  window.resize(1300, 800); QApplication::processEvents(); plot->replot();
  near(plot->xAxis->range().size() / plot->axisRect()->width(), plot->yAxis->range().size() / plot->axisRect()->height());
  const auto xRange = plot->xAxis->range(), yRange = plot->yAxis->range();
  file_action(window, "mapSavePdf", {directory + "/period-map.pdf"});
  check(QFileInfo(directory + "/period-map.pdf").size() > 1000 &&
        plot->xAxis->range() == xRange && plot->yAxis->range() == yRange, "PDF export changed the map view");
  window.setData({}, {});
  check(plot->property("phaseTensorCount").toUInt() == 0 && selectedCount() == 0 && !mask->isEnabled(),
        "Empty survey retained old map data or selections");

  // Close frequencies must not cause the first tolerance match to be masked instead.
  write(observedPath, rows("A", "1") + rows("A", "1.0005"));
  auto closeSurvey = std::make_shared<MTSurveyData>();
  closeSurvey->load_from_native_responses(observedPath.toStdString());
  auto &closeStation = closeSurvey->get_station_data("A");
  closeStation.set_position({{55., 9., 100.}});
  window.setData(closeSurvey, {});
  period->setCurrentIndex(0);
  near(period->currentData().toDouble(), 1. / 1.0005);
  widget<QCheckBox>(window, "mapSelectStations")->setChecked(true);
  const auto closeCenter = widget<QCPItemLine>(*plot, "inductionReal_A")->start->pixelPosition();
  map_gesture(plot, closeCenter, closeCenter);
  group->setCurrentIndex(0); mask->click();
  const int exactIndex = MTMapData::nearest_period(closeStation.frequencies(), 1. / 1.0005, 0.);
  const int neighborIndex = MTMapData::nearest_period(closeStation.frequencies(), 1., 0.);
  check(!closeStation.tipper_mask()[0][exactIndex] && closeStation.tipper_mask()[0][neighborIndex],
        "Map masked a nearby frequency instead of the selected frequency");
}

void resampling_project_workflow(const QString &directory)
{
  auto survey = std::make_shared<MTSurveyData>();
  survey->load_from_gofem((directory + "/resampling-input.gofem").toStdString());
  for(const auto &name: survey->get_stations_names()) survey->get_station_data(name).set_position({{55., 9., 100.}});
  const auto original = directory + "/layout-original.mtd", saved = directory + "/layout-copy.mtd";
  {
    std::ofstream file(original.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(file); archive << survey << std::map<std::string, MTSurveyData>{};
  }
  MainWindow window;
  file_action(window, "actionLoad_project", {original}); select_station(window, "A");
  widget<QAction>(window, "actionPeriod_layout")->trigger();
  auto *layout = widget<QDialog>(window, "periodLayoutWindow");
  widget<QComboBox>(*layout, "layoutGridMode")->setCurrentIndex(2);
  widget<QLineEdit>(*layout, "layoutCustomPeriods")->setText("1 2 4");
  widget<QDoubleSpinBox>(*layout, "layoutMaximumRatio")->setValue(4.);
  widget<QPushButton>(*layout, "layoutPreview")->click();
  file_action(*layout, "layoutExportReport", {directory + "/resampling-audit.csv"});
  QFile audit(directory + "/resampling-audit.csv"); check(audit.open(QIODevice::ReadOnly), "Audit CSV missing");
  const auto text = audit.readAll();
  check(text.contains("\"A\",Zxy,2,\"Interpolated\",1,4"), "CSV omitted interpolation source brackets");
  widget<QPushButton>(*layout, "layoutApply")->click(); QApplication::processEvents();
  MainWindow *copy = nullptr;
  for(auto *top: QApplication::topLevelWidgets())
    if(top->objectName() == "resampledSurveyWindow") copy = qobject_cast<MainWindow *>(top);
  check(copy, "Apply did not open a separate survey window");
  select_station(*copy, "A");
  near(value(*widget<QCustomPlot>(*copy, "plot12"), 1, 2.), std::atan2(3., 2.) * 180. / std::acos(-1.));
  check(widget<QCustomPlot>(window, "plot12")->graph(1)->data()->findBegin(2.)->key != 2., "Resampling modified the original window");
  file_action(*copy, "actionSave_project", {saved});
  copy->close(); QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
  MainWindow restored;
  file_action(restored, "actionLoad_project", {saved}); select_station(restored, "A");
  near(value(*widget<QCustomPlot>(restored, "plot12"), 1, 2.), std::atan2(3., 2.) * 180. / std::acos(-1.));
  widget<QAction>(restored, "actionPeriod_layout")->trigger();
  auto *restoredLayout = widget<QDialog>(restored, "periodLayoutWindow");
  check(!widget<QPushButton>(*restoredLayout, "layoutOpenSource")->isHidden(), "Saved project lost source access");
  widget<QPushButton>(*restoredLayout, "layoutOpenSource")->click(); QApplication::processEvents();
  copy = nullptr;
  for(auto *top: QApplication::topLevelWidgets())
    if(top->objectName() == "resampledSurveyWindow") copy = qobject_cast<MainWindow *>(top);
  check(copy, "Could not reopen the saved source survey");
  select_station(*copy, "A");
  near(value(*widget<QCustomPlot>(*copy, "plot12"), 1, 100.), std::atan2(9., 8.) * 180. / std::acos(-1.));
  copy->close(); QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
}

struct LegacyAxisV5 {
  bool automatic = true;
  double lower = 0., upper = 1.;
  template<class Archive> void serialize(Archive &ar, const unsigned int) { ar & automatic & lower & upper; }
};

struct LegacyOptionsV5 {
  bool phaseWrap = true, showStationNames = true, tipperArrows = false;
  std::vector<LegacyAxisV5> axes = std::vector<LegacyAxisV5>(4);
  std::array<bool, 4> impedance{{false, true, true, false}}, tipper{{true, false, true, false}};
  template<class Archive> void serialize(Archive &ar, const unsigned int) {
    ar & phaseWrap & showStationNames & tipperArrows & axes & impedance & tipper;
  }
};

void workflow(const QString &directory)
{
  const QString observations = directory + "/observed.gofem";
  write(observations,
    "RealZxy 1 Plane_wave S01 .002 .0001\nImagZxy 1 Plane_wave S01 .002 .0001\n"
    "RealZxy 2 Plane_wave S01 .003 .0001\nImagZxy 2 Plane_wave S01 .003 .0001\n"
    "RealZxy 4 Plane_wave S01 .004 .0001\nImagZxy 4 Plane_wave S01 .004 .0001\n"
    "RealZxy 8 Plane_wave S01 .006 .0001\nImagZxy 8 Plane_wave S01 .006 .0001\n"
    "RealZxy 1 Plane_wave S02 .002 .0001\nImagZxy 1 Plane_wave S02 .002 .0001\n");
  auto survey = std::make_shared<MTSurveyData>();
  survey->load_from_gofem(observations.toStdString());
  survey->get_station_data("S01").set_position({{55., 12., 100.}});
  survey->get_station_data("S02").set_position({{55.1, 12.1, 100.}});
  const QString project = directory + "/survey.mtd";
  {
    std::ofstream output(project.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(output);
    archive << survey << std::map<std::string, MTSurveyData>{};
  }
  const QString first = directory + "/inversion_predicted_iter0000.txt";
  const QString second = directory + "/inversion_predicted_iter0001.gofem";
  write(first, "1 S01 impedance_xy real .002 .0001\n1 S01 impedance_xy imag .002 .0001\n"
               "4 S01 impedance_xy real .004 .0001\n4 S01 impedance_xy imag .004 .0001\n");
  write(second, "# frequency_hz receiver observable component value error\n"
                "1 S01 apparent_resistivity_xy value 10 1\n"
                "4 S01 apparent_resistivity_xy value 40 1\n"
                "2 S01 phase_xy value 60 2\n"
                "1 S01 phase_yx value -135 2\n4 S01 phase_yx value -120 2\n"
                "1 S01 phase_tensor_xx value .5 .1\n4 S01 phase_tensor_xx value .75 .1\n"
                "1 S01 induction_x real 0 .01\n4 S01 induction_x real .2 .01\n");
  // A single selection may mix formats; they share every subsequent control.
  {
    std::ifstream input(second.toStdString());
    const auto rows = NativeMT::read_observations(input);
    input.close();
    write(second, gofem_text(rows));
  }
  MainWindow window;
  check(!widget<QCheckBox>(window, "stationMapNames")->isChecked(), "New sessions should start with map names hidden");
  file_action(window, "actionLoad_project", {project});
  select_station(window, "S01");
  auto *stationMap = widget<QCustomPlot>(window, "mapPlot");
  const auto firstPosition = stationMap->graph(0)->data()->constBegin()->value;
  widget<QCheckBox>(window, "stationMapNames")->setChecked(true);
  check(stationMap->itemCount() == 2 && widget<QAction>(window, "actionShow_station_names")->isChecked(),
        "Map name checkbox and menu action are out of sync");
  widget<QAction>(window, "actionShow_station_names")->setChecked(false);
  check(stationMap->itemCount() == 0 && !widget<QCheckBox>(window, "stationMapNames")->isChecked(),
        "Names did not disappear from station map");
  near(stationMap->graph(0)->data()->constBegin()->value, firstPosition);
  widget<QCheckBox>(window, "stationMapNames")->setChecked(true);
  file_action(window, "actionLoad_responses", {first, second});
  auto *list = widget<QListWidget>(window, "responsesList");
  check(list->count() == 2 && list->currentRow() == 1, "Iteration files were not loaded/selected");
  widget<QPushButton>(window, "openPeriodMaps")->click();
  auto *maps = widget<QDialog>(window, "periodMapWindow");
  check(maps->isVisible() && widget<QComboBox>(*maps, "mapDataset")->count() == 3, "Maps button did not open loaded datasets");
  maps->hide();
  widget<QAction>(window, "actionPeriod_maps")->trigger();
  check(maps->isVisible(), "Data menu did not reopen period maps");
  maps->hide();
  auto &rho = *widget<QCustomPlot>(window, "plot11");
  auto &phase = *widget<QCustomPlot>(window, "plot12");
  auto &tipper = *widget<QCustomPlot>(window, "plot21");
  auto &tensor = *widget<QCustomPlot>(window, "plot22");
  const auto beforeReportRange = rho.xAxis->range();
  const auto beforeReportValue = value(rho, 9, 1.);
  export_report(window, directory + "/survey-report.pdf");
  near(rho.xAxis->range().lower, beforeReportRange.lower); near(rho.xAxis->range().upper, beforeReportRange.upper);
  near(value(rho, 9, 1.), beforeReportValue);
  const std::array<SurveyReport::Options::Axis, 4> guiAxes{{
    {false, 1., 1000.}, {false, -10., 100.}, {false, -.4, .4}, {false, -2., 2.}}};
  main_axis_ranges(window, false, guiAxes);
  export_report(window, directory + "/survey-fixed-ranges.pdf", &guiAxes);
  main_axis_ranges(window, true, guiAxes);
  for(auto *plot: {&rho, &phase, &tipper, &tensor})
    for(int i = 8; i < 12; ++i) {
      check(plot->graph(i)->lineStyle() == QCPGraph::lsLine, "Responses must be lines");
      check(plot->graph(i)->scatterStyle().shape() == QCPScatterStyle::ssNone,
            "Responses must not have observation symbols");
      check(plot->graph(i)->selectable() == QCP::stNone, "Responses must not be maskable");
    }
  near(value(rho, 9, 1.), 10.); near(value(rho, 9, .25), 40.);
  check(std::isnan(value(rho, 9, .5)), "Missing response should break the line");
  near(value(phase, 9, .5), 60.); near(value(phase, 10, 1.), 45.);
  near(value(tipper, 8, 1.), 0.); near(value(tensor, 8, .25), .75);
  check(tipper.graph(9)->data()->isEmpty() && tensor.graph(9)->data()->isEmpty(),
        "Absent response component produced a curve");
  widget<QAction>(window, "actionFit_statistics")->trigger();
  auto *analysis = widget<QDialog>(window, "fitStatisticsWindow");
  auto *checks = widget<QTreeWidget>(*analysis, "fitResponses");
  auto *fitPeriods = widget<QCustomPlot>(*analysis, "fitPeriods");
  check(checks->topLevelItemCount() == 2 && fitPeriods->graphCount() == 2, "Missing response comparisons");
  near(checks->topLevelItem(0)->text(1).toDouble(), 0.);
  checks->topLevelItem(0)->setCheckState(0, Qt::Unchecked);
  check(fitPeriods->graphCount() == 1, "Unchecked response remained in fit plot");
  checks->topLevelItem(1)->setCheckState(0, Qt::Unchecked);
  check(fitPeriods->graphCount() == 0, "Unchecking all responses left fit curves");
  checks->topLevelItem(0)->setCheckState(0, Qt::Checked);
  checks->topLevelItem(1)->setCheckState(0, Qt::Checked);
  auto *distributions = widget<QCustomPlot>(*analysis, "fitHistogram");
  for(int i = 0; i < distributions->graphCount(); ++i) {
    double count = 0;
    for(const auto &point: *distributions->graph(i)->data()) count += point.value;
    const auto matched = checks->topLevelItem(i)->text(2).section(" / ", 0, 0).toDouble();
    near(count, matched);
  }
  auto *normalization = widget<QComboBox>(*analysis, "fitErrorSource");
  normalization->setCurrentIndex(1);
  check(checks->topLevelItem(0)->text(2) == "4 / 4", "Normalization lost paired scalars");
  normalization->setCurrentIndex(0);
  // A map edit immediately refreshes both the main curves and the open fit window.
  maps->show(); QApplication::processEvents();
  auto *mapPlot = widget<QCustomPlot>(*maps, "periodMapPlot");
  widget<QComboBox>(*maps, "mapPeriod")->setCurrentIndex(3); // 1 s (four sorted periods)
  near(widget<QComboBox>(*maps, "mapPeriod")->currentData().toDouble(), 1.);
  widget<QCheckBox>(*maps, "mapSelectStations")->setChecked(true);
  const auto stationPoint = widget<QCPGraph>(*mapPlot, "mapStations")->data()->constBegin();
  const QPointF stationPixel(mapPlot->xAxis->coordToPixel(stationPoint->key), mapPlot->yAxis->coordToPixel(stationPoint->value));
  map_gesture(mapPlot, stationPixel, stationPixel);
  widget<QComboBox>(*maps, "mapMaskGroup")->setCurrentIndex(2);
  const auto responsePointCount = rho.graph(9)->data()->size();
  widget<QPushButton>(*maps, "mapMask")->click();
  check(checks->topLevelItem(0)->text(2) == "2 / 4", "Map mask did not refresh statistics");
  check(rho.graph(5)->data()->size() == 1 && phase.graph(5)->data()->size() == 1,
        "Map mask did not refresh observed curves");
  check(rho.graph(9)->data()->size() == responsePointCount, "Map mask modified computed curves");
  widget<QPushButton>(*maps, "mapUnmask")->click();
  check(checks->topLevelItem(0)->text(2) == "4 / 4" && rho.graph(5)->data()->isEmpty(),
        "Map unmask did not restore plots/statistics");
  maps->hide();
  statistics_controls(*analysis);
  const auto pairedCount = checks->topLevelItem(1)->text(2);
  for(int mask: {6, 9, 1, 2, 4, 8, 0, 15}) {
    set_visibility(rho, mask);
    check_visibility(rho, mask);
    for(auto *plot: {&phase, &tensor}) check_visibility(*plot, 15);
    check_visibility(tipper, 15);
  }
  rho.graph(1)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  auto *xy = component_box(rho, 1);
  const QPointF position(6, xy->height() / 2);
  QMouseEvent press(QEvent::MouseButtonPress, position, Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
  QMouseEvent release(QEvent::MouseButtonRelease, position, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(xy, &press);
  QApplication::sendEvent(xy, &release);
  check(!xy->isChecked(), "Clicking a legend checkbox did not hide the component");
  check(rho.graph(1)->selection().isEmpty(), "Hiding a component left maskable selections");
  QKeyEvent spacePress(QEvent::KeyPress, Qt::Key_Space, Qt::NoModifier);
  QKeyEvent spaceRelease(QEvent::KeyRelease, Qt::Key_Space, Qt::NoModifier);
  QApplication::sendEvent(xy, &spacePress);
  QApplication::sendEvent(xy, &spaceRelease);
  check(xy->isChecked() && rho.graph(1)->visible(), "Keyboard cannot restore a hidden component");
  set_visibility(rho, 6);
  set_visibility(phase, 2);
  set_visibility(tensor, 9);
  widget<QAction>(window, "actionShow_error_bars")->setChecked(false);
  check_visibility(rho, 6, false); check_visibility(phase, 2, false); check_visibility(tensor, 9, false);
  widget<QAction>(window, "actionShow_error_bars")->setChecked(true);
  check_visibility(rho, 6); check_visibility(phase, 2); check_visibility(tensor, 9);
  for(int mask: {5, 10, 3, 12, 1, 2, 4, 8, 0, 15}) {
    set_visibility(tipper, mask);
    check_visibility(tipper, mask);
  }
  set_visibility(tipper, 1);
  widget<QAction>(window, "actionTipper_arrows")->setChecked(true);
  check(tipper.graph(0)->visible() && !tipper.graph(2)->visible() && tipper.legend->itemCount() == 2 &&
        tipper.graph(0)->name() == "Real Tzx", "Tipper arrow projection is not identified");
  check(tipper.findChildren<QCheckBox *>().size() == 2, "Arrow legend should have one checkbox per vector");
  bool projectedResponse = false;
  for(int i = 0; i < tipper.itemCount(); ++i)
    if(auto *arrow = dynamic_cast<QCPItemLine *>(tipper.item(i))) {
      near(arrow->end->coords().y(), 0.);
      if(std::abs(arrow->end->coords().x() - 18.) < 1e-10) projectedResponse = true;
    }
  check(projectedResponse, "Selected tipper response projection is missing");
  component_box(tipper, 0)->click();
  check(tipper.itemCount() == 0 && tipper.legend->itemCount() == 2, "Hidden arrows should leave usable legend controls");
  component_box(tipper, 0)->click();
  check(tipper.graph(0)->visible() && tipper.graph(0)->name() == "Real Tzx", "Restoring arrows lost the selected projection");
  widget<QAction>(window, "actionTipper_arrows")->setChecked(false);
  check_visibility(tipper, 1);
  set_visibility(tipper, 5);
  check_visibility(tipper, 5);
  // Refresh the statistics to ensure component visibility never changes masks.
  normalization->setCurrentIndex(1);
  normalization->setCurrentIndex(0);
  check(checks->topLevelItem(1)->text(2) == pairedCount, "Component visibility changed fit coverage");
  const auto pdf = directory + "/fit.pdf";
  file_action(*analysis, "fitSavePdf", {pdf});
  check(QFileInfo(pdf).size() > 1000, "Fit PDF was not written");
  // Checkbox choices survive refresh and statistics react to station masks.
  checks->topLevelItem(1)->setCheckState(0, Qt::Unchecked);
  auto *observedStations = widget<QListWidget>(window, "stationList");
  observedStations->currentItem()->setCheckState(Qt::Unchecked);
  check(checks->topLevelItem(0)->text(2) == "0 / 4" &&
        checks->topLevelItem(1)->checkState(0) == Qt::Unchecked, "Fit refresh ignored station mask or checkbox state");
  observedStations->currentItem()->setCheckState(Qt::Checked);
  check(checks->topLevelItem(0)->text(2) == "4 / 4", "Fit refresh did not restore active station");
  widget<QAction>(window, "actionPhase_wrap")->setChecked(false);
  near(value(phase, 10, 1.), -135.);
  list->setCurrentRow(0);
  near(value(rho, 9, 1.), 1.0132118364233778);
  check(tipper.graph(8)->data()->isEmpty(), "Previous iteration's tipper was left on the plot");
  select_station(window, "S02");
  for(auto *plot: {&rho, &phase, &tipper, &tensor})
    for(int i = 8; i < 12; ++i)
      check(plot->graph(i)->data()->isEmpty(), "Unmatched station retained computed curves");
  select_station(window, "S01");
  check_visibility(rho, 6);
  check_visibility(tipper, 5);
  write(first, "1 S01 apparent_resistivity_xy value 99 1\n4 S01 apparent_resistivity_xy value 100 1\n");
  file_action(window, "actionLoad_responses", {first});
  check(list->count() == 2, "Reloading an iteration duplicated its list item");
  near(value(rho, 9, 1.), 99.);

  // The response table is stored in the project, independently of source files.
  widget<QAction>(window, "actionSave_project")->trigger();
  write(first, "invalid source file\n");
  MainWindow restored;
  file_action(restored, "actionLoad_project", {project});
  check(widget<QCheckBox>(restored, "stationMapNames")->isChecked() &&
        widget<QCustomPlot>(restored, "mapPlot")->itemCount() == 2, "Project lost station-name visibility");
  auto *restoredList = widget<QListWidget>(restored, "responsesList");
  check(restoredList->count() == 2, "Project reload did not restore response list");
  restoredList->setCurrentRow(0);
  select_station(restored, "S01");
  check_visibility(*widget<QCustomPlot>(restored, "plot11"), 6);
  check_visibility(*widget<QCustomPlot>(restored, "plot12"), 2);
  check_visibility(*widget<QCustomPlot>(restored, "plot21"), 5);
  check_visibility(*widget<QCustomPlot>(restored, "plot22"), 9);
  near(value(*widget<QCustomPlot>(restored, "plot11"), 9, 1.), 99.);
  widget<QAction>(restored, "actionFit_statistics")->trigger();
  auto *restoredAnalysis = widget<QDialog>(restored, "fitStatisticsWindow");
  auto *restoredChecks = widget<QTreeWidget>(*restoredAnalysis, "fitResponses");
  check(restoredChecks->topLevelItem(0)->text(2) == "2 / 2", "Project lost native scalar response metadata");
  // Existing GoFEM import remains available and selects the loaded response.
  file_action(restored, "actionLoad_responses", {observations});
  check(restoredList->count() == 3, "GoFEM response import regressed");
  near(value(*widget<QCustomPlot>(restored, "plot11"), 9, 1.), 1.0132118364233778);
  auto checkGoFEM = [&] {
    auto *tree = widget<QTreeWidget>(*restoredAnalysis, "fitResponses");
    QTreeWidgetItem *gofem = nullptr;
    for(int i = 0; i < tree->topLevelItemCount(); ++i)
      if(tree->topLevelItem(i)->data(0, Qt::UserRole).toString() == observations) gofem = tree->topLevelItem(i);
    check(gofem && gofem->text(2) == "10 / 10", "GoFEM statistics did not use the original scalar rows");
    near(gofem->data(1, Qt::UserRole).toDouble(), 0.);
    for(int i = 8; i < 12; ++i)
      check(widget<QCustomPlot>(restored, "plot21")->graph(i)->data()->isEmpty(),
            "GoFEM file without tippers displayed zero response lines");
    gofem->setCheckState(0, Qt::Unchecked);
    check(widget<QCustomPlot>(*restoredAnalysis, "fitPeriods")->graphCount() == 2,
          "GoFEM statistics visibility control failed");
    gofem->setCheckState(0, Qt::Checked);
  };
  checkGoFEM();
  widget<QAction>(restored, "actionSave_project")->trigger();
  write(observations, "source no longer available\n");
  file_action(restored, "actionLoad_project", {project});
  select_station(restored, "S01");
  restoredList->setCurrentRow(2);
  checkGoFEM();

  // Projects saved with the shared toolbar controls retain their visibility.
  const auto legacyProject = directory + "/toolbar-v5.mtd";
  {
    std::ofstream output(legacyProject.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(output);
    const std::uint32_t magic = 0x45444954, version = 5;
    archive << magic << version << survey << std::map<std::string, MTSurveyData>{} << LegacyOptionsV5{};
  }
  file_action(restored, "actionLoad_project", {legacyProject});
  select_station(restored, "S01");
  for(const auto *name: {"plot11", "plot12", "plot22"})
    check_visibility(*widget<QCustomPlot>(restored, name), 6);
  check_visibility(*widget<QCustomPlot>(restored, "plot21"), 5);
}

void real_run(const QString &project, const QString &directory, const QString &screenshot, bool mixed)
{
  MainWindow window;
  file_action(window, "actionLoad_project", {project});
  for(const auto *name: {"plot11", "plot12", "plot21", "plot22"})
    for(int i = 0; i < 4; ++i)
      check(component_box(*widget<QCustomPlot>(window, name), i)->isChecked(),
            "Older projects should show all components by default");
  QStringList paths;
  for(const auto &file: QDir(directory).entryList({"*_predicted_iter*.txt"}, QDir::Files, QDir::Name))
    paths.push_back(QDir(directory).filePath(file));
  check(!paths.isEmpty(), "No inversion files found");
  QTemporaryDir converted;
  if(mixed) {
    check(converted.isValid(), "Cannot create mixed-format fixture directory");
    for(int i = 0; i < paths.size(); ++i) {
      const auto path = converted.filePath(QFileInfo(paths[i]).fileName());
      if(i % 2) {
        std::ifstream input(paths[i].toStdString());
        write(path, gofem_text(NativeMT::read_observations(input)));
      } else check(QFile::copy(paths[i], path), "Cannot copy native fixture");
      paths[i] = path;
    }
    std::cout << "Testing one selection containing alternating native and GoFEM rows, all with .txt filenames.\n";
  }
  file_action(window, "actionLoad_responses", paths);
  auto *stations = widget<QListWidget>(window, "stationList");
  auto *responses = widget<QListWidget>(window, "responsesList");
  check(responses->count() >= paths.size(), "Real iterations failed to import");
  auto &rho = *widget<QCustomPlot>(window, "plot11");
  unsigned plotted = 0;
  for(int i = 0; i < stations->count(); ++i) {
    stations->setCurrentRow(i);
    if(!rho.graph(9)->data()->isEmpty()) ++plotted;
  }
  check(plotted > 0, "No real station names matched");
  select_station(window, "LDG2010_01");
  // Independent conversion of the supplied run's latest impedance row.
  std::ifstream input(paths.back().toStdString());
  double real = 0, imag = 0;
  for(const auto &row: MTResponseData::read(input))
    if(row.receiver == "LDG2010_01" && row.frequency == .01) {
      if(row.type == RealZxy) real = row.value;
      if(row.type == ImagZxy) imag = row.value;
    }
  const double expected = (real*real + imag*imag) / (8e-7 * std::acos(-1.) * std::acos(-1.) * .01);
  near(value(rho, 9, 100.), expected);
  if(!screenshot.isEmpty()) {
    window.showNormal();
    window.resize(1500, 1000);
    widget<QCheckBox>(window, "stationMapNames")->setChecked(false);
    QApplication::processEvents();
    check(window.grab().save(screenshot), "Could not save response preview");
    set_visibility(rho, 6);
    set_visibility(*widget<QCustomPlot>(window, "plot12"), 6);
    set_visibility(*widget<QCustomPlot>(window, "plot22"), 6);
    set_visibility(*widget<QCustomPlot>(window, "plot21"), 5);
    QApplication::processEvents();
    check(window.grab().save(screenshot + ".components.png"), "Could not save component preview");
    const auto checkboxPosition = component_box(rho, 0)->geometry();
    check(rho.savePdf(screenshot + ".legend.pdf", 640, 420), "Could not export checkbox legend");
    check(component_box(rho, 0)->geometry() == checkboxPosition, "Export displaced the on-screen legend controls");
    check(rho.toPixmap(640, 420).save(screenshot + ".legend.png"), "Could not export raster checkbox legend");
  }
  widget<QAction>(window, "actionFit_statistics")->trigger();
  auto *analysis = widget<QDialog>(window, "fitStatisticsWindow");
  auto *checks = widget<QTreeWidget>(*analysis, "fitResponses");
  check(checks->topLevelItemCount() >= paths.size(), "Real statistics are missing iterations");
  // Independently check every loaded nRMS against the solver's history table.
  std::ifstream history(QDir(directory).filePath("inversion_history.csv").toStdString());
  std::map<int, double> expectedRms;
  std::string historyLine;
  std::getline(history, historyLine);
  while(std::getline(history, historyLine)) {
    std::istringstream stream(historyLine);
    std::vector<std::string> columns;
    std::string column;
    while(std::getline(stream, column, ',')) columns.push_back(column);
    if(columns.size() > 7) expectedRms[std::stoi(columns[0])] = std::stod(columns[7]);
  }
  for(int i = 0; i < checks->topLevelItemCount(); ++i) {
    const auto *item = checks->topLevelItem(i);
    const auto filename = QFileInfo(item->data(0, Qt::UserRole).toString()).baseName();
    bool ok = false;
    const int iteration = filename.section("iter", -1).toInt(&ok);
    if(ok && expectedRms.count(iteration)) near(item->data(1, Qt::UserRole).toDouble(), expectedRms.at(iteration));
  }
  std::cout << "Initial/latest nRMS: " << checks->topLevelItem(0)->text(1).toStdString() << "/"
            << checks->topLevelItem(checks->topLevelItemCount()-1)->text(1).toStdString() << "; "
            << checks->topLevelItem(checks->topLevelItemCount()-1)->text(2).toStdString() << " scalars.\n";
  if(!screenshot.isEmpty()) {
    analysis->resize(1600, 1000);
    QApplication::processEvents();
    check(analysis->grab().save(screenshot + ".statistics.png"), "Could not save statistics preview");
    range_dialog(*analysis, [&](QDialog &dialog) {
      check(dialog.grab().save(screenshot + ".ranges.png"), "Could not save range controls preview");
      dialog.reject();
    });
    widget<QComboBox>(*analysis, "fitColorMap")->setCurrentText("Thermal");
    widget<QCheckBox>(*analysis, "fitMapRangeAuto")->setChecked(false);
    widget<QDoubleSpinBox>(*analysis, "fitMapRangeMin")->setValue(.5);
    widget<QDoubleSpinBox>(*analysis, "fitMapRangeMax")->setValue(3.);
    widget<QCheckBox>(*analysis, "fitHeatRangeAuto")->setChecked(false);
    widget<QDoubleSpinBox>(*analysis, "fitHeatRangeMax")->setValue(5.);
    widget<QTabWidget>(*analysis, "fitTabs")->setCurrentIndex(2);
    QApplication::processEvents();
    check(analysis->grab().save(screenshot + ".spatial.png"), "Could not save spatial preview");
    widget<QTabWidget>(*analysis, "fitSpatialTabs")->setCurrentIndex(1);
    QApplication::processEvents();
    check(analysis->grab().save(screenshot + ".heatmaps.png"), "Could not save heatmap preview");
    file_action(*analysis, "fitSavePdf", {screenshot + ".statistics.pdf"});
  }
  widget<QPushButton>(window, "openPeriodMaps")->click();
  auto *maps = widget<QDialog>(window, "periodMapWindow");
  auto *mapPlot = widget<QCustomPlot>(*maps, "periodMapPlot");
  check(mapPlot->property("phaseTensorCount").toUInt() > 0 && mapPlot->property("realVectorCount").toUInt() > 0,
        "Real survey has no mapped tensors or vectors");
  std::cout << "Observed period map: " << widget<QLabel>(*maps, "mapSummary")->text().toStdString() << '\n';
  if(!screenshot.isEmpty()) {
    maps->resize(1200, 1000); QApplication::processEvents();
    check(maps->grab().save(screenshot + ".period-map.png"), "Could not save observed period map");
    file_action(*maps, "mapSavePdf", {screenshot + ".period-map.pdf"});
  }
  auto *mapDataset = widget<QComboBox>(*maps, "mapDataset");
  mapDataset->setCurrentIndex(mapDataset->count() - 1);
  check(mapPlot->property("phaseTensorCount").toUInt() > 0, "Real inversion response has no mapped tensors");
  std::cout << "Response period map: " << widget<QLabel>(*maps, "mapSummary")->text().toStdString() << '\n';
  if(!screenshot.isEmpty()) {
    widget<QCheckBox>(*maps, "mapStationNames")->setChecked(true);
    widget<QCheckBox>(*maps, "mapRealArrows")->setChecked(false);
    widget<QComboBox>(*maps, "mapColorBy")->setCurrentIndex(2);
    widget<QComboBox>(*maps, "mapColorMap")->setCurrentText("Polar");
    widget<QPushButton>(*maps, "mapNextPeriod")->click();
    QApplication::processEvents();
    check(maps->grab().save(screenshot + ".response-map.png"), "Could not save response period map");
  }
  std::cout << paths.size() << " real iterations loaded; " << plotted << " stations plotted.\n";
  if(!screenshot.isEmpty()) export_report(window, screenshot + ".survey.pdf");
  widget<QAction>(window, "actionPeriod_layout")->trigger();
  auto *layout = widget<QDialog>(window, "periodLayoutWindow");
  widget<QPushButton>(*layout, "layoutPreview")->click();
  check(widget<QPushButton>(*layout, "layoutApply")->isEnabled(), "Real survey resampling preview has no usable samples");
  std::cout << "Period layout: " << widget<QLabel>(*layout, "layoutSummary")->text().toStdString() << '\n';
  if(!screenshot.isEmpty()) {
    layout->resize(1300, 950); QApplication::processEvents();
    check(layout->grab().save(screenshot + ".period-layout.png"), "Could not save source layout preview");
    widget<QTabWidget>(*layout, "layoutTabs")->setCurrentIndex(1); QApplication::processEvents();
    check(layout->grab().save(screenshot + ".resampling.png"), "Could not save resampling preview");
    widget<QTabWidget>(*layout, "layoutTabs")->setCurrentIndex(2); QApplication::processEvents();
    check(layout->grab().save(screenshot + ".coverage.png"), "Could not save coverage preview");
  }
}
}

int main(int argc, char **argv)
{
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("EDIToolsTests");
  QCoreApplication::setApplicationName("ResponsePlots");
  QTemporaryDir directory;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
  try {
    check(directory.isValid(), "No temporary directory");
    if(argc >= 3) real_run(argv[1], argv[2], argc >= 4 ? argv[3] : QString(), argc >= 5 && QString(argv[4]) == "mixed");
    else { statistics(directory.path()); gofem_responses(directory.path()); period_maps(directory.path()); period_resampling_tests(directory.path()); resampling_project_workflow(directory.path()); survey_reports(directory.path()); workflow(directory.path()); }
    std::cout << "Response import, curves and project workflow checks passed.\n";
  } catch(const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
