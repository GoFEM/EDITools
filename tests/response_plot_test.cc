#include "include/MapBackground.h"
#include <QDataStream>
#include "mainwindow.h"
#include "include/NativeMT.h"
#include "include/FitStatistics.h"
#include "PeriodMapWindow.h"
#include "PeriodLayoutWindow.h"
#include "include/SurveyReport.h"
#include "include/EDIPeriodMerge.h"
#include "EDIImportDialog.h"
#include "include/ApparentResistivityPlot.h"
#include "include/PhaseTensorPlot.h"
#include "include/PhasePlot.h"
#include "include/TipperPlot.h"
#include "include/StationListDelegate.h"
#include "include/PlotColorMaps.h"
#include <QTabBar>
#include <QSplitter>
#include <QApplication>
#include <QFileDialog>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QEventLoop>
#include <QLineEdit>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QTreeWidget>
#include <QTabWidget>
#include <QSettings>
#include <QTemporaryDir>
#include <QTimer>
#include <QInputDialog>
#include <QRegularExpression>
#include <QProcess>
#include <QStandardPaths>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
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
void edi_import_interpolation_tests(const QString &directory);
void group_mask_tests(const QString &directory);

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
void write_edi(const QString &path, const std::vector<double> &frequencies)
{
  std::ostringstream edi; edi << std::setprecision(17);
  edi << ">HEAD\nDATAID=fixture\nLAT=55\nLONG=9\nELEV=100\nEMPTY=1e32\n>=MTSECT\n>FREQ // " << frequencies.size() << '\n';
  for(double f: frequencies) edi << f << ' ';
  edi << '\n';
  unsigned component = 0;
  for(const auto *name: {"ZXXR", "ZXXI", "ZXYR", "ZXYI", "ZYXR", "ZYXI", "ZYYR", "ZYYI"}) {
    edi << '>' << name << " // " << frequencies.size() << '\n';
    for(unsigned i = 0; i < frequencies.size(); ++i) edi << ++component << ' ';
    edi << '\n';
  }
  for(const auto *name: {"ZXX.VAR", "ZXY.VAR", "ZYX.VAR", "ZYY.VAR", "TXR.EXP", "TXI.EXP", "TYR.EXP", "TYI.EXP", "TXVAR.EXP", "TYVAR.EXP"}) {
    edi << '>' << name << " // " << frequencies.size() << '\n';
    for(unsigned i = 0; i < frequencies.size(); ++i) edi << ".01 ";
    edi << '\n';
  }
  edi << ">END\n"; write(path, edi.str());
}
void edi_import_tests(const QString &directory)
{
  edi_import_interpolation_tests(directory);
  const auto a = directory + "/EDI_A.edi", b = directory + "/EDI_B.edi", c = directory + "/EDI_C.edi";
  write_edi(a, {1., .1}); write_edi(b, {1. / 1.0005, .1}); write_edi(c, {1. / 1.0013});
  std::vector<std::string> paths{a.toStdString(), b.toStdString(), c.toStdString()};
  MTSurveyData source("Import test"); source.load_from_edi(paths);
  const auto names = source.get_stations_names(); const std::set<std::string> all(names.begin(), names.end());
  auto plan = EDIPeriodMerge::analyze(source, all, .001);
  check(plan.before == 4 && plan.after == 3 && plan.mergeable == 1 && plan.changes.size() == 1,
        "EDI clustering must use the full group span, not chained nearest neighbours");
  auto merged = source;
  merged.get_station_data("EDI_B").set_data_mask(RealZxx, 1. / 1.0005, false);
  const auto masks = merged.get_station_data("EDI_B").impedance_mask();
  double z, error, rho, rhoError, result, resultError;
  source.get_station_data("EDI_B").scalar_value(RealZxy, 0, z, error);
  source.get_station_data("EDI_B").scalar_value(RhoZxy, 0, rho, rhoError);
  EDIPeriodMerge::apply(merged, plan);
  const auto &m = merged.get_station_data("EDI_B");
  check(m.frequencies()[0] == 1. && m.frequencies().size() == 2 && m.impedance_mask() == masks, "Merging changed samples or masks");
  check(m.scalar_value(RealZxy, 0, result, resultError) && result == z && resultError == error, "Merging averaged impedance or errors");
  check(m.scalar_value(RhoZxy, 0, result, resultError) && std::abs(result / rho - 1. / 1.0005) < 1e-12, "Resistivity did not follow the assigned period");
  check(source.get_station_data("EDI_B").frequencies()[0] != 1., "Merge mutated the original import");
  check(EDIPeriodMerge::analyze(source, all, 0.).changes.empty(), "Zero tolerance must preserve nearby periods");
  const auto anchored = EDIPeriodMerge::analyze(source, {"EDI_A", "EDI_C"}, .001);
  check(anchored.changes.size() == 1 && anchored.changes[0].station == "EDI_A" && anchored.changes[0].after == 1. / 1.0005, "Existing survey period was not retained as anchor");
  const auto conflict = EDIPeriodMerge::analyze(source, {"EDI_C"}, .005);
  check(conflict.blocked == 1 && conflict.changes.empty(), "Conflicting existing periods were merged");
  const auto collisionPath = directory + "/EDI_collision.edi";
  write_edi(collisionPath, {1., 1. / 1.0004}); paths = {collisionPath.toStdString()}; merged = source; merged.load_from_edi(paths);
  const auto collision = EDIPeriodMerge::analyze(merged, {"EDI_collision"}, .001);
  check(collision.blocked == 1 && collision.changes.empty(), "Two samples from the same station were collapsed");
  // The threshold is strict, including at a binary-exact boundary.
  const auto edgePath = directory + "/EDI_edge.edi";
  write_edi(edgePath, {1., 1. / 1.03125}); paths = {edgePath.toStdString()}; MTSurveyData edge; edge.load_from_edi(paths);
  const auto boundary = EDIPeriodMerge::analyze(edge, {"EDI_edge"}, .03125);
  check(boundary.groups.size()==2 && boundary.mergeable==0 && boundary.changes.empty(), "Exact tolerance boundary was included");

  // Already aligned surveys still need a populated review, not an empty table.
  const auto alignedPath=directory+"/EDI_aligned.edi"; write_edi(alignedPath,{1.,.1});
  MTSurveyData aligned; aligned.load_from_edi({a.toStdString(),alignedPath.toStdString()});
  const auto alignedNames=aligned.get_stations_names();
  QTimer::singleShot(0,[] {
    auto *dialog=qobject_cast<QDialog *>(QApplication::activeModalWidget()); check(dialog,"Import review missing");
    auto *periods=widget<QTableWidget>(*dialog,"ediPeriodGroups");
    check(periods->rowCount()==2,"Exact shared periods missing from import review");
    for(int r=0;r<2;++r) {
      check(periods->item(r,3)->text()=="1" && periods->item(r,4)->text()=="2","Shared-period counts incorrect");
      check(periods->item(r,5)->text()=="Keep (already shared)","Aligned periods offered a merge");
    }
    check(widget<QLabel>(*dialog,"ediPeriodStatus")->text().contains("No close periods to merge"),"No-merge explanation missing");
    dialog->accept();
  });
  check(reviewEDIImport(nullptr,aligned,{alignedNames.begin(),alignedNames.end()},{}),"Aligned import could not proceed");
  check(aligned.get_station_data("EDI_A").frequencies()==std::vector<double>({1.,.1}),"Aligned review changed frequencies");

  MainWindow window;
  auto import = [&](const QStringList &files, bool accept, int periodMode) {
    QTimer timer; bool reviewed = false;
    QObject::connect(&timer, &QTimer::timeout, [&] {
      for(auto *top: QApplication::topLevelWidgets()) if(auto *name = qobject_cast<QInputDialog *>(top)) {
        if(name->isVisible()) name->accept();
      }
      auto *dialog = window.findChild<QDialog *>("ediImportDialog");
      if(!dialog || !dialog->isVisible()) return;
      check(widget<QComboBox>(*dialog, "ediPeriodMode")->currentIndex() == 0, "EDI period changes must be opt-in");
      check(widget<QDoubleSpinBox>(*dialog, "ediPeriodTolerance")->value() == .1, "Default EDI tolerance must be 0.1%");
      check(widget<QLabel>(*dialog, "ediImportSummary")->text().contains("distinct periods"), "Import period statistics missing");
      widget<QComboBox>(*dialog, "ediPeriodMode")->setCurrentIndex(periodMode);
      reviewed = true; timer.stop(); if(accept) dialog->accept(); else dialog->reject();
    });
    timer.start(0); file_action(window, "actionLoad_EDI", files); check(reviewed, "EDI import did not show review dialog");
  };
  import({a, b}, false, true);
  check(widget<QListWidget>(window, "stationList")->count() == 0, "Canceled import changed the survey");
  import({a, b}, true, true);
  check(widget<QListWidget>(window, "stationList")->count() == 2, "Accepted EDI stations missing");
  auto *list = widget<QListWidget>(window, "stationList"); list->setCurrentRow(1);
  const auto graph = widget<QCustomPlot>(window, "plot11")->graph(1)->data();
  check(graph->size() == 2 && graph->constBegin()->key == 1., "GUI import did not apply the period merge");
  const auto title = window.windowTitle(); import({c}, false, true);
  check(list->count() == 2 && window.windowTitle() == title && list->currentRow() == 1, "Canceled append changed the open survey");
  import({c}, true, false); list->setCurrentRow(2);
  check(list->count() == 3 && widget<QCustomPlot>(window, "plot11")->graph(1)->data()->constBegin()->key == 1.0013,
        "Keep-original import changed a period");
  const auto interpolatedPath = directory + "/EDI_interpolated.edi";
  write_edi(interpolatedPath, {2., 1., .5, .1});
  import({interpolatedPath}, false, 2);
  check(list->count() == 3, "Canceled interpolation added a station");
  import({interpolatedPath}, true, 2);
  auto matches = list->findItems("EDI_interpolated", Qt::MatchExactly);
  check(matches.size() == 1, "Interpolated station missing"); list->setCurrentItem(matches.front());
  auto data = widget<QCustomPlot>(window, "plot11")->graph(1)->data();
  check(data->size() == 3 && std::abs(data->constBegin()->key - 1.) < 1e-12 &&
        std::abs((data->constBegin()+1)->key - 1.0013) < 1e-12, "Import did not use the existing-period union");
  widget<QAction>(window, "actionUndo")->trigger(); check(list->count() == 3, "Cannot undo interpolated import");
  widget<QAction>(window, "actionRedo")->trigger(); check(list->count() == 4, "Cannot redo interpolated import");
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
                   const std::array<SurveyReport::Options::Axis, 4> *expectedAxes = nullptr,
                   const std::function<void(QDialog &)> &configure = {}, int expectedPages = -1)
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
    check(widget<QCheckBox>(*dialog, "reportOverview")->toolTip().contains("overview"), "Report overview tooltip missing");
    auto *save = widget<QDialogButtonBox>(*dialog, "")->button(QDialogButtonBox::Save);
    widget<QCheckBox>(*dialog, "reportOverview")->setChecked(false);
    widget<QCheckBox>(*dialog, "reportStationPages")->setChecked(false);
    check(!save->isEnabled(), "Empty report sections accepted");
    widget<QCheckBox>(*dialog, "reportPTMaps")->setChecked(true);
    widget<QLineEdit>(*dialog, "reportMapPeriods")->setText("bad");
    check(!save->isEnabled(), "Invalid map period accepted");
    widget<QLineEdit>(*dialog, "reportMapPeriods")->setText("1, 10, 1");
    check(save->isEnabled() && widget<QLabel>(*dialog, "reportPageCount")->text() == "2 pages", "Map page count did not deduplicate periods");
    widget<QCheckBox>(*dialog, "reportInductionMaps")->setChecked(true);
    check(widget<QLabel>(*dialog, "reportPageCount")->text() == "2 pages", "Combined map layers added separate pages");
    widget<QCheckBox>(*dialog, "reportInductionMaps")->setChecked(false);
    widget<QCheckBox>(*dialog, "reportPTMaps")->setChecked(false);
    widget<QCheckBox>(*dialog, "reportOverview")->setChecked(true);
    widget<QCheckBox>(*dialog, "reportStationPages")->setChecked(true);
    if(configure) configure(*dialog);
    handled = true; timer.stop(); dialog->accept();
  });
  // Exercise the default suffix when the remaining filename has no other dots.
  const auto enteredPath = QFileInfo(path).completeBaseName().contains('.') ? path : path.chopped(4);
  timer.start(10); file_action(window, "actionExport_survey_report", {enteredPath});
  check(handled, "Report options did not open");
  check(!expectedAxes || inspected > 0, "GUI fixed ranges were not checked during export");
  report_pages(path, expectedPages < 0 ? stations->count() + 1 : expectedPages);
  check(stations->currentRow() == stationRow && responses->currentRow() == responseRow, "Report changed current selection");
}
void main_axis_ranges(MainWindow &window, bool autoscale, const std::array<SurveyReport::Options::Axis, 4> &axes)
{
  widget<QAction>(window, "actionPlot_axis_ranges")->trigger();
  check(!QApplication::activeModalWidget(), "Axis settings opened a modal window");
  for(unsigned i = 0; i < 2; ++i) {
    widget<QComboBox>(window, "viewRangeAxis")->setCurrentIndex(i + 1);
    widget<QLineEdit>(window, "viewRangeMin")->setText(QString::number(axes[i].lower));
    widget<QLineEdit>(window, "viewRangeMax")->setText(QString::number(axes[i].upper));
    widget<QPushButton>(window, "viewRangeApply")->click();
    widget<QCheckBox>(window, "viewRangeAuto")->setChecked(autoscale);
  }
}

void map_backgrounds(const QString &directory)
{
  struct EnvironmentGuard {
    QByteArray previous = qgetenv("EDITOOLS_MAP_DATA_DIR");
    ~EnvironmentGuard() { if(previous.isNull()) qunsetenv("EDITOOLS_MAP_DATA_DIR"); else qputenv("EDITOOLS_MAP_DATA_DIR", previous); }
  } restoreEnvironment;
  qputenv("EDITOOLS_MAP_DATA_DIR", (directory + "/no-map-data").toUtf8());
  check(MapBackground::availableLayers() == 0, "Missing external data should disable geography");
  QWidget parent;
  auto *menu = MapBackground::button(&parent, 15, [](int) {})->menu();
  for(auto *action: menu->actions()) if(action->isCheckable() && action->data().toInt() <= 8) check(!action->isEnabled() && !action->isChecked(), "Unavailable layer action enabled");
  // Synthetic geometry keeps tests independent of optional third-party datasets.
  const auto maps = directory + "/external-maps"; check(QDir().mkpath(maps), "Cannot create test maps");
  QByteArray raw; QDataStream output(&raw, QIODevice::WriteOnly);
  output.setByteOrder(QDataStream::LittleEndian); output.setFloatingPointPrecision(QDataStream::SinglePrecision);
  output << quint32(1) << quint32(3) << float(-113.1) << float(37.9) << float(-113.) << float(38.) << float(-112.9) << float(38.1);
  for(const auto *name: {"coast", "countries", "states", "rivers", "lakes"}) {
    QFile file(maps + "/" + name + ".bin"); check(file.open(QIODevice::WriteOnly), "Cannot write test map"); file.write(qCompress(raw));
  }
  const auto corrupt = directory + "/corrupt-maps"; check(QDir().mkpath(corrupt), "Cannot create corrupt fixture");
  QFile broken(corrupt + "/coast.bin"); check(broken.open(QIODevice::WriteOnly), "Cannot create broken geometry");
  broken.write(qCompress(raw.left(8))); broken.close();
  qputenv("EDITOOLS_MAP_DATA_DIR", corrupt.toUtf8());
  check(MapBackground::availableLayers() == 0, "Truncated external data accepted");
  qputenv("EDITOOLS_MAP_DATA_DIR", maps.toUtf8());
  check(MapBackground::availableLayers() == 15, "External geometry not detected");

  QCustomPlot plot; plot.resize(800, 600);
  auto *background = new MapBackground(&plot);
  plot.xAxis->setRange(-125., -65.); plot.yAxis->setRange(25., 55.);
  auto render = [&] { return plot.toPixmap(800, 600).toImage(); };
  check(background->layers() == 15, "Available offline layers are not enabled by default");
  auto *defaults = MapBackground::button(&parent, background->layers(), [](int) {})->menu();
  for(auto *action: defaults->actions()) if(action->isCheckable())
    check(action->isChecked() == (action->data().toInt() <= 8), "Default layer controls differ from the map or enable online tiles");
  background->setLayers(0);
  const auto empty = render(); const auto x = plot.xAxis->range(), y = plot.yAxis->range();
  for(int layer: {1, 2, 4, 8}) {
    background->setLayers(layer);
    check(render() != empty, "External geography layer rendered no lines");
    check(plot.xAxis->range() == x && plot.yAxis->range() == y && plot.graphCount() == 0 && plot.itemCount() == 0,
          "Geography changed plot ranges or data items");
  }
  background->setLayers(0); check(render() == empty, "Disabling layers left geography behind");
  // A known state-boundary vertex must align with the station projection in
  // degrees, UTM metres, and centred UTM kilometres, including a southern zone.
  QFile file(maps + "/states.bin"); check(file.open(QIODevice::ReadOnly), "Missing test states");
  QDataStream input(qUncompress(file.readAll())); input.setByteOrder(QDataStream::LittleEndian); input.setFloatingPointPrecision(QDataStream::SinglePrecision);
  quint32 count; input >> count; std::array<double, 3> location{{NAN, NAN, 0.}};
  for(quint32 i = 0; i < count && !std::isfinite(location[0]); ++i) {
    quint32 size; input >> size;
    for(quint32 j = 0; j < size; ++j) {
      float lon, lat; input >> lon >> lat;
      if(lat > 36. && lat < 40. && lon > -115. && lon < -111.) location = {{lat, lon, 0.}};
    }
  }
  check(std::isfinite(location[0]), "No test boundary vertex");
  for(int mode = 0; mode < 5; ++mode) {
    SurveyCoordinates coordinates;
    // Continental surveys can place a valid station more than 30 degrees from
    // the common zone's central meridian; its background must still align.
    if(mode) coordinates = SurveyCoordinates::calculate({location}, mode == 4 ? 19 : 12, mode != 3, mode >= 2);
    const double scale = mode >= 2 ? .001 : 1.;
    const auto point = coordinates.transform({location}).front();
    const double extent = mode == 0 ? .015 : (mode >= 2 ? 1. : 1000.);
    plot.xAxis->setRange(point[1] * scale - extent, point[1] * scale + extent);
    plot.yAxis->setRange(point[0] * scale - extent, point[0] * scale + extent);
    background->setCoordinates(coordinates, scale); background->setLayers(MapBackground::States);
    plot.clearPlottables(); plot.clearItems(); const auto image = render();
    const int px = qRound(plot.xAxis->coordToPixel(point[1] * scale)), py = qRound(plot.yAxis->coordToPixel(point[0] * scale));
    bool hit = false;
    for(int dx = -3; dx <= 3; ++dx) for(int dy = -3; dy <= 3; ++dy) {
      const auto color = image.pixelColor(px + dx, py + dy);
      hit |= color.blue() > color.red() + 3 && color.red() < 230;
    }
    check(hit, "Boundary vertex misaligned with coordinate frame/origin/units");
  }
  check(plot.savePdf(directory + "/geography.pdf"), "Geography PDF failed");
}

void checkNear(double actual, double expected);

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
  auto selected = options; selected.overview = false; selected.stationPages = false;
  selected.phaseTensorMaps = true; selected.inductionMaps = true; selected.mapPeriods = {1., .5, 1.};
  selected.mapData = SurveyReport::Options::MapData::Both; selected.mapSettings.imaginary = true;
  unsigned mapPages = 0;
  check(SurveyReport::write_pdf(directory + "/report-maps.pdf", survey, &response, selected, [&](unsigned done, unsigned total) {
    check(total == 4, "Map periods should deduplicate and combine layers for each dataset");
    if(done) {
      ++mapPages; QCustomPlot *plot = nullptr;
      for(auto *top: QApplication::topLevelWidgets()) if(top->objectName() == "periodMapWindow") plot = top->findChild<QCustomPlot *>("periodMapPlot");
      check(plot, "Map renderer missing");
      check(plot->findChild<QCPItemLine *>("mapReferenceArrow"), "Combined map omitted the induction reference");
    }
    return true;
  }), "Map-only export failed");
  check(mapPages == 4, "Map-only progress missing pages"); report_pages(directory + "/report-maps.pdf", 4);
  // Inspect actual PDF text and metadata when Poppler tools are available.
  // Include all page types and paths from both platforms, regardless of host OS.
  const auto pdfText = QStandardPaths::findExecutable("pdftotext"), pdfInfo = QStandardPaths::findExecutable("pdfinfo");
  if(!pdfText.isEmpty() && !pdfInfo.isEmpty()) {
    auto privacy = selected; privacy.overview = true; privacy.stationPages = true; privacy.mapPeriods = {1.};
    for(const auto &label: QStringList{"/home/private-report-user/results/response.dat",
        "C:\\Users\\private-report-user\\results\\response.dat", "\\\\private-report-host\\results\\response.dat"}) {
      privacy.responseName = label;
      const auto path = directory + "/private-report-output.pdf";
      check(SurveyReport::write_pdf(path, survey, &response, privacy), "Privacy report export failed");
      for(const auto &tool: {pdfText, pdfInfo}) {
        QProcess process;
        process.start(tool, tool == pdfText ? QStringList{path, "-"} : QStringList{path});
        check(process.waitForFinished(10000) && process.exitStatus() == QProcess::NormalExit && process.exitCode() == 0, "Cannot inspect PDF privacy");
        const auto content = process.readAllStandardOutput();
        check(!content.contains("private-report-") && !content.contains("Source:") && !content.contains(directory.toUtf8()) &&
              !content.contains("/home/") && !content.contains("C:\\Users") && !content.contains("/mnt/"), "Local path or machine details leaked into PDF");
        if(tool == pdfText) check(content.contains("response.dat"), "Public response label missing from PDF");
      }
    }
  }
  selected.mapPeriods = {0.}; bool invalid = false;
  try { SurveyReport::write_pdf(enabled, survey, &response, selected); } catch(const std::invalid_argument &) { invalid = true; }
  check(invalid && read_report(enabled) == original, "Invalid map selection replaced destination");
  selected.mapPeriods = {1.}; invalid = false;
  try { SurveyReport::write_pdf(enabled, survey, nullptr, selected); } catch(const std::invalid_argument &) { invalid = true; }
  check(invalid, "Response maps accepted absent response");
  auto arrows = fixedOptions; arrows.overview = false; arrows.tipperArrows = true;
  arrows.components[2][1] = false; arrows.components[2][3] = false;
  check(SurveyReport::write_pdf(directory + "/report-arrows.pdf", survey, &response, arrows), "Repeated arrow station pages failed");
  report_pages(directory + "/report-arrows.pdf", 2);
  auto lines = options; lines.overview = false;
  check(SurveyReport::write_pdf(directory + "/report-lines.pdf", survey, nullptr, lines), "Station-only line export failed");
  report_pages(directory + "/report-lines.pdf", 1);
  const auto gapsPath = directory + "/tipper-gaps.gofem";
  write(gapsPath, "RealTzx 1 Plane_wave GAP .1 .01\nRealTzx 2 Plane_wave GAP .2 .02\nRealTzx 4 Plane_wave GAP .3 .03\n");
  MTSurveyData gaps; gaps.load_from_gofem(gapsPath.toStdString()); gaps.get_station_data("GAP").set_data_mask(RealTzx, 2., false);
  bool checkedGaps = false;
  check(SurveyReport::write_pdf(directory + "/tipper-gaps.pdf", gaps, nullptr, lines, [&](unsigned done, unsigned) {
    if(!done) return true;
    for(auto *top: QApplication::topLevelWidgets()) if(top->objectName() == "surveyReportPlot2") {
      auto *plot = dynamic_cast<QCustomPlot *>(top); const auto data = plot->graph(0)->data();
      check(data->size() == 3 && std::isnan(data->at(1)->value), "Masked tipper period did not break its line");
      for(int i = 0; i < plot->plottableCount(); ++i) {
        auto *bars = dynamic_cast<QCPErrorBars *>(plot->plottable(i));
        if(bars && bars->dataPlottable() == plot->graph(0)) {
          check(bars->data()->size() == 3 && std::isnan(bars->data()->at(1).errorPlus), "Mask gap shifted error bars");
          checkNear(bars->data()->at(0).errorPlus, .03); checkNear(bars->data()->at(2).errorPlus, .01); checkedGaps = true;
        }
      }
    }
    return true;
  }), "Tipper gap export failed"); check(checkedGaps, "Tipper line uncertainties were not inspected");
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
void checkNear(double actual, double expected)
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
  const auto boxes = plot.findChildren<QCheckBox *>();
  for(auto *box: boxes) {
    const int component = box->objectName().mid(QStringLiteral("componentVisibility").size()).toInt();
    const int bits = boxes.size() == 1 ? 15 : boxes.size() == 2 ? 3 << component : 1 << component;
    if(box->isChecked() != bool(mask & bits)) box->click();
  }
}

void check_visibility(QCustomPlot &plot, int mask, bool errorBars = true)
{
  const int checkboxCount = plot.findChildren<QCheckBox *>().size();
  const bool arrows = plot.graph(0)->property("tipperArrowRole").toInt() == 1;
  const bool ellipses = checkboxCount == 1;
  const int visibleMask = arrows ? (mask & 3 ? 1 : 0) | (mask & 12 ? 4 : 0) : ellipses ? (mask ? 1 : 0) : mask;
  for(int i = 0; i < plot.graphCount(); ++i) {
    check(plot.graph(i)->visible() == (bool(visibleMask & (1 << (i % 4))) && (!(arrows || ellipses) || i < 8)), "Component visibility is inconsistent");
    check(plot.graph(i)->selectable() == (i < 8 && plot.graph(i)->visible()
          ? QCP::stMultipleDataRanges : QCP::stNone), "Hidden or calculated curves are selectable");
  }
  int errorCount = 0;
  for(int i = 0; i < plot.plottableCount(); ++i)
    if(auto *bars = dynamic_cast<QCPErrorBars *>(plot.plottable(i))) {
      ++errorCount;
      check(bars->visible() == (errorBars && !arrows && !ellipses && bars->dataPlottable()->visible()), "Hidden component has visible error bars");
    }
  check(errorCount == 8, "Error bars were lost");
  const auto legends = plot.findChildren<QCPLegend *>();
  int legendItems = 0;
  for(auto *legend: legends) legendItems += legend->itemCount();
  check(legendItems == checkboxCount && checkboxCount == (arrows ? 2 : ellipses ? 1 : 4),
        "Each legend item must keep exactly one visibility checkbox");
  for(int i = 0; i < 4; ++i) {
    if((arrows && i % 2) || (ellipses && i != 0)) continue;
    auto *box = component_box(plot, i);
    check(box->isChecked() == bool(visibleMask & (1 << i)), "Legend checkbox is out of sync");
    check(!box->isHidden() && box->isEnabled() && std::any_of(legends.begin(), legends.end(),
          [box](QCPLegend *legend) { return legend->outerRect().contains(box->geometry()); }),
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
    checkNear(periods->xAxis->range().lower, .1); checkNear(periods->xAxis->range().upper, 10.);
    checkNear(periods->yAxis->range().lower, .5); checkNear(periods->yAxis->range().upper, 3.);
    auto *histogram = widget<QCustomPlot>(analysis, "fitHistogram");
    checkNear(histogram->xAxis->range().lower, -2.); checkNear(histogram->xAxis->range().upper, 2.);
    checkNear(widget<QCustomPlot>(analysis, "fitOverall")->yAxis->range().upper, 4.);
    checkNear(widget<QCustomPlot>(analysis, "fitComponents")->yAxis->range().upper, 4.);
    checkNear(widget<QCustomPlot>(analysis, "fitStations")->yAxis->range().upper, 5.);
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
  checkNear(checks->topLevelItem(1)->data(1, Qt::UserRole).toDouble(), originalRms);
  check(checks->topLevelItem(1)->text(2) == originalCount, "Display controls changed the comparison counts");
  auto *histogram = widget<QCustomPlot>(analysis, "fitHistogram");
  double retained = 0.;
  for(const auto &point: *histogram->graph(1)->data()) retained += point.value;
  checkNear(retained, originalCount.section(" / ", 0, 0).toDouble());
  range_dialog(analysis, [](QDialog &dialog) {
    widget<QPushButton>(dialog, "fitAxesAutoscale")->click();
    dialog.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->click();
  });
  checkNear(periods->yAxis->range().lower, 0.);
  check(periods->yAxis->range().upper != 3., "Autoscale failed to clear fixed limits");
  widget<QCheckBox>(analysis, "fitMapRangeAuto")->setChecked(true);
  widget<QCheckBox>(analysis, "fitHeatRangeAuto")->setChecked(true);
  auto *scale = dynamic_cast<QCPColorScale *>(widget<QCustomPlot>(analysis, "fitMapA")->plotLayout()->element(1, 1));
  checkNear(scale->dataRange().lower, 0.);
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
  checkNear(result.total.rms(), 5.);
  check(result.residuals == std::vector<double>({5., 5.}), "Scalar response errors or residual sign changed");
  checkNear(result.periods.at(1.).rms(), 5.);
  checkNear(result.components.at("Zxy").rms(), 5.);
  checkNear(result.stations.at("S01").rms(), 5.);
  checkNear(result.station_periods.at("S01").at(1.).rms(), 5.);
  result = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Observed);
  checkNear(result.total.rms(), std::sqrt(250.));
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
  checkNear(result.residuals.at(0), -.5);
  check(result.periods.count(1.) == 1, "Rounded frequency did not align with observed periods");
  write(observedPath, "1 S01 phase_xy value 45 2\n");
  write(responsePath, "1 S01 phase_xy value 405 2\n");
  observed.load_from_native_responses(observedPath.toStdString());
  response.load_from_native_responses(responsePath.toStdString());
  checkNear(FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response).total.rms(), 0.);
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
  checkNear(response.response_observations()[0].error, .4);
  checkNear(response.response_observations()[3].error, .2);
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
  checkNear(fit.total.rms(), 5.);
  check(fit.residuals == std::vector<double>({5., -5., 5.}), "GoFEM merged real/imaginary response errors");
  fit = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Observed);
  check(fit.total.count == 4, "Zero-error GoFEM response cannot use observed errors");
  checkNear(fit.total.rms(), std::sqrt(175.));
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., false);
  fit = FitStatistics::compare(observed, response, FitStatistics::ErrorSource::Response);
  check(fit.total.count == 1, "GoFEM statistics ignored observed masks");
  observed.get_station_data("S01").set_data_mask(RealZxy, 1., true);

  // Both missing values and individual scalar errors must survive a project archive.
  std::stringstream saved;
  { boost::archive::binary_oarchive archive(saved); archive << response; }
  MTSurveyData restored;
  { boost::archive::binary_iarchive archive(saved); archive >> restored; }
  checkNear(FitStatistics::compare(observed, restored, FitStatistics::ErrorSource::Response).total.rms(), 5.);
  check(restored.response_observations().size() == 4 &&
        !restored.get_station_data("S01").scalar_value(ImagTzx, 0, scalar, error), "Project lost sparse GoFEM data");

  // Reloads replace all prior rows; malformed input leaves the current data intact.
  write(path, "RealTzy 1e-8 Plane_wave station:03 0 .01\nRealTzy 1.005e-8 Plane_wave station:03 .1 .02\n");
  response.load_from_gofem(path.toStdString());
  check(response.n_stations() == 1 && !response.is_station_present("S01") && response.response_observations().size() == 2,
        "GoFEM reload retained old response data");
  check(response.get_station_data("station:03").scalar_value(RealTzy, 1, scalar, error), "GoFEM collapsed distinct small frequencies");
  checkNear(scalar, .1);
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
    checkNear(scalar, row.value); checkNear(scalar, expectedValue); checkNear(error, expectedError);
  }
  const auto &preserved = response.response_observations();
  for(unsigned i = 0; i < rows.size(); ++i) {
    checkNear(preserved[i].value, rows[i].value); checkNear(preserved[i].error, rows[i].error);
  }
}

void linked_mask_tests(const QString &directory)
{
  const auto path = directory + "/linked.edi"; write_edi(path, {1., .1});
  auto survey = std::make_shared<MTSurveyData>(); std::vector<std::string> paths{path.toStdString()}; survey->load_from_edi(paths);
  const auto original = survey->get_station_data("linked");
  const std::array<RealDataType, 4> zTypes{{RealZxx, RealZxy, RealZyx, RealZyy}}, ptTypes{{PTxx, PTxy, PTyx, PTyy}};
  // Exercise real rectangle selection: the plotting library tests hidden
  // graphs too unless they are explicitly made unselectable.
  for(bool linked: {false, true}) {
    auto edited = original;
    QCustomPlot plot; ApparentResistivityPlot handler(&plot);
    handler.set_link_tensor_masks(linked);
    handler.set_observed_data(edited);
    plot.resize(800, 600); plot.show(); QApplication::processEvents();
    set_visibility(plot, 3); // XX and XY only
    auto selectPeriod = [&](double period) {
      plot.replot();
      const auto rect = plot.axisRect()->rect();
      map_gesture(&plot, QPointF(plot.xAxis->coordToPixel(period / 1.1), rect.top() + 2),
                  QPointF(plot.xAxis->coordToPixel(period * 1.1), rect.bottom() - 2), true);
      check(plot.graph(0)->selected() || plot.graph(4)->selected(), "Rectangle missed visible XX data");
      for(int c: {2, 3, 6, 7}) check(!plot.graph(c)->selected(), "Rectangle selected a hidden component");
    };
    selectPeriod(1.); handler.maskSelectedData();
    for(unsigned c = 0; c < 4; ++c) {
      check(edited.impedance_mask()[c][0] == (c >= 2) && edited.impedance_mask()[c][1],
            "Rectangle mask changed hidden components or another period");
      check(edited.phase_tensor_mask()[c][0] == !linked, "Rectangle mask lost PT linkage");
    }
    selectPeriod(10.); handler.maskSelectedData();
    selectPeriod(1.); handler.unmaskSelectedData();
    selectPeriod(1.); handler.invMaskSelectedData();
    for(unsigned c = 0; c < 4; ++c)
      check(edited.impedance_mask()[c][0] && edited.impedance_mask()[c][1] == (c >= 2),
            "Repeated or inverted masking changed hidden components");
    check(edited.tipper_mask() == original.tipper_mask(), "Rectangle masking changed tippers");
    set_visibility(plot, 15);
    for(int c = 0; c < 8; ++c) check(plot.graph(c)->selectable() == QCP::stMultipleDataRanges,
                                    "Shown observations did not become selectable again");
    for(int c = 8; c < plot.graphCount(); ++c) check(plot.graph(c)->selectable() == QCP::stNone,
                                                  "Response curves became selectable");
  }
  // Real rectangle selections in both TF columns, including hidden paired
  // diagonals and masked observations, must use the same component/period plan.
  for(bool phase: {false, true}) for(bool linked: {false, true}) for(int off: {1, 2}) {
    auto edited = original;
    QCustomPlot plot; plot.resize(1000, 600);
    std::unique_ptr<MTDataPlot> handler(phase ? static_cast<MTDataPlot *>(new PhasePlot(&plot))
                                            : new ApparentResistivityPlot(&plot));
    handler->set_split_components(true); handler->set_link_tensor_masks(linked);
    handler->set_observed_data(edited); handler->set_predicted_data(original);
    plot.show(); QApplication::processEvents();
    const int diagonal = off == 1 ? 0 : 3;
    handler->set_mask_preview({{zTypes[off], 1.}, {zTypes[diagonal], 1.}});
    const auto markers = plot.findChildren<QCPItemTracer *>("autoMaskMarker");
    check(markers.size() == 2, "TF preview omitted the paired diagonal");
    std::set<QCPAxis *> markerAxes;
    for(auto *marker: markers) {
      markerAxes.insert(marker->position->keyAxis());
      check(marker->clipAxisRect() == marker->position->keyAxis()->axisRect(), "TF marker clipped to the wrong column");
      check(marker->position->keyAxis() == plot.graph(off)->keyAxis() ||
            marker->position->keyAxis() == plot.graph(diagonal)->keyAxis(), "TF marker uses an unrelated axis");
    }
    check(markerAxes.size() == 2, "TF previews did not appear in both columns");
    handler->set_mask_preview({});
    set_visibility(plot, 1 << off);
    auto selectPeriod = [&](int component, double period) {
      plot.replot();
      auto *axis = plot.graph(component)->keyAxis();
      const auto rect = axis->axisRect()->rect();
      map_gesture(&plot, QPointF(axis->coordToPixel(period / 1.1), rect.top() + 2),
                  QPointF(axis->coordToPixel(period * 1.1), rect.bottom() - 2), true);
      check(plot.graph(component)->selected() || plot.graph(component + 4)->selected(), "TF rectangle missed its component");
      for(int c = 0; c < 12; ++c) if(c != component && c != component + 4)
        check(!plot.graph(c)->selected(), "TF rectangle selected another column, hidden component or response");
    };
    selectPeriod(off, 1.); handler->maskSelectedData();
    for(unsigned c = 0; c < 4; ++c) {
      check(edited.impedance_mask()[c][0] == (c != off && c != diagonal) && edited.impedance_mask()[c][1],
            "TF mask missed a hidden paired diagonal or changed another component/period");
      check(edited.phase_tensor_mask()[c][0] == !linked && edited.phase_tensor_mask()[c][1], "TF mask ignored Z/PT linkage");
    }
    selectPeriod(off, 1.); handler->unmaskSelectedData();
    check(edited.impedance_mask() == original.impedance_mask() && edited.phase_tensor_mask() == original.phase_tensor_mask(),
          "TF unmask did not restore the paired diagonal and linked tensor");
    selectPeriod(off, 1.); handler->invMaskSelectedData();
    for(unsigned c = 0; c < 4; ++c) {
      check(edited.impedance_mask()[c][0] && edited.impedance_mask()[c][1] == (c != off && c != diagonal), "TF inverted mask missed the paired diagonal");
      check(edited.phase_tensor_mask()[c][0] && edited.phase_tensor_mask()[c][1] == !linked, "TF inverted mask lost Z/PT linkage");
    }
    check(edited.tipper_mask() == original.tipper_mask(), "TF impedance edit changed tipper masks");
    edited = original; handler->set_observed_data(edited); set_visibility(plot, 1 << diagonal);
    selectPeriod(diagonal, 1.); handler->maskSelectedData();
    check(!edited.impedance_mask()[diagonal][0] && edited.impedance_mask()[off][0], "Diagonal rectangle masked its off-diagonal partner");
    // Conflicting selections in the two columns: a mask must win for the
    // shared diagonal, regardless of component order or displayed quantity.
    edited = original; handler->set_observed_data(edited); set_visibility(plot, 15); plot.deselectAll();
    plot.graph(diagonal)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
    plot.graph(off)->setSelection(QCPDataSelection(QCPDataRange(1, 2))); handler->invMaskSelectedData();
    check(!edited.impedance_mask()[diagonal][0] && !edited.impedance_mask()[diagonal][1] &&
          !edited.impedance_mask()[off][0] && edited.impedance_mask()[off][1], "Conflicting TF inversion depends on component order");
  }
  {
    auto edited = original;
    QCustomPlot plot; TipperPlot handler(&plot); handler.set_arrow_mode(true); handler.set_observed_data(edited);
    plot.graph(0)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
    plot.graph(2)->setSelection(QCPDataSelection(QCPDataRange(1, 2))); handler.invMaskSelectedData();
    for(const auto &mask: edited.tipper_mask()) check(!mask[0] && !mask[1], "Inverted real/imaginary arrow masks conflict");
    check(edited.impedance_mask() == original.impedance_mask() && edited.phase_tensor_mask() == original.phase_tensor_mask(),
          "Inverted tipper arrows changed Z/PT masks");
  }
  // Z/PT edits must preserve both enabled and already-masked tipper samples.
  for(bool linked: {false, true}) for(auto type: {RealZxx, RealZxy, RealZyx, RealZyy,
      ImagZxx, ImagZxy, ImagZyx, ImagZyy, RhoZxx, RhoZxy, RhoZyx, RhoZyy,
      PhsZxx, PhsZxy, PhsZyx, PhsZyy, PTxx, PTxy, PTyx, PTyy}) {
    auto edited = original;
    edited.set_data_mask(RealTzx, .1, false);
    const auto before = edited.tipper_mask();
    for(bool enabled: {false, true}) {
      edited.set_data_mask(type, 1., enabled, linked);
      check(edited.tipper_mask() == before, "Point Z/PT edit changed tipper masks");
      edited.mask_type(type, enabled, linked);
      check(edited.tipper_mask() == before, "Whole-component Z/PT edit changed tipper masks");
    }
  }
  // Numerically check the dependency: perturb each real Z entry independently
  // in a nonsingular, non-2D tensor. Every PT entry must respond.
  for(auto zType: zTypes) {
    auto perturbed = original; double value, error; original.scalar_value(zType, 0, value, error);
    perturbed.set_data(1., {zType}, {value + 1e-6 * std::abs(value)}, {error}); perturbed.calculate_phase_tensor();
    for(auto ptType: ptTypes) {
      double before, after; original.scalar_value(ptType, 0, before, error); perturbed.scalar_value(ptType, 0, after, error);
      check(std::abs(after - before) > 1e-9, "Expected general 3D PT dependency on every real Z entry");
    }
  }
  auto station = original;
  station.set_data_mask(PhsZxy, 1., false, true);
  check(!station.impedance_mask()[1][0] && station.impedance_mask()[0][0] && station.impedance_mask()[2][0], "Linking recursively changed other Z components");
  for(unsigned c = 0; c < 4; ++c) check(!station.phase_tensor_mask()[c][0] && station.phase_tensor_mask()[c][1] && station.tipper_mask()[c][0], "Z mask affected wrong PT period or tippers");
  station.set_data_mask(RhoZxy, 1., true, true);
  station.set_data_mask(PTyx, 1., false, true);
  check(!station.phase_tensor_mask()[2][0] && station.phase_tensor_mask()[0][0] && station.phase_tensor_mask()[3][0], "PT linkage was recursive");
  for(unsigned c = 0; c < 4; ++c) check(!station.impedance_mask()[c][0] && station.impedance_mask()[c][1], "PT entry did not link to all real-Z dependencies");
  station.set_data_mask(PTyx, 1., true, true);
  station.mask_type(PhsZxy, false, true);
  for(unsigned c = 0; c < 4; ++c) for(unsigned f = 0; f < 2; ++f) check(!station.phase_tensor_mask()[c][f], "Whole-component mask did not link PT");
  station = original; station.set_data_mask(ImagTzx, 1., false, true);
  check(station.impedance_mask() == original.impedance_mask() && station.phase_tensor_mask() == original.phase_tensor_mask(), "Tipper masks must remain independent");
  station = original;
  QCustomPlot inversePlot; ApparentResistivityPlot inverse(&inversePlot); inverse.set_observed_data(station);
  inversePlot.graph(1)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  inversePlot.graph(2)->setSelection(QCPDataSelection(QCPDataRange(1, 2))); inverse.invMaskSelectedData();
  for(unsigned c = 0; c < 4; ++c) check(!station.phase_tensor_mask()[c][0] && !station.phase_tensor_mask()[c][1], "Inverted linked mask depends on component iteration order");
  check(station.tipper_mask() == original.tipper_mask(), "Inverted Z selection changed tipper masks");

  const auto project = directory + "/linked.mtd";
  { std::ofstream file(project.toStdString(), std::ios::binary); boost::archive::binary_oarchive archive(file); archive << survey << std::map<std::string, MTSurveyData>{}; }
  MainWindow window; file_action(window, "actionLoad_project", {project}); select_station(window, "linked");
  auto *setting = widget<QAction>(window, "actionLink_tensor_masks"); check(setting->isChecked(), "Linking should be enabled by default");
  auto &rho = *widget<QCustomPlot>(window, "plot11"), &phase = *widget<QCustomPlot>(window, "plot12"), &pt = *widget<QCustomPlot>(window, "plot22");
  auto &tipper = *widget<QCustomPlot>(window, "plot21");
  auto checkTippers = [&] {
    for(unsigned c: {0u, 2u})
      check(tipper.graph(c + 4)->data()->isEmpty() && tipper.graph(c)->data()->size() == 2,
            "Z/PT plot masking changed displayed tipper samples");
  };
  auto action = [&](QCustomPlot &plot, int graph, const QString &name) {
    plot.deselectAll(); plot.graph(graph)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
    QMetaObject::invokeMethod(&plot, "customContextMenuRequested", Qt::DirectConnection, Q_ARG(QPoint, QPoint(100, 100)));
    auto *menu = qobject_cast<QMenu *>(QApplication::activePopupWidget()); check(menu, "Mask menu did not open");
    QAction *chosen = nullptr; for(auto *item: menu->actions()) if(item->text() == name) chosen = item;
    check(chosen, "Mask action missing"); menu->hide(); chosen->trigger(); plot.deselectAll();
  };
  const auto periodRange = rho.xAxis->range(); action(phase, 1, "Mask");
  checkTippers();
  check(rho.graph(5)->data()->size() == 1 && rho.graph(4)->data()->size() == 1 && rho.graph(7)->data()->isEmpty(), "Phase masking missed its paired diagonal or changed the other pair");
  check(pt.graph(4)->data()->size() == 1, "Linked PT ellipses did not refresh immediately");
  check(rho.xAxis->range() == periodRange, "Masking changed zoom limits");
  action(phase, 5, "Unmask"); action(pt, 0, "Mask");
  checkTippers();
  for(unsigned c = 0; c < 4; ++c) check(rho.graph(c + 4)->data()->size() == 1 && phase.graph(c + 4)->data()->size() == 1, "PT masking did not refresh all related Z plots");
  check(pt.graph(4)->data()->size() == 1, "Ellipse selection did not mask the whole tensor"); action(pt, 4, "Unmask");
  setting->setChecked(false); action(rho, 1, "Mask");
  for(unsigned c = 0; c < 4; ++c) check(pt.graph(c + 4)->data()->isEmpty(), "Disabling linked masks still changed PT");
  check(rho.graph(4)->dataCount() == 1 && rho.graph(5)->dataCount() == 1, "Unlinked plot edit omitted paired diagonal");
  widget<QAction>(window, "actionUndo")->trigger();
  check(rho.graph(4)->data()->isEmpty() && rho.graph(5)->data()->isEmpty(), "One undo did not restore both impedance masks");
  widget<QAction>(window, "actionRedo")->trigger();
  check(rho.graph(4)->dataCount() == 1 && rho.graph(5)->dataCount() == 1, "Redo omitted paired diagonal");
  auto *stations = widget<QListWidget>(window, "stationList");
  bool componentEdited = false;
  QTimer::singleShot(0, [&] {
    auto *menu = qobject_cast<QMenu *>(QApplication::activePopupWidget());
    if(!menu) return;
    for(auto *action: menu->findChildren<QAction *>()) if(action->data().isValid() && action->data().toInt() == RealZyx) {
      action->trigger(); componentEdited = true; break;
    }
    menu->close();
  });
  QMetaObject::invokeMethod(stations, "customContextMenuRequested", Qt::DirectConnection,
                            Q_ARG(QPoint, stations->visualItemRect(stations->currentItem()).center()));
  check(componentEdited && rho.graph(6)->dataCount() == 2 && rho.graph(7)->dataCount() == 2 &&
        rho.graph(4)->dataCount() == 1 && pt.graph(4)->data()->isEmpty(), "Station component mask lost pairing or changed unrelated masks");
  checkTippers();
  widget<QAction>(window, "actionUndo")->trigger();
  check(rho.graph(6)->data()->isEmpty() && rho.graph(7)->data()->isEmpty() && rho.graph(4)->dataCount() == 1,
        "Station component undo did not restore the whole pair");
  widget<QAction>(window, "actionSave_project")->trigger();
  MainWindow restored; file_action(restored, "actionLoad_project", {project}); select_station(restored, "linked");
  auto *restoredRho = widget<QCustomPlot>(restored, "plot11");
  check(restoredRho->axisRectCount() == 2 && restoredRho->graph(4)->dataCount() == 1 && restoredRho->graph(5)->dataCount() == 1,
        "Project round trip lost TF layout or paired masks");
}

void period_maps(const QString &directory)
{
  MTMapData::PhaseTensor tensor;
  const double degrees = 180. / std::acos(-1.);
  check(MTMapData::phase_tensor({{2., 0., 0., 1.}}, tensor), "Diagonal tensor rejected");
  checkNear(tensor.phiMax, std::atan(2.) * degrees); checkNear(tensor.phiMin, 45.);
  checkNear(tensor.axisRatio, .5); checkNear(tensor.azimuth, 0.); checkNear(tensor.skew, 0.);
  check(MTMapData::phase_tensor({{1.75, std::sqrt(3.) / 4., std::sqrt(3.) / 4., 1.25}}, tensor), "Rotated tensor rejected");
  checkNear(tensor.azimuth, 30.); checkNear(tensor.axisRatio, .5);
  check(MTMapData::phase_tensor({{2., 1., 0., 1.}}, tensor), "Non-symmetric tensor rejected");
  // The ellipse's major direction is the eigenvector of Phi * transpose(Phi).
  checkNear(tensor.azimuth, .5 * std::atan2(2., 4.) * degrees);
  check(MTMapData::phase_tensor({{2., 0., 0., -1.}}, tensor), "Negative determinant tensor rejected");
  checkNear(tensor.phiMin, -45.); checkNear(tensor.axisRatio, .5);
  check(MTMapData::phase_tensor({{1., 0., 0., 1.}}, tensor), "Isotropic tensor rejected");
  checkNear(tensor.axisRatio, 1.);
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

  SurveyReport::Options reportOptions;
  reportOptions.overview = false; reportOptions.stationPages = false;
  reportOptions.phaseTensorMaps = true; reportOptions.inductionMaps = true;
  reportOptions.mapData = SurveyReport::Options::MapData::Both; reportOptions.mapPeriods = {1., 10.};
  reportOptions.mapSettings.imaginary = true; reportOptions.responseName = "GoFEM";
  check(SurveyReport::write_pdf(directory + "/report-map-values.pdf", *survey, &responses.begin()->second, reportOptions, [&](unsigned done, unsigned total) {
    check(total == 4, "Wrong combined map page count"); if(!done) return true;
    QCustomPlot *page = nullptr;
    for(auto *top: QApplication::topLevelWidgets()) if(top->objectName() == "periodMapWindow") page = top->findChild<QCustomPlot *>("periodMapPlot");
    check(page, "Missing report map");
    const unsigned tensors[]{2, 1, 1, 0}, vectors[]{2, 1, 1, 0};
    check(page->property("phaseTensorCount").toUInt() == tensors[done - 1] && page->property("realVectorCount").toUInt() == vectors[done - 1] &&
          page->property("imagVectorCount").toUInt() == vectors[done - 1], "Report map changed data, period or selected layers");
    return true;
  }), "Map data report failed");

  auto station = survey->get_station_data("A");
  const auto index = MTMapData::nearest_period(station.frequencies(), 1., 0.);
  std::array<double, 2> vector;
  check(MTMapData::phase_tensor(station, index, tensor), "Impedance-derived tensor rejected");
  checkNear(tensor.axisRatio, .5); checkNear(tensor.azimuth, 0.);
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
  checkNear(vector[0], 0.); checkNear(vector[1], 0.);
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
  checkNear(period->currentData().toDouble(), 1.);
  check(period->count() == 2 && dataset->count() == 2 && plot->property("phaseTensorCount").toUInt() == 2,
        "Default period map coverage is incorrect");
  auto *ellipse = widget<QCPCurve>(*plot, "phaseTensor_A");
  checkNear(ellipse->property("axisRatio").toDouble(), .5);
  auto *arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  const QPointF center = arrow->start->coords(), delta = arrow->end->coords() - center;
  const double size = widget<QDoubleSpinBox>(window, "mapArrowSize")->value();
  checkNear(delta.x(), -.1 * size); checkNear(delta.y(), -.2 * size);
  const auto majorPoint = ellipse->data()->constBegin();
  checkNear(majorPoint->key, center.x());
  checkNear(majorPoint->value - center.y(), .5 * widget<QDoubleSpinBox>(window, "mapEllipseSize")->value());
  check(!plot->findChild<QCPItemLine *>("inductionReal_C"), "Incomplete station vector displayed");
  widget<QComboBox>(window, "mapArrowConvention")->setCurrentIndex(1);
  arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  checkNear((arrow->end->coords() + delta - center).manhattanLength(), 0.);
  widget<QDoubleSpinBox>(window, "mapArrowSize")->setValue(size * 2.);
  arrow = widget<QCPItemLine>(*plot, "inductionReal_A");
  checkNear((arrow->end->coords() + 2. * delta - center).manhattanLength(), 0.);
  widget<QCheckBox>(window, "mapImagArrows")->setChecked(true);
  arrow = widget<QCPItemLine>(*plot, "inductionImag_A");
  checkNear(arrow->end->coords().x() - center.x(), -.1 * size * 2.);
  checkNear(arrow->end->coords().y() - center.y(), .05 * size * 2.);
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
  checkNear(period->currentData().toDouble(), 10.);
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
  check(!observedA.impedance_mask()[0][index] && !observedA.phase_tensor_mask()[0][index] && observedA.impedance_mask()[0][otherIndex],
        "Default linked PT map masking did not affect Z at the same period");
  unmask->click();
  window.setLinkTensorMasks(false); mask->click();
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
  widget<QCheckBox>(window, "mapRealArrows")->setChecked(true);
  group->setCurrentIndex(2); mask->click();
  check(plot->findChild<QCPItemLine *>("inductionReal_A") && plot->findChild<QCPItemLine *>("inductionImag_A"),
        "Z/PT map masking removed independent induction vectors");
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
  checkNear((arrow->end->coords() - arrow->start->coords() - expected).manhattanLength(), 0.);
  check(scale->dataRange() == QCPRange(10., 60.), "Refreshing maps reset color limits");
  window.resize(1300, 800); QApplication::processEvents(); plot->replot();
  checkNear(plot->xAxis->range().size() / plot->axisRect()->width(), plot->yAxis->range().size() / plot->axisRect()->height());
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
  for(const auto &point: *widget<QCPCurve>(*plot, "phaseTensor_A")->data())
    check(plot->xAxis->range().contains(point.key) && plot->yAxis->range().contains(point.value),
          "Single-station fit clipped its ellipse");
  period->setCurrentIndex(0);
  checkNear(period->currentData().toDouble(), 1. / 1.0005);
  widget<QCheckBox>(window, "mapSelectStations")->setChecked(true);
  const auto closeCenter = widget<QCPItemLine>(*plot, "inductionReal_A")->start->pixelPosition();
  map_gesture(plot, closeCenter, closeCenter);
  group->setCurrentIndex(0); mask->click();
  const int exactIndex = MTMapData::nearest_period(closeStation.frequencies(), 1. / 1.0005, 0.);
  const int neighborIndex = MTMapData::nearest_period(closeStation.frequencies(), 1., 0.);
  check(!closeStation.tipper_mask()[0][exactIndex] && closeStation.tipper_mask()[0][neighborIndex],
        "Map masked a nearby frequency instead of the selected frequency");

  // A finite but extreme vector must remain inspectable without collapsing
  // the entire survey to a dot in the interactive map or report.
  write(observedPath, rows("West", "1") + rows("Center", "1") + rows("East", "1"));
  auto wide = std::make_shared<MTSurveyData>();
  wide->load_from_native_responses(observedPath.toStdString());
  wide->get_station_data("West").set_position({{42., -124., 100.}});
  wide->get_station_data("Center").set_position({{35., -96., 100.}});
  wide->get_station_data("East").set_position({{45., -68., 100.}});
  wide->get_station_data("Center").set_data(1., {RealTzx, RealTzy}, {2e6, -1e6}, {.01, .01});
  PeriodMapWindow wideMap(nullptr); wideMap.setData(wide, {}); wideMap.show(); QApplication::processEvents();
  auto *widePlot = widget<QCustomPlot>(wideMap, "periodMapPlot");
  const auto checkExtent = [](QCustomPlot *map) {
    check(map->xAxis->range().size() < 20000. && map->yAxis->range().size() < 20000.,
          "Outlier glyph expanded the map beyond the survey");
    for(const auto &point: *widget<QCPGraph>(*map, "mapStations")->data())
      check(map->xAxis->range().contains(point.key) && map->yAxis->range().contains(point.value),
            "Fitted map clipped a station");
  };
  checkExtent(widePlot);
  auto *outlier = widget<QCPItemLine>(*widePlot, "inductionReal_Center");
  check(outlier->property("north").toDouble() == 2e6, "Map silently changed the outlier value");
  widget<QPushButton>(wideMap, "mapFitView")->click(); checkExtent(widePlot);
  SurveyReport::Options wideReport; wideReport.overview = false; wideReport.stationPages = false;
  wideReport.phaseTensorMaps = true; wideReport.inductionMaps = true; wideReport.mapPeriods = {1.};
  check(SurveyReport::write_pdf(directory + "/wide-map.pdf", *wide, nullptr, wideReport, [&](unsigned done, unsigned) {
    if(done) for(auto *top: QApplication::topLevelWidgets())
      if(auto *map = top->findChild<QCustomPlot *>("periodMapPlot")) checkExtent(map);
    return true;
  }), "Continental map report failed");
  widget<QCheckBox>(wideMap, "mapSelectStations")->setChecked(true);
  const auto outlierCenter = outlier->start->pixelPosition();
  map_gesture(widePlot, outlierCenter, outlierCenter);
  widget<QPushButton>(wideMap, "mapMask")->click();
  check(!wide->get_station_data("Center").tipper_mask()[0][0] &&
        wide->get_station_data("West").tipper_mask()[0][0], "Outlier could not be masked from its map location");

  // Scalar maps use local complex impedances, not means of rho or phase.
  auto &centerStation = wide->get_station_data("Center");
  centerStation.set_data(1., {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy},
                        {.2, .1, 2., 3., -4., -2., .7, -.2}, std::vector<double>(8, .01));
  wideMap.setData(wide, {});
  auto *colorBy = widget<QComboBox>(wideMap, "mapColorBy");
  auto *wideScale = dynamic_cast<QCPColorScale *>(widePlot->plotLayout()->element(1, 1));
  using Complex = std::complex<double>;
  const Complex xx(.2, .1), xy(2., 3.), yx(-4., -2.), yy(.7, -.2);
  const std::array<Complex, 4> expectedZ{{xy, -yx, (xy - yx) / 2., std::sqrt((xx*xx + xy*xy + yx*yx + yy*yy) / 2.)}};
  for(int mode = 3; mode < 11; ++mode) {
    colorBy->setCurrentIndex(mode);
    const bool rho = mode < 7;
    const auto z = expectedZ[(mode - 3) % 4];
    const double expected = rho ? std::norm(z) / (8e-7 * std::pow(std::acos(-1.), 2)) : std::arg(z) * degrees;
    auto *marker = widget<QCPItemTracer>(*widePlot, "impedanceMarker_Center");
    checkNear(marker->property("value").toDouble(), expected);
    check(wideScale->visible() && wideScale->dataScaleType() == (rho ? QCPAxis::stLogarithmic : QCPAxis::stLinear),
          "Scalar map used the wrong color scale");
    check(widget<QCPCurve>(*widePlot, "phaseTensor_Center")->brush().style() == Qt::NoBrush,
          "Scalar map left a conflicting filled tensor color");
  }
  colorBy->setCurrentIndex(7); // XY phase remains visible without YX or diagonals.
  centerStation.set_data_mask(RealZyx, 1., false, false); wideMap.setData(wide, {});
  check(widePlot->findChild<QCPItemTracer *>("impedanceMarker_Center"), "XY map required unrelated components");
  for(int mode: {8, 9, 10}) {
    colorBy->setCurrentIndex(mode);
    check(!widePlot->findChild<QCPItemTracer *>("impedanceMarker_Center"), "Scalar map included a masked required component");
  }
  centerStation.set_data_mask(RealZyx, 1., true, false);
  centerStation.set_data_mask(RealZxx, 1., false, false); wideMap.setData(wide, {});
  check(!widePlot->findChild<QCPItemTracer *>("impedanceMarker_Center"), "SSQ replaced a masked diagonal with zero");
  colorBy->setCurrentIndex(9);
  check(widePlot->findChild<QCPItemTracer *>("impedanceMarker_Center"), "Berdichevsky map required diagonal data");
  centerStation.set_data_mask(RealZxx, 1., true, false);
  widget<QCheckBox>(wideMap, "mapPhaseTensors")->setChecked(false);
  colorBy->setCurrentIndex(3);
  check(wideScale->visible() && !widePlot->findChild<QCPCurve *>("phaseTensor_Center"),
        "Scalar colors depended on ellipse visibility");
  widget<QDoubleSpinBox>(wideMap, "mapArrowWidth")->setValue(.4);
  widget<QDoubleSpinBox>(wideMap, "mapColorMin")->setValue(1000.);
  widget<QDoubleSpinBox>(wideMap, "mapColorMax")->setValue(1e7);
  const auto settings = wideMap.displaySettings();
  wideReport.useMapSettings = true; wideReport.mapSettings = settings;
  check(SurveyReport::write_pdf(directory + "/scalar-map.pdf", *wide, nullptr, wideReport, [&](unsigned done, unsigned) {
    if(done) for(auto *top: QApplication::topLevelWidgets()) if(top != &wideMap && top != &window)
      if(auto *map = top->findChild<QCustomPlot *>("periodMapPlot")) {
        const auto *scale = dynamic_cast<QCPColorScale *>(map->plotLayout()->element(1, 1));
        check(scale && scale->dataScaleType() == QCPAxis::stLogarithmic && scale->dataRange() == QCPRange(1000., 1e7),
              "Report lost scalar map range or logarithmic scale");
        checkNear(widget<QCPItemLine>(*map, "inductionReal_West")->pen().widthF(), .4);
        check(map->property("impedanceMarkerCount").toUInt() == 3, "Report lost scalar markers");
      }
    return true;
  }), "Scalar map report failed");
  colorBy->setCurrentIndex(0);
  check(!widePlot->findChild<QCPItemTracer *>("impedanceMarker_Center") && !wideScale->visible(),
        "Switching to hidden tensors retained scalar markers or legend");
}

void station_selection_edits(const QString &directory)
{
  group_mask_tests(directory);
  auto survey = std::make_shared<MTSurveyData>();
  const auto data = directory + "/delete-stations.gofem", project = directory + "/delete-stations.mtd";
  std::string rows;
  for(const auto *name : {"A", "B", "C"})
    rows += std::string("RealZxy 1 Plane_wave ") + name + " .002 .0001\nImagZxy 1 Plane_wave " + name + " .003 .0001\n";
  write(data, rows); survey->load_from_gofem(data.toStdString());
  for(const auto &name : survey->get_stations_names()) survey->get_station_data(name).set_position({{55., 9., 100.}});
  {
    std::ofstream file(project.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(file); archive << survey << std::map<std::string, MTSurveyData>{};
  }
  MainWindow window;
  file_action(window, "actionLoad_project", {project});
  auto *list = widget<QListWidget>(window, "stationList");
  auto *remove = widget<QAction>(window, "actionDelete_stations");
  auto *undo = widget<QAction>(window, "actionUndo");
  auto *redo = widget<QAction>(window, "actionRedo");
  auto *map = widget<QCustomPlot>(window, "mapPlot");
  list->setCurrentRow(1); list->clearSelection();
  remove->trigger();
  check(list->count() == 3 && !undo->isEnabled(), "Empty selection deleted the current station");
  list->item(0)->setSelected(true); list->item(2)->setSelected(true);
  remove->trigger();
  check(list->count() == 1 && list->item(0)->text() == "B" && map->graph(0)->dataCount() == 1,
        "Delete did not remove exactly the selected stations and update the map");
  undo->trigger();
  check(list->count() == 3 && map->graph(0)->dataCount() == 3 && !undo->isEnabled(),
        "Batch deletion was not restored in one undo step");
  redo->trigger();
  check(list->count() == 1 && list->item(0)->text() == "B", "Redo did not repeat batch deletion");
  list->selectAll(); remove->trigger();
  check(list->count() == 0 && map->graph(0)->dataCount() == 0 &&
        !widget<QAction>(window, "actionExport_native_MT")->isEnabled(), "Deleting the last station left stale survey state");
  remove->trigger(); undo->trigger();
  check(list->count() == 1 && list->item(0)->text() == "B", "Empty deletion added an undo step");
  undo->trigger(); list->selectAll(); remove->trigger();
  check(list->count() == 0, "Delete all left stations behind");
  undo->trigger(); check(list->count() == 3, "Undo did not restore all deleted stations");
  list->clearSelection();
  auto *groupMask = widget<QAction>(window, "actionGroup_masking");
  check(!groupMask->isEnabled(), "Group masking enabled without selected stations");
  list->item(0)->setSelected(true); list->item(2)->setSelected(true);
  QTimer::singleShot(0, [&] {
    auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
    check(dialog && dialog->objectName() == "groupMaskDialog", "Group masking action did not open dialog");
    check(widget<QCheckBox>(*dialog, "groupMaskLink")->isChecked(), "Group dialog ignored the main linkage setting");
    widget<QCheckBox>(*dialog, "groupMaskLink")->setChecked(false);
    widget<QCheckBox>(*dialog, "groupMaskZxy")->setChecked(true);
    check(widget<QTableWidget>(*dialog, "groupMaskPreview")->rowCount() == 2, "Group dialog ignored selected station scope");
    const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
    if(!artifacts.isEmpty()) { QApplication::processEvents(); dialog->grab().save(artifacts + "/group-mask.png"); }
    widget<QPushButton>(*dialog, "groupMaskApply")->click();
  });
  groupMask->trigger();
  check(list->selectedItems().size() == 2 && widget<QAction>(window, "actionLink_tensor_masks")->isChecked(),
        "Group masking changed selection or global linkage");
  const auto masked = [&](MainWindow &view, const QString &name) {
    select_station(view, name); return widget<QCustomPlot>(view, "plot11")->graph(5)->dataCount();
  };
  check(masked(window, "A") == 1 && masked(window, "C") == 1 && masked(window, "B") == 0, "Group edit did not refresh the correct station plots");
  undo->trigger();
  check(!undo->isEnabled() && masked(window, "A") == 0 && masked(window, "C") == 0, "Group masking was not one Undo step");
  redo->trigger();
  check(masked(window, "A") == 1 && masked(window, "C") == 1, "Group mask redo failed");
  widget<QAction>(window, "actionSave_project")->trigger();
  MainWindow restored; file_action(restored, "actionLoad_project", {project});
  check(masked(restored, "A") == 1 && masked(restored, "C") == 1 && masked(restored, "B") == 0, "Project lost group masking");
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
  checkNear(value(*widget<QCustomPlot>(*copy, "plot12"), 1, 2.), std::atan2(3., 2.) * 180. / std::acos(-1.));
  check(widget<QCustomPlot>(window, "plot12")->graph(1)->data()->findBegin(2.)->key != 2., "Resampling modified the original window");
  file_action(*copy, "actionSave_project", {saved});
  copy->close(); QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
  MainWindow restored;
  file_action(restored, "actionLoad_project", {saved}); select_station(restored, "A");
  checkNear(value(*widget<QCustomPlot>(restored, "plot12"), 1, 2.), std::atan2(3., 2.) * 180. / std::acos(-1.));
  widget<QAction>(restored, "actionPeriod_layout")->trigger();
  auto *restoredLayout = widget<QDialog>(restored, "periodLayoutWindow");
  check(!widget<QPushButton>(*restoredLayout, "layoutOpenSource")->isHidden(), "Saved project lost source access");
  widget<QPushButton>(*restoredLayout, "layoutOpenSource")->click(); QApplication::processEvents();
  copy = nullptr;
  for(auto *top: QApplication::topLevelWidgets())
    if(top->objectName() == "resampledSurveyWindow") copy = qobject_cast<MainWindow *>(top);
  check(copy, "Could not reopen the saved source survey");
  select_station(*copy, "A");
  checkNear(value(*widget<QCustomPlot>(*copy, "plot12"), 1, 100.), std::atan2(9., 8.) * 180. / std::acos(-1.));
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
struct LegacyOptionsV6 {
  bool phaseWrap = true, showStationNames = false, tipperArrows = false;
  std::vector<LegacyAxisV5> axes = std::vector<LegacyAxisV5>(4);
  std::array<std::array<bool, 4>, 4> visible{{{{false, true, true, false}}, {{true, false, false, true}},
                                            {{true, false, true, false}}, {{false, true, false, true}}}};
  template<class Archive> void serialize(Archive &ar, const unsigned int) {
    ar & phaseWrap & showStationNames & tipperArrows & axes & visible;
  }
};

struct LegacyOptionsV7 : LegacyOptionsV6 {
  bool linkTensorMasks = true;
  template<class Archive> void serialize(Archive &ar, const unsigned int version) {
    LegacyOptionsV6::serialize(ar,version); ar & linkTensorMasks;
  }
};
struct LegacyOptionsV8 : LegacyOptionsV7 {
  bool rhoPlus = true;
  double obsoleteFloor = 0.;
  unsigned density = 16;
  template<class Archive> void serialize(Archive &ar, const unsigned int version) {
    LegacyOptionsV7::serialize(ar, version); ar & rhoPlus & obsoleteFloor & density;
  }
};
void tensor_strip()
{
  MTStationData station;
  station.set_size(2, true); station.set_frequencies({1., .1});
  for(double f: station.frequencies()) {
    station.set_data(f, {PTxx, PTxy, PTyx, PTyy}, {2., .2, .2, 1.}, {.01, .01, .01, .01});
    station.set_data(f, {RealTzx, RealTzy, ImagTzx, ImagTzy}, {.1, .2, .3, .4}, {.01, .01, .01, .01});
  }
  QCustomPlot plot; plot.resize(800, 160); plot.show();
  PhaseTensorPlot handler(&plot); handler.set_link_tensor_masks(false);
  handler.set_ellipse_mode(true); handler.set_observed_data(station);
  check(plot.graph(0)->data()->size() == 2 && plot.graph(4)->data()->isEmpty(), "Explicit PT ellipse anchors missing");
  const auto before = plot.toPixmap().toImage();
  plot.graph(0)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  handler.maskSelectedData();
  check(plot.graph(4)->data()->size() == 1, "Masked ellipse disappeared instead of becoming selectable grey data");
  for(const auto &mask: station.phase_tensor_mask()) check(mask[0] && !mask[1], "Ellipse masking missed a PT component");
  for(const auto &mask: station.tipper_mask()) check(mask[0] && mask[1], "Ellipse masking affected tippers");
  plot.graph(0)->setSelection(QCPDataSelection());
  plot.graph(4)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  handler.unmaskSelectedData();
  handler.set_predicted_data(station);
  check(plot.toPixmap().toImage() != before, "Response ellipse outlines did not render");
  handler.clear_predicted_data(); handler.clear_observed_data();
  check(plot.graph(0)->data()->isEmpty() && plot.graph(4)->data()->isEmpty(), "Clearing a station left ellipse anchors");
  handler.set_ellipse_mode(false); handler.set_observed_data(station);
  checkNear(value(plot, 0, 1.), 2.); checkNear(value(plot, 3, 1.), 1.);
  check(plot.legend->itemCount() == 4, "Scalar PT controls were not restored");
}

void tf_view_controls(MainWindow &window, const QString &directory, const QString &project)
{
  select_station(window, "GUIDE");
  auto *heading = widget<QLabel>(window, "tfStationName");
  check(heading->text() == "GUIDE" && heading->textFormat() == Qt::PlainText, "TF heading is not the current station");
  auto *stations = widget<QListWidget>(window, "stationList");
  check(stations->currentItem()->data(StationListDelegate::PeriodCount).toUInt() == 12, "Missing station period count");
  widget<QAction>(window, "actionPlot_axis_ranges")->trigger();
  check(!QApplication::activeModalWidget() && widget<QTabBar>(window, "ribbonTabs")->currentIndex() == 6,
        "View controls must be in the ribbon");
  auto *rho = widget<QCustomPlot>(window, "plot11"), *phase = widget<QCustomPlot>(window, "plot12");
  auto *tipper = widget<QCustomPlot>(window, "plot21"), *tensor = widget<QCustomPlot>(window, "plot22");
  const auto setRange = [&](int index, double minimum, double maximum) {
    widget<QComboBox>(window, "viewRangeAxis")->setCurrentIndex(index);
    widget<QLineEdit>(window, "viewRangeMin")->setText(QString::number(minimum));
    widget<QLineEdit>(window, "viewRangeMax")->setText(QString::number(maximum));
    widget<QPushButton>(window, "viewRangeApply")->click();
  };
  setRange(0, .03, 30.); setRange(1, 10., 500.); setRange(3, .1, 20.);
  setRange(2, -150., 100.); setRange(4, -170., 170.);
  check(rho->yAxis->range() == QCPRange(10, 500) && rho->graph(0)->valueAxis()->range() == QCPRange(.1, 20),
        "TF columns did not retain independent ranges");
  for(auto *plot: {rho, phase, tipper, tensor}) check(plot->xAxis->range() == QCPRange(.03, 30), "Ribbon period range is not shared");
  setRange(0, -1, 30);
  check(rho->xAxis->range() == QCPRange(.03, 30) && !widget<QLineEdit>(window, "viewRangeMin")->property("validInput").toBool(),
        "Invalid logarithmic range changed the view");
  setRange(0, 30, 1);
  check(rho->xAxis->range() == QCPRange(.03, 30), "Reversed limits changed the range");
  setRange(0, .03, 30);
  rho->graph(1)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  const double observed = rho->graph(1)->data()->constBegin()->value;
  const auto points = rho->graph(1)->data()->size();
  const auto tensorBefore = tensor->toPixmap().toImage();
  const auto setColor = [&](int component, const QString &color) {
    widget<QComboBox>(window, "viewColorTarget")->setCurrentIndex(component);
    auto *editor = widget<QComboBox>(window, "viewColor"); editor->setEditText(color);
    QMetaObject::invokeMethod(editor->lineEdit(), "editingFinished");
  };
  setColor(1, "#15803d"); setColor(4, "#dc2626"); setColor(5, "#7c3aed"); setColor(7, "#dc2626");
  setColor(1, "invalid-color");
  check(rho->graph(1)->scatterStyle().pen().color() == QColor("#15803d"), "Invalid color replaced the applied color");
  setColor(1, "#15803d");
  widget<QDoubleSpinBox>(window, "viewMarkerSize")->setValue(9.);
  widget<QDoubleSpinBox>(window, "viewLineWidth")->setValue(2.);
  widget<QDoubleSpinBox>(window, "viewArrowScale")->setValue(120.);
  widget<QDoubleSpinBox>(window, "viewEllipseSize")->setValue(40.);
  widget<QComboBox>(window, "viewColorMap")->setCurrentText("Thermal");
  widget<QCheckBox>(window, "viewColorMapReverse")->setChecked(true);
  widget<QDoubleSpinBox>(window, "viewColorMin")->setValue(10.);
  widget<QDoubleSpinBox>(window, "viewColorMax")->setValue(70.);
  widget<QAction>(window, "actionShow_error_bars")->setChecked(false);
  QApplication::processEvents();
  check(rho->graph(1)->selection().dataPointCount() == 1 && rho->graph(1)->data()->size() == points &&
        rho->graph(1)->data()->constBegin()->value == observed && rho->graph(5)->data()->isEmpty(),
        "Appearance changes altered selection, observations or masks");
  check(rho->xAxis->range() == QCPRange(.03, 30) && rho->yAxis->range() == QCPRange(10, 500), "Appearance changed zoom");
  for(auto *plot: {rho, phase}) {
    check(plot->graph(1)->scatterStyle().size() == 9. && plot->graph(9)->pen().color() == QColor("#15803d") &&
          plot->graph(9)->pen().widthF() == 2. && plot->graph(5)->scatterStyle().pen().color() == QColor("#dc2626"),
          "Colors/sizes did not reach observations, masked points and responses");
    for(auto *bars: plot->findChildren<QCPErrorBars *>()) if(bars->dataPlottable() == plot->graph(1))
      check(!bars->visible() && bars->pen().color() == QColor("#15803d"), "Error bars ignored view settings");
  }
  check(widget<QCPCurve>(*rho, "rhoPlusXY")->pen().color() == QColor("#15803d"), "Rho+ guide color is stale");
  check(tipper->graph(0)->scatterStyle().pen().color() == QColor("#7c3aed"), "Tipper marker color is stale");
  bool reference = false;
  for(auto *arrow: tipper->findChildren<QCPItemLine *>())
    if(arrow->objectName() == "tipperReferenceArrow") {
      reference = true;
      checkNear(QLineF(arrow->start->pixelPosition(), arrow->end->pixelPosition()).length(), 120.);
      check(arrow->pen().color() == QColor("#7c3aed") && arrow->pen().widthF() == 2., "Tipper reference did not follow appearance");
    }
  check(reference, "Tipper reference missing after restyling");
  auto *scale = tensor->findChild<QCPColorScale *>();
  check(scale && scale->dataRange() == QCPRange(10, 70) && scale->gradient() == PlotColorMaps::gradient("Thermal").inverted(),
        "Tensor color scale does not match ribbon");
  check(tensor->toPixmap().toImage() != tensorBefore, "Ellipse appearance did not repaint");
  select_station(window, "OUTLIER"); select_station(window, "GUIDE");
  check(heading->text() == "GUIDE" && rho->xAxis->range() == QCPRange(.03, 30) &&
        rho->graph(0)->valueAxis()->range() == QCPRange(.1, 20), "Station navigation lost heading or fixed ranges");
  // Masking remains connected after recoloring and resizing; the paired diagonal follows XY.
  for(int i = 0; i < 8; ++i) rho->graph(i)->setSelection(QCPDataSelection());
  rho->graph(1)->setSelection(QCPDataSelection(QCPDataRange(0, 1)));
  widget<QMenu>(*rho, "")->actions().front()->trigger();
  check(rho->graph(4)->data()->size() == 1 && rho->graph(5)->data()->size() == 1, "Restyling broke paired masking");
  widget<QAction>(window, "actionUndo")->trigger();
  check(rho->graph(4)->data()->isEmpty() && rho->graph(5)->data()->isEmpty(), "Restyled masks cannot be undone");
  widget<QAction>(window, "actionSave_project")->trigger();
  {
    MainWindow restored; file_action(restored, "actionLoad_project", {project}); select_station(restored, "GUIDE");
    auto *plot = widget<QCustomPlot>(restored, "plot11");
    check(plot->xAxis->range() == QCPRange(.03, 30) && plot->graph(0)->valueAxis()->range() == QCPRange(.1, 20) &&
          plot->graph(1)->scatterStyle().size() == 9. && plot->graph(9)->pen().color() == QColor("#15803d"), "Project lost TF appearance/ranges");
    check(!widget<QAction>(restored, "actionShow_error_bars")->isChecked() &&
          widget<QComboBox>(restored, "viewColorMap")->currentText() == "Thermal" &&
          widget<QCheckBox>(restored, "viewColorMapReverse")->isChecked() &&
          widget<QDoubleSpinBox>(restored, "viewArrowScale")->value() == 120., "Project lost strip settings");
  }
  const auto imagePath = directory + "/tf-heading.png";
  file_action(window, "actionSave_as_PDF", {imagePath});
  const QImage exported(imagePath);
  check(!exported.isNull() && exported.height() > widget<QSplitter>(window, "splitter_4")->height(), "Plot image omitted station heading");
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if(!artifacts.isEmpty()) {
    window.showNormal(); window.resize(1400, 1000); QApplication::processEvents();
    check(window.grab().save(artifacts + "/view-ribbon.png"), "Could not save ribbon preview");
    window.resize(1000, 800); QApplication::processEvents();
    auto *referenceArrow = widget<QCPItemLine>(*tipper, "tipperReferenceArrow");
    check(tipper->axisRect()->rect().contains(referenceArrow->start->pixelPosition().toPoint()) &&
          tipper->axisRect()->rect().contains(referenceArrow->end->pixelPosition().toPoint()), "Narrow view clipped the reference arrow");
    check(window.grab().save(artifacts + "/view-ribbon-narrow.png"), "Could not save narrow ribbon preview");
  }
  widget<QAction>(window, "actionDecimate")->trigger();
  check(stations->currentItem()->data(StationListDelegate::PeriodCount).toUInt() == 6, "Decimation left stale period counts");
  widget<QAction>(window, "actionUndo")->trigger();
  check(stations->currentItem()->data(StationListDelegate::PeriodCount).toUInt() == 12, "Undo did not restore period counts");
  QTimer nameTimer;
  QObject::connect(&nameTimer, &QTimer::timeout, [&] {
    if(auto *dialog = qobject_cast<QInputDialog *>(QApplication::activeModalWidget())) {
      dialog->setTextValue("<renamed>"); nameTimer.stop(); dialog->accept();
    }
  });
  nameTimer.start(0); QMetaObject::invokeMethod(&window, "renameStation");
  check(heading->text() == "<renamed>", "Rename did not refresh the plain-text heading");
  widget<QAction>(window, "actionUndo")->trigger();
  check(heading->text() == "GUIDE", "Undo rename left a stale station heading");
  widget<QAction>(window, "actionFit_TF_axes")->trigger();
  check(rho->xAxis->range() != QCPRange(.03, 30) && widget<QCheckBox>(window, "viewRangeAuto")->isChecked(), "Fit axes did not restore automatic limits");
  widget<QAction>(window, "actionReset_TF_appearance")->trigger();
  check(rho->graph(1)->scatterStyle().size() == 5 && rho->graph(9)->pen().color() == PlotColors::componentColors()[1], "Reset appearance failed");
  stations->setCurrentRow(-1);
  check(heading->text() == "No station selected", "Empty selection retained previous station heading");
}

void rho_plus_overlay(const QString &directory)
{
  const auto input=directory+"/rho-guide.gofem", project=directory+"/rho-guide.mtd";
  std::ostringstream text; text << std::setprecision(17);
  const double pi=std::acos(-1.);
  for(const auto &name:{std::string("GUIDE"),std::string("OUTLIER")}) for(int i=0;i<12;++i) {
    const double f=std::pow(10.,-2.+4.*i/11), z=std::sqrt(pi*f*4e-7*pi*100.);
    for(const auto &part:{"Real","Imag"}) for(const auto &component:{"xy","yx"}) {
      double value=component==std::string("xy")?z:-z;
      if(name=="OUTLIER" && i==5 && component==std::string("xy") && part==std::string("Imag")) value=-value;
      text<<part<<'Z'<<component<<' '<<f<<" Plane_wave "<<name<<' '<<value<<' '<<.02*z<<'\n';
    }
    text << "RealZxx " << f << " Plane_wave " << name << ' ' << .1*z << ' ' << .01*z << '\n'
         << "ImagZxx " << f << " Plane_wave " << name << ' ' << .2*z << ' ' << .01*z << '\n'
         << "RealZyy " << f << " Plane_wave " << name << ' ' << .15*z << ' ' << .01*z << '\n'
         << "ImagZyy " << f << " Plane_wave " << name << ' ' << -.1*z << ' ' << .01*z << '\n';
    for(const auto &type: {"RealTzx", "RealTzy", "ImagTzx", "ImagTzy"})
      text << type << ' ' << f << " Plane_wave " << name << ' ' << .1*std::sin(i*.3) << " .01\n";
  }
  write(input,text.str());
  auto survey=std::make_shared<MTSurveyData>(); survey->load_from_gofem(input.toStdString());
  survey->get_station_data("GUIDE").set_position({55., 12., 100.});
  survey->get_station_data("OUTLIER").set_position({55.1, 12.1, 100.});
  {
    std::ofstream file(project.toStdString(),std::ios::binary); boost::archive::binary_oarchive archive(file);
    const std::uint32_t magic=0x45444954, version=7;
    LegacyOptionsV7 options; for(auto &visible:options.visible) visible.fill(true);
    archive << magic << version << survey << std::map<std::string,MTSurveyData>{} << options;
  }
  auto settle=[] { QEventLoop loop; QTimer::singleShot(300,&loop,&QEventLoop::quit); loop.exec(); };
  MainWindow window; file_action(window,"actionLoad_project",{project}); select_station(window,"GUIDE");
  auto *toggle=widget<QAction>(window,"actionShow_rho_plus"); check(!toggle->isChecked(),"Legacy project enabled rho+ guides");
  auto *rho=widget<QCustomPlot>(window,"plot11"), *phase=widget<QCustomPlot>(window,"plot12");
  auto *xy=widget<QCPCurve>(*rho,"rhoPlusXY"), *yxPhase=widget<QCPCurve>(*phase,"rhoPlusYX");
  rho->xAxis->setRange(.02,20.); rho->yAxis->setRange(10.,1000.);
  const auto xRange=rho->xAxis->range(), yRange=rho->yAxis->range();
  toggle->setChecked(true); settle();
  check(xy->visible() && xy->data()->size()==201,"Rho+ guide missing");
  check(xy->pen().color()==PlotColors::componentColors()[1] && xy->pen().widthF()<1. && xy->pen().style()==Qt::DashLine,"Guide style differs from component");
  check(xy->selectable()==QCP::stNone && rho->legend->itemCount()==2 && rho->graphCount()==12,"Guide affects masking or legend layout");
  check(rho->xAxis->range()==xRange && rho->yAxis->range()==yRange,"Fitting changed zoom");
  check(std::abs(xy->data()->constBegin()->value-100.)<.01 &&
        std::abs(yxPhase->data()->constBegin()->value-45.)<.01,"Incorrect half-space guide");
  const double originalStart=xy->data()->constBegin()->key;
  file_action(window,"actionLoad_responses",{input}); settle();
  check(!rho->graph(9)->data()->isEmpty() && xy->visible(),"Loaded response and rho+ cannot coexist");
  check(!window.findChild<QAction *>("actionTransfer_function_view") && !window.findChild<QAction *>("actionTipper_arrows"),
        "The main view still offers alternative layouts");
  check(rho->axisRectCount() == 2 && phase->axisRectCount() == 2, "TF columns missing");
  check(rho->graph(0)->keyAxis() != rho->xAxis && rho->graph(1)->keyAxis() == rho->xAxis,
        "TF components assigned to wrong columns");
  check(rho->graph(8)->keyAxis() == rho->graph(0)->keyAxis() && xy->visible(), "TF responses or guides missing");
  {
    QSettings().setValue("plots/transferFunctionView", false);
    MainWindow remembered;
    check(widget<QCustomPlot>(remembered, "plot11")->axisRectCount() == 2 &&
          widget<QCustomPlot>(remembered, "plot21")->legend->itemCount() == 2 &&
          widget<QCustomPlot>(remembered, "plot22")->legend->itemCount() == 1,
          "An old preference disabled the default TF layout");
  }
  auto *rightAxis = rho->graph(0)->keyAxis();
  rightAxis->setRange(.1, 10.); settle();
  check(phase->xAxis->range() == rightAxis->range() &&
        widget<QCustomPlot>(window, "plot21")->xAxis->range() == rightAxis->range(), "TF period axes are not linked");
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if(!artifacts.isEmpty()) {
    window.showNormal(); window.resize(1400, 1000); QMetaObject::invokeMethod(&window, "updatePlots"); settle();
    check(window.grab().save(artifacts + "/transfer-function-view.png"), "Could not save TF view preview");
  }
  auto *menu=widget<QMenu>(*rho,"");
  rho->graph(0)->setSelection(QCPDataSelection(QCPDataRange(0,1)));
  menu->actions().front()->trigger(); settle();
  check(rho->graph(4)->data()->size() == 1 && rho->graph(5)->data()->isEmpty() &&
        rho->graph(6)->data()->isEmpty() && rho->graph(7)->data()->isEmpty(),
        "Masking a diagonal in TF view affected another impedance component");
  rho->graph(0)->setSelection(QCPDataSelection());
  rho->graph(4)->setSelection(QCPDataSelection(QCPDataRange(0,1)));
  menu->actions()[2]->trigger(); settle();
  rho->graph(4)->setSelection(QCPDataSelection());
  rho->graph(1)->setSelection(QCPDataSelection(QCPDataRange(0,1)));
  menu->actions().front()->trigger(); settle();
  check(xy->data()->constBegin()->key>originalStart && rho->graph(5)->data()->size()==1,"Masking did not refit the guide");
  check(!rho->graph(9)->data()->isEmpty(),"Masking cleared loaded responses");
  rho->graph(5)->setSelection(QCPDataSelection(QCPDataRange(0,1)));
  menu->actions()[2]->trigger(); settle();
  checkNear(xy->data()->constBegin()->key,originalStart);
  component_box(*rho,1)->setChecked(false);
  check(!xy->visible() && widget<QCPCurve>(*phase,"rhoPlusXY")->visible(),"Guide ignored independent component visibility");
  component_box(*rho,1)->setChecked(true);
  widget<QAction>(window,"actionPhase_wrap")->setChecked(false);
  check(std::abs(yxPhase->data()->constBegin()->value+135.)<.01,"Guide ignored phase wrapping");
  widget<QAction>(window,"actionImpedance_validity")->trigger();
  auto *validity=widget<QDialog>(window,"impedanceValidityWindow");
  check(!validity->findChild<QDoubleSpinBox *>("validityErrorFloor"),"Rho+ still has an independent error floor");
  widget<QSpinBox>(*validity,"validityGridDensity")->setValue(16); settle();
  validity->hide();
  widget<QAction>(window,"actionSave_project")->trigger();
  MainWindow restored; file_action(restored,"actionLoad_project",{project}); select_station(restored,"GUIDE"); settle();
  check(widget<QAction>(restored,"actionShow_rho_plus")->isChecked(),"Project lost rho+ visibility");
  widget<QAction>(restored,"actionImpedance_validity")->trigger();
  auto *restoredValidity=widget<QDialog>(restored,"impedanceValidityWindow");
  check(!restoredValidity->findChild<QDoubleSpinBox *>("validityErrorFloor"),"Project restored an obsolete floor control");
  check(widget<QSpinBox>(*restoredValidity,"validityGridDensity")->value()==16,"Project lost rho+ grid density");
  auto *list=widget<QListWidget>(window,"stationList"); list->currentItem()->setCheckState(Qt::Unchecked); settle();
  check(!xy->visible(),"Disabled station retained a guide");
  list->currentItem()->setCheckState(Qt::Checked); settle(); check(xy->visible(),"Re-enabled station lost its guide");
  toggle->setChecked(false); check(!xy->visible() && !yxPhase->visible(),"Overlay switch ignored");
  select_station(window,"OUTLIER");
  widget<QAction>(window,"actionImpedance_validity")->trigger();
  auto *preview=widget<QCheckBox>(*validity,"validityPreviewMasks");
  auto *apply=widget<QPushButton>(*validity,"validityApplyMasks");
  auto *undo=widget<QPushButton>(*validity,"validityUndoMasks");
  preview->setChecked(true); settle(); check(apply->isEnabled(),"Main-window masking preview missing");
  widget<QAction>(window,"actionLink_tensor_masks")->setChecked(false);
  check(!apply->isEnabled(),"Link setting left a stale masking preview");
  preview->setChecked(true); settle(); apply->click();
  check(rho->graph(5)->data()->size()>0 && phase->graph(5)->data()->size()>0,"Auto-masking did not refresh main plots");
  check(undo->isEnabled(),"Main-window refresh lost undo");
  undo->click(); check(rho->graph(5)->data()->isEmpty() && phase->graph(5)->data()->isEmpty(),"Undo did not refresh main plots");
  tf_view_controls(window, directory, project);
}

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
  const auto stationItem = [&] { return widget<QListWidget>(window, "stationList")->findItems("S01", Qt::MatchExactly).front(); };
  check(stationItem()->data(StationListDelegate::PeriodCount).toUInt() == 4 &&
        stationItem()->data(StationListDelegate::Coordinates).toString() == "Lat 55.0000°  Lon 12.0000°",
        "Station list omitted period count or geographic coordinates");
  check(stationItem()->text() == "S01" && widget<QLabel>(window, "tfStationName")->text() == "S01", "Station metadata changed identity");
  auto *stationMap = widget<QCustomPlot>(window, "mapPlot");
  const auto firstPosition = stationMap->graph(0)->data()->constBegin()->value;
  auto *mapMode = widget<QComboBox>(window, "mapCoordinateSwitch");
  auto *centerUtm = widget<QCheckBox>(window, "mapCenterUtm");
  auto *originInfo = widget<QLabel>(window, "mapCoordinateInfo");
  check(!centerUtm->isEnabled(), "Geographic map enabled UTM centering");
  mapMode->setCurrentIndex(1);
  check(centerUtm->isEnabled() && !centerUtm->isChecked(), "UTM centering unavailable");
  const auto absolute = SurveyCoordinates::calculate(survey->geographic_locations(), 33, true, false).transform(survey->geographic_locations());
  const auto coordinateText = [](const std::array<double, 3> &p) { return QString("E %1  N %2 m").arg(p[1], 0, 'f', 1).arg(p[0], 0, 'f', 1); };
  check(stationItem()->data(StationListDelegate::Coordinates).toString() == coordinateText(absolute.front()) &&
        stationItem()->toolTip().contains("UTM 33N"), "Station list does not follow UTM map coordinates");
  const auto expectedCoordinates = SurveyCoordinates::calculate(survey->geographic_locations(), 33, true, true);
  centerUtm->setChecked(true);
  check(stationItem()->data(StationListDelegate::Coordinates).toString() ==
          coordinateText(expectedCoordinates.transform(survey->geographic_locations()).front()), "Station list did not subtract UTM origin");
  const QString expectedOrigin = QString::number(expectedCoordinates.origin_easting, 'f', 3);
  check(originInfo->text().contains(expectedOrigin) &&
        originInfo->text().contains(QString::number(expectedCoordinates.origin_northing, 'f', 3)),
        "Centered map omitted actual UTM origin");
  check(std::abs(stationMap->graph(0)->data()->constBegin()->value) < 10000., "Map did not subtract origin");
  widget<QAction>(window,"actionUndo")->trigger();
  check(!centerUtm->isChecked() && originInfo->text().contains("Origin E: 0.000 m"), "Shared undo did not restore map origin");
  widget<QAction>(window,"actionRedo")->trigger();
  check(centerUtm->isChecked() && originInfo->text().contains(expectedOrigin), "Shared redo did not restore map origin");
  centerUtm->setChecked(false);
  check(originInfo->text().contains("Origin E: 0.000 m") &&
        std::abs(stationMap->graph(0)->data()->constBegin()->value) > 100000., "Absolute UTM not restored");
  centerUtm->setChecked(true);
  mapMode->setCurrentIndex(0);
  check(stationItem()->data(StationListDelegate::Coordinates).toString() == "Lat 55.0000°  Lon 12.0000°",
        "Switching to Lat/Lon left projected station coordinates");
  check(!centerUtm->isEnabled() && centerUtm->isChecked() && originInfo->text().contains(expectedOrigin),
        "Geographic view lost shared survey origin");
  widget<QCheckBox>(window, "stationMapNames")->setChecked(true);
  check(stationMap->itemCount() == 2 && widget<QAction>(window, "actionShow_station_names")->isChecked(),
        "Map name checkbox and menu action are out of sync");
  widget<QAction>(window, "actionShow_station_names")->setChecked(false);
  check(stationMap->itemCount() == 0 && !widget<QCheckBox>(window, "stationMapNames")->isChecked(),
        "Names did not disappear from station map");
  checkNear(stationMap->graph(0)->data()->constBegin()->value, firstPosition);
  widget<QCheckBox>(window, "stationMapNames")->setChecked(true);
  file_action(window, "actionLoad_responses", {first, second});
  auto *list = widget<QListWidget>(window, "responsesList");
  check(list->count() == 2 && list->currentRow() == 1, "Iteration files were not loaded/selected");
  widget<QPushButton>(window, "openPeriodMaps")->click();
  auto *maps = widget<QDialog>(window, "periodMapWindow");
  check(maps->isVisible() && widget<QComboBox>(*maps, "mapDataset")->count() == 3, "Maps button did not open loaded datasets");
  check(widget<QCheckBox>(*maps, "mapStationNames")->isChecked(), "Period map did not inherit main map station names");
  widget<QCheckBox>(*maps, "mapStationNames")->setChecked(false);
  maps->hide();
  widget<QAction>(window, "actionPeriod_maps")->trigger();
  check(maps->isVisible(), "Data menu did not reopen period maps");
  check(!widget<QCheckBox>(*maps, "mapStationNames")->isChecked(), "Reopening period map reset its display choices");
  maps->hide();
  auto &rho = *widget<QCustomPlot>(window, "plot11");
  auto &phase = *widget<QCustomPlot>(window, "plot12");
  auto &tipper = *widget<QCustomPlot>(window, "plot21");
  auto &tensor = *widget<QCustomPlot>(window, "plot22");
  const auto beforeReportRange = rho.xAxis->range();
  const auto beforeReportValue = value(rho, 9, 1.);
  export_report(window, directory + "/survey-report.pdf");
  checkNear(rho.xAxis->range().lower, beforeReportRange.lower); checkNear(rho.xAxis->range().upper, beforeReportRange.upper);
  checkNear(value(rho, 9, 1.), beforeReportValue);
  const std::array<SurveyReport::Options::Axis, 4> guiAxes{{
    {false, 1., 1000.}, {false, -10., 100.}, {true, -.4, .4}, {true, -2., 2.}}};
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
  checkNear(value(rho, 9, 1.), 10.); checkNear(value(rho, 9, .25), 40.);
  check(std::isnan(value(rho, 9, .5)), "Missing response should break the line");
  checkNear(value(phase, 9, .5), 60.); checkNear(value(phase, 10, 1.), 45.);
  check(tipper.graph(8)->data()->isEmpty() && tensor.graph(8)->data()->isEmpty(),
        "TF strips retained scalar response curves");
  check(tipper.graph(9)->data()->isEmpty() && tensor.graph(9)->data()->isEmpty(),
        "Absent response component produced a curve");
  widget<QAction>(window, "actionFit_statistics")->trigger();
  auto *analysis = widget<QDialog>(window, "fitStatisticsWindow");
  auto *checks = widget<QTreeWidget>(*analysis, "fitResponses");
  auto *fitPeriods = widget<QCustomPlot>(*analysis, "fitPeriods");
  check(checks->topLevelItemCount() == 2 && fitPeriods->graphCount() == 2, "Missing response comparisons");
  checkNear(checks->topLevelItem(0)->text(1).toDouble(), 0.);
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
    checkNear(count, matched);
  }
  auto *normalization = widget<QComboBox>(*analysis, "fitErrorSource");
  normalization->setCurrentIndex(1);
  check(checks->topLevelItem(0)->text(2) == "4 / 4", "Normalization lost paired scalars");
  normalization->setCurrentIndex(0);
  // A map edit immediately refreshes both the main curves and the open fit window.
  maps->show(); QApplication::processEvents();
  auto *mapPlot = widget<QCustomPlot>(*maps, "periodMapPlot");
  widget<QComboBox>(*maps, "mapPeriod")->setCurrentIndex(3); // 1 s (four sorted periods)
  checkNear(widget<QComboBox>(*maps, "mapPeriod")->currentData().toDouble(), 1.);
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
  check(tipper.graph(0)->visible() && !tipper.graph(2)->visible() && tipper.legend->itemCount() == 2 &&
        tipper.graph(0)->name() == "Real", "Tipper vector visibility is not identified");
  check(tipper.findChildren<QCheckBox *>().size() == 2, "Arrow legend should have one checkbox per vector");
  component_box(tipper, 0)->click();
  check(tipper.itemCount() == 0 && tipper.legend->itemCount() == 2, "Hidden arrows should leave usable legend controls");
  component_box(tipper, 0)->click();
  check(tipper.graph(0)->visible() && tipper.graph(0)->name() == "Real", "Restoring arrows lost the selected vector");
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
  checkNear(value(phase, 10, 1.), -135.);
  list->setCurrentRow(0);
  checkNear(value(rho, 9, 1.), 1.0132118364233778);
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
  checkNear(value(rho, 9, 1.), 99.);

  // The response table is stored in the project, independently of source files.
  widget<QAction>(window, "actionLink_tensor_masks")->setChecked(false);
  widget<QAction>(window, "actionSave_project")->trigger();
  write(first, "invalid source file\n");
  MainWindow restored;
  file_action(restored, "actionLoad_project", {project});
  check(!widget<QAction>(restored, "actionLink_tensor_masks")->isChecked(), "Project lost independent masking option");
  check(widget<QCheckBox>(restored, "mapCenterUtm")->isChecked() &&
        widget<QLabel>(restored, "mapCoordinateInfo")->text().contains(expectedOrigin),
        "Project reload lost centered origin");
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
  checkNear(value(*widget<QCustomPlot>(restored, "plot11"), 9, 1.), 99.);
  widget<QAction>(restored, "actionFit_statistics")->trigger();
  auto *restoredAnalysis = widget<QDialog>(restored, "fitStatisticsWindow");
  auto *restoredChecks = widget<QTreeWidget>(*restoredAnalysis, "fitResponses");
  check(restoredChecks->topLevelItem(0)->text(2) == "2 / 2", "Project lost native scalar response metadata");
  // Existing GoFEM import remains available and selects the loaded response.
  file_action(restored, "actionLoad_responses", {observations});
  check(restoredList->count() == 3, "GoFEM response import regressed");
  checkNear(value(*widget<QCustomPlot>(restored, "plot11"), 9, 1.), 1.0132118364233778);
  auto checkGoFEM = [&] {
    auto *tree = widget<QTreeWidget>(*restoredAnalysis, "fitResponses");
    QTreeWidgetItem *gofem = nullptr;
    for(int i = 0; i < tree->topLevelItemCount(); ++i)
      if(tree->topLevelItem(i)->data(0, Qt::UserRole).toString() == observations) gofem = tree->topLevelItem(i);
    check(gofem && gofem->text(2) == "10 / 10", "GoFEM statistics did not use the original scalar rows");
    checkNear(gofem->data(1, Qt::UserRole).toDouble(), 0.);
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
  check(widget<QAction>(restored, "actionLink_tensor_masks")->isChecked(), "Older projects should default to linked masking");
  select_station(restored, "S01");
  for(const auto *name: {"plot11", "plot12", "plot22"})
    check_visibility(*widget<QCustomPlot>(restored, name), 6);
  check_visibility(*widget<QCustomPlot>(restored, "plot21"), 5);
  const auto version6Project = directory + "/independent-legends-v6.mtd";
  auto *projectedTipper = widget<QCustomPlot>(restored, "plot21");
  check(projectedTipper->graph(0)->name() == "Real Tzx" && projectedTipper->graph(2)->name() == "Imaginary Tzx",
        "Legacy component choices lost their arrow projection");
  file_action(restored, "actionLoad_responses", {second});
  bool projectedResponse = false;
  for(int i = 0; i < projectedTipper->itemCount(); ++i)
    if(auto *arrow = dynamic_cast<QCPItemLine *>(projectedTipper->item(i)))
      if(arrow->pen().style() == Qt::DashLine && std::abs(arrow->end->coords().x() - 18.) < 1e-10) {
        checkNear(arrow->end->coords().y(), 0.); projectedResponse = true;
      }
  check(projectedResponse, "Legacy projected response arrows are missing");
  {
    std::ofstream output(version6Project.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(output);
    const std::uint32_t magic = 0x45444954, version = 6;
    archive << magic << version << survey << std::map<std::string, MTSurveyData>{} << LegacyOptionsV6{};
  }
  widget<QAction>(restored, "actionLink_tensor_masks")->setChecked(false);
  file_action(restored, "actionLoad_project", {version6Project}); select_station(restored, "S01");
  check(widget<QAction>(restored, "actionLink_tensor_masks")->isChecked(), "Version 6 projects should default to linked masks");
  check_visibility(*widget<QCustomPlot>(restored, "plot11"), 6); check_visibility(*widget<QCustomPlot>(restored, "plot12"), 9);
  check_visibility(*widget<QCustomPlot>(restored, "plot21"), 5); check_visibility(*widget<QCustomPlot>(restored, "plot22"), 10);
  const auto version8Project = directory + "/tf-view-v8.mtd";
  {
    std::ofstream output(version8Project.toStdString(), std::ios::binary);
    boost::archive::binary_oarchive archive(output);
    const std::uint32_t magic = 0x45444954, version = 8;
    LegacyOptionsV8 options;
    options.axes[0] = {false, 2., 2000.}; options.axes[1] = {false, -160., 160.};
    archive << magic << version << survey << std::map<std::string, MTSurveyData>{} << options;
  }
  file_action(restored, "actionLoad_project", {version8Project}); select_station(restored, "S01");
  auto *oldRho = widget<QCustomPlot>(restored, "plot11");
  check(oldRho->yAxis->range() == QCPRange(2, 2000) && oldRho->graph(0)->valueAxis()->range() == QCPRange(2, 2000),
        "Version 8 fixed ranges were not migrated to both columns");
  check(widget<QAction>(restored, "actionShow_rho_plus")->isChecked() &&
        widget<QDoubleSpinBox>(restored, "viewMarkerSize")->value() == 5. &&
        widget<QComboBox>(restored, "viewColorMap")->currentText() == "Viridis", "Version 8 plot options were not preserved");
}

void real_run(const QString &project, const QString &directory, const QString &screenshot, bool mixed)
{
  MainWindow window;
  file_action(window, "actionLoad_project", {project});
  for(const auto *name: {"plot11", "plot12", "plot21", "plot22"})
    for(auto *box: widget<QCustomPlot>(window, name)->findChildren<QCheckBox *>())
      check(box->isChecked(),
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
  checkNear(value(rho, 9, 100.), expected);
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
    if(ok && expectedRms.count(iteration)) checkNear(item->data(1, Qt::UserRole).toDouble(), expectedRms.at(iteration));
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
  if(!screenshot.isEmpty()) {
    export_report(window, screenshot + ".survey.pdf", nullptr, [](QDialog &dialog) {
      widget<QComboBox>(dialog, "reportTipperStyle")->setCurrentIndex(1);
      for(int mask: {1, 2, 4, 8}) widget<QAction>(dialog, QString("mapLayer%1").arg(mask).toUtf8().constData())->setChecked(true);
    });
    export_report(window, screenshot + ".selected-maps.pdf", nullptr, [](QDialog &dialog) {
      widget<QCheckBox>(dialog, "reportStationPages")->setChecked(false);
      widget<QCheckBox>(dialog, "reportPTMaps")->setChecked(true);
      widget<QCheckBox>(dialog, "reportInductionMaps")->setChecked(true);
      widget<QCheckBox>(dialog, "reportRealArrows")->setChecked(true);
      widget<QCheckBox>(dialog, "reportImagArrows")->setChecked(true);
      widget<QComboBox>(dialog, "reportMapData")->setCurrentIndex(2);
      widget<QLineEdit>(dialog, "reportMapPeriods")->setText("1, 10");
      for(int mask: {1, 2, 4, 8}) widget<QAction>(dialog, QString("mapLayer%1").arg(mask).toUtf8().constData())->setChecked(true);
    }, 5);
  }
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
    else { tensor_strip(); station_selection_edits(directory.path()); rho_plus_overlay(directory.path()); map_backgrounds(directory.path()); linked_mask_tests(directory.path()); edi_import_tests(directory.path()); statistics(directory.path()); gofem_responses(directory.path()); period_maps(directory.path()); period_resampling_tests(directory.path()); resampling_project_workflow(directory.path()); survey_reports(directory.path()); workflow(directory.path()); }
    std::cout << "Response import, curves and project workflow checks passed.\n";
  } catch(const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
