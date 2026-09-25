#include "PeriodLayoutWindow.h"
#include "include/HelpButton.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>

using namespace PeriodResampling;
namespace {
double number(const QLineEdit *edit) {
  bool ok = false;
  const double value = edit->text().trimmed().toDouble(&ok);
  if(!ok) throw std::invalid_argument("Enter periods as numbers in seconds (scientific notation is supported).");
  return value;
}
QColor color(Status status) {
  switch(status) {
    case Status::Existing: return QColor("#2563eb");
    case Status::Interpolated: return QColor("#16a34a");
    case Status::Derived: return QColor("#0d9488");
    case Status::Masked: return QColor("#b45309");
    default: return QColor("#999999");
  }
}
QString csv(QString value) { value.replace('"', "\"\""); return '"' + value + '"'; }
QString statusLabel(Status status) {
  switch(status) {
    case Status::Existing: return QObject::tr("Original");
    case Status::Interpolated: return QObject::tr("Estimated");
    case Status::Derived: return QObject::tr("Derived");
    case Status::Masked: return QObject::tr("Masked");
    case Status::Outside: return QObject::tr("Out of range");
    case Status::Gap: return QObject::tr("Large gap");
    case Status::InvalidError: return QObject::tr("Invalid error");
    default: return QObject::tr("Missing");
  }
}
}

PeriodLayoutWindow::PeriodLayoutWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> create)
  : QDialog(parent, Qt::Window), createSurvey(std::move(create))
{
  setObjectName("periodLayoutWindow"); setWindowTitle(tr("Period layout and resampling")); resize(1200, 900);
  auto *layout = new QVBoxLayout(this);
  auto *controls = new QGridLayout;
  gridMode = new QComboBox(this); gridMode->setObjectName("layoutGridMode");
  gridMode->addItems({tr("Logarithmic grid"), tr("Reference station"), tr("Custom periods")});
  controls->addWidget(UiHelp::label(this, tr("Grid:"), tr("Choose a logarithmic grid, the periods of a reference station, or a custom list in seconds. Custom periods accept spaces, commas or semicolons and are sorted and deduplicated. A logarithmic grid includes both limits; its last interval may be shorter."), "layoutGridHelp"), 0, 0); controls->addWidget(gridMode, 0, 1);
  minimum = new QLineEdit("0.01", this); minimum->setObjectName("layoutMinimum");
  maximum = new QLineEdit("1000", this); maximum->setObjectName("layoutMaximum");
  density = new QSpinBox(this); density->setObjectName("layoutDensity"); density->setRange(1, 100); density->setValue(6);
  controls->addWidget(new QLabel(tr("Min (s):"), this), 0, 2); controls->addWidget(minimum, 0, 3);
  controls->addWidget(new QLabel(tr("Max (s):"), this), 0, 4); controls->addWidget(maximum, 0, 5);
  controls->addWidget(UiHelp::label(this, tr("Pts/decade:"), tr("Number of grid intervals per tenfold increase in period. Higher density creates more estimates, not more independent measurements."), "layoutDensityHelp"), 0, 6); controls->addWidget(density, 0, 7);
  reference = new QComboBox(this); reference->setObjectName("layoutReference");
  controls->addWidget(new QLabel(tr("Reference:"), this), 1, 0); controls->addWidget(reference, 1, 1);
  custom = new QLineEdit(this); custom->setObjectName("layoutCustomPeriods");
  custom->setPlaceholderText(tr("e.g. 0.1, 1, 10"));
  controls->addWidget(new QLabel(tr("Periods (s):"), this), 1, 2); controls->addWidget(custom, 1, 3, 1, 5);
  gap = new QDoubleSpinBox(this); gap->setObjectName("layoutMaximumRatio");
  gap->setDecimals(3); gap->setRange(1.001, 10000.); gap->setValue(2.); gap->setKeyboardTracking(false);
  gap->setToolTip(tr("Maximum ratio of upper to lower bounding periods. A limit of 2 allows 1–2 s, but blocks 1–10 s, including targets near either endpoint."));
  controls->addWidget(UiHelp::label(this, tr("Max gap ratio:"), gap->toolTip(), "layoutGapHelp"), 2, 0, 1, 3); controls->addWidget(gap, 2, 3);
  barriers = new QCheckBox(tr("Respect masks"), this);
  barriers->setObjectName("layoutMaskBarriers"); barriers->setChecked(true);
  auto *maskRow = new QHBoxLayout; maskRow->addWidget(barriers);
  maskRow->addWidget(UiHelp::button(this, tr("Mask barriers"), tr("Masked samples block interpolation across them. Turn this off to bridge masks within the gap limit. Exact masked samples remain masked in either mode."), "layoutMasksHelp")); maskRow->addStretch();
  controls->addLayout(maskRow, 2, 4, 1, 4);
  layout->addLayout(controls);
  auto *filters = new QHBoxLayout;
  filters->addWidget(new QLabel(tr("Component:"), this));
  group = new QComboBox(this); group->setObjectName("layoutGroup"); group->addItems({tr("Impedance"), tr("Tipper"), tr("Phase tensor")}); filters->addWidget(group);
  component = new QComboBox(this); component->setObjectName("layoutComponent"); filters->addWidget(component);
  component->setToolTip(tr("Impedance and tipper availability requires both real and imaginary parts. Matching partial source samples are retained, but do not form interpolation endpoints."));
  filters->addWidget(UiHelp::button(this, tr("Preview component"), tr("These filters change the preview only. All components are resampled. Impedance and tipper endpoints need both real and imaginary parts; exact partial samples are retained."), "layoutComponentHelp")); filters->addStretch();
  filters->addWidget(UiHelp::label(this, tr("Method"), tr("Real and imaginary values are interpolated linearly against log(period), with no extrapolation. Derived quantities are recomputed from impedance; tensor-only stations use direct tensor interpolation. Standard errors use weighted endpoint errors. Interpolated samples share source information and are not independent measurements."), "layoutMethodHelp"));
  filters->addWidget(UiHelp::label(this, tr("Coverage"), tr("Coverage counts enabled stations with a complete component at each exact period. Source shows the original periods; Original on grid counts original samples at target periods; Usable on grid includes interpolated estimates. An increase means more estimates on the chosen grid, not new measurements or a wider valid period range.\nTarget markers: Original = retained sample; Estimated = interpolation; Derived = calculated from impedance. Hover over a target for its source periods or the reason it was skipped. Counts in the preview summary span all complex and tensor components."), "layoutCoverageHelp"));
  layout->addLayout(filters);
  auto *tabs = new QTabWidget(this); tabs->setObjectName("layoutTabs");
  auto addPlot = [&](const char *name) {
    auto *plot = new QCustomPlot(this); plot->setObjectName(name);
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    plot->xAxis->setLabel(tr("Period (s)")); plot->legend->setVisible(true);
    plot->legend->setBrush(QColor(255, 255, 255, 230));
    plot->axisRect()->insetLayout()->take(plot->legend);
    plot->plotLayout()->insertRow(0); plot->plotLayout()->addElement(0, 0, plot->legend);
    plot->legend->setFillOrder(QCPLayoutGrid::foColumnsFirst); plot->legend->setWrap(4);
    plot->legend->setBorderPen(Qt::NoPen); plot->plotLayout()->setRowStretchFactor(0, .001);
    return plot;
  };
  sourcePlot = addPlot("layoutSourcePlot"); targetPlot = addPlot("layoutTargetPlot"); coverage = addPlot("layoutCoveragePlot");
  tabs->addTab(sourcePlot, tr("Source layout")); tabs->addTab(targetPlot, tr("Target preview")); tabs->addTab(coverage, tr("Coverage"));
  counts = new QTableWidget(this); counts->setObjectName("layoutCounts"); counts->setColumnCount(7);
  counts->setHorizontalHeaderLabels({tr("Station"), tr("Source"), tr("Gaps"), tr("Original"), tr("Estimated"), tr("Missing"), tr("Enabled")});
  const QStringList headerHelp{tr("Station name."), tr("Active source samples for the selected component."), tr("Adjacent valid source periods exceeding the maximum gap ratio."), tr("Existing samples matching a target period."), tr("Interpolated values or values derived from impedance."), tr("Unavailable, masked or incomplete target samples."), tr("Whether this station is enabled; only enabled stations contribute to coverage.")};
  for(int i = 0; i < headerHelp.size(); ++i) counts->horizontalHeaderItem(i)->setToolTip(headerHelp[i]);
  counts->setEditTriggers(QAbstractItemView::NoEditTriggers); counts->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  tabs->addTab(counts, tr("Station summary")); layout->addWidget(tabs, 1);
  sourceSummary = new QLabel(this); sourceSummary->setObjectName("layoutSourceSummary"); sourceSummary->setWordWrap(true); layout->addWidget(sourceSummary);
  summary = new QLabel(this); summary->setObjectName("layoutSummary"); summary->setWordWrap(true); layout->addWidget(summary);
  history = new QLabel(this); history->setObjectName("layoutHistory"); history->setWordWrap(true);
  auto *actions = new QHBoxLayout;
  auto *previewButton = new QPushButton(tr("Preview"), this); previewButton->setObjectName("layoutPreview"); actions->addWidget(previewButton);
  report = new QPushButton(tr("Preview CSV…"), this); report->setObjectName("layoutExportReport"); actions->addWidget(report);
  sourceReport = new QPushButton(tr("Saved audit…"), this); sourceReport->setObjectName("layoutExportSourceReport"); actions->addWidget(sourceReport);
  openSource = new QPushButton(tr("Open source"), this); openSource->setObjectName("layoutOpenSource"); actions->addWidget(openSource);
  actions->addStretch();
  actions->addWidget(history);
  apply = new QPushButton(tr("Create survey"), this); apply->setObjectName("layoutApply"); actions->addWidget(apply);
  actions->addWidget(UiHelp::button(this, tr("Survey copies and audit"), tr("Create survey opens an independent resampled survey in a new window. The original remains open and its file is unchanged. Saving the new project embeds a source snapshot and audit; no separate backup file is created automatically. Open source reopens that snapshot. Preview CSV exports the current calculation; Saved audit exports the saved resampling operation, before later edits."), "layoutSourceHelp"));
  auto *close = new QPushButton(tr("Close"), this); actions->addWidget(close); layout->addLayout(actions);
  apply->setToolTip(tr("Open an independent survey in a new window. Save it as a new project; the original remains open and a source snapshot is retained in the new project."));
  connect(close, &QPushButton::clicked, this, &QDialog::hide);
  connect(previewButton, &QPushButton::clicked, this, &PeriodLayoutWindow::preview);
  connect(apply, &QPushButton::clicked, this, [this] { if(result && createSurvey) createSurvey(std::make_shared<MTSurveyData>(*result)); });
  connect(report, &QPushButton::clicked, this, [this] { if(result) exportReport(result->resampling_info()); });
  connect(sourceReport, &QPushButton::clicked, this, [this] { if(survey) exportReport(survey->resampling_info()); });
  connect(openSource, &QPushButton::clicked, this, [this] {
    if(survey && survey->resampling_source() && createSurvey) createSurvey(std::make_shared<MTSurveyData>(*survey->resampling_source()));
  });
  connect(gridMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { invalidate(); });
  connect(reference, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { invalidate(); });
  for(auto *edit: {minimum, maximum, custom}) connect(edit, &QLineEdit::textChanged, this, [this, edit] { edit->setToolTip(edit->text()); invalidate(); });
  connect(density, QOverload<int>::of(&QSpinBox::valueChanged), this, [this] { invalidate(); });
  connect(gap, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this] { invalidate(); });
  connect(barriers, &QCheckBox::toggled, this, [this] { invalidate(); });
  connect(group, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] {
    const QSignalBlocker block(component); component->clear();
    const unsigned first = group->currentIndex() == 0 ? 0 : (group->currentIndex() == 1 ? 4 : 6);
    const unsigned last = first == 0 ? 4 : (first == 4 ? 6 : 10);
    for(unsigned c = first; c < last; ++c) component->addItem(component_name(c), c);
    if(first == 0) component->setCurrentIndex(1);
    updatePlots();
  });
  connect(component, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] { updatePlots(); });
  for(unsigned c = 0; c < 4; ++c) component->addItem(component_name(c), c);
  component->setCurrentIndex(1);
  for(auto *plot: {sourcePlot, targetPlot}) connect(plot, &QCustomPlot::mouseMove, this, [this, plot](QMouseEvent *event) {
    if(!survey) return;
    const auto names = survey->get_stations_names();
    const int row = std::lround(plot->yAxis->pixelToCoord(event->pos().y()));
    if(row < 0 || row >= static_cast<int>(names.size())) { plot->setToolTip({}); return; }
    const double period = plot->xAxis->pixelToCoord(event->pos().x());
    QString text = QString::fromStdString(names[row]);
    double distance = 12.;
    const unsigned c = component->currentData().toUInt();
    if(plot == targetPlot && result) {
      for(const auto &r: result->resampling_info().records) if(r.station == names[row] && r.component == c) {
        const double dx = std::abs(plot->xAxis->coordToPixel(r.period) - event->pos().x());
        if(dx < distance) {
          distance = dx; text = tr("%1 · %2 s · %3").arg(QString::fromStdString(r.station)).arg(r.period, 0, 'g', 10).arg(status_name(r.status));
          if(r.lower > 0.) text += tr("\nSource periods: %1 – %2 s").arg(r.lower, 0, 'g', 10).arg(r.upper, 0, 'g', 10);
        }
      }
    } else text += tr(" · %1 s").arg(period, 0, 'g', 6);
    plot->setToolTip(text);
  });
  invalidate();
}

void PeriodLayoutWindow::setData(const std::shared_ptr<MTSurveyData> &data)
{
  const bool changed = survey != data;
  survey = data;
  const auto previous = reference->currentText();
  {
    const QSignalBlocker block(reference); reference->clear();
    if(survey) for(const auto &name: survey->get_stations_names()) reference->addItem(QString::fromStdString(name));
    reference->setCurrentIndex(std::max(0, reference->findText(previous)));
  }
  if(changed && survey) {
    std::vector<double> periods;
    for(double p: survey->get_unique_periods()) if(std::isfinite(p) && p > 0.) periods.push_back(p);
    if(!periods.empty()) {
      const QSignalBlocker lo(minimum), hi(maximum);
      minimum->setText(QString::number(periods.front(), 'g', 17)); maximum->setText(QString::number(periods.back(), 'g', 17));
      for(auto *edit: {minimum, maximum}) { edit->setCursorPosition(0); edit->setToolTip(edit->text()); }
    }
  }
  const bool resampled = survey && survey->resampling_source();
  openSource->setVisible(resampled); sourceReport->setVisible(resampled);
  history->setVisible(resampled);
  history->setText(resampled ? tr("Source retained") : QString());
  history->setToolTip(resampled ? tr("Saved grid: %1 periods; maximum gap ratio %2.")
                     .arg(survey->resampling_info().periods.size()).arg(survey->resampling_info().maximumRatio) : QString());
  invalidate();
}

std::vector<double> PeriodLayoutWindow::targetPeriods() const
{
  if(!survey) throw std::invalid_argument("Load a survey first.");
  if(gridMode->currentIndex() == 0) return logarithmic_grid(number(minimum), number(maximum), density->value());
  std::vector<double> periods;
  if(gridMode->currentIndex() == 1) {
    if(!survey->is_station_present(reference->currentText().toStdString())) throw std::invalid_argument("Select a reference station.");
    for(double f: survey->get_station_data(reference->currentText().toStdString()).frequencies()) periods.push_back(1. / f);
  } else {
    for(const auto &word: custom->text().split(QRegularExpression("[\\s,;]+"), Qt::SkipEmptyParts)) {
      bool ok; const double p = word.toDouble(&ok);
      if(!ok) throw std::invalid_argument("Custom periods must be numbers separated by spaces, commas or semicolons.");
      periods.push_back(p);
    }
  }
  return validate_periods(periods);
}

void PeriodLayoutWindow::invalidate()
{
  const int mode = gridMode->currentIndex();
  minimum->setEnabled(mode == 0); maximum->setEnabled(mode == 0); density->setEnabled(mode == 0);
  reference->setEnabled(mode == 1); custom->setEnabled(mode == 2);
  result.reset(); apply->setEnabled(false); report->setEnabled(false);
  summary->setText(tr("Click Preview to check the grid."));
  updatePlots();
}

void PeriodLayoutWindow::preview()
{
  result.reset(); apply->setEnabled(false); report->setEnabled(false);
  try {
    const auto periods = targetPeriods();
    result = resample(*survey, periods, {gap->value(), barriers->isChecked()});
    unsigned existing = 0, interpolated = 0, unavailable = 0;
    for(const auto &r: result->resampling_info().records) if(survey->is_active(r.station)) {
      if(r.status == Status::Existing) ++existing;
      else if(available(r.status)) ++interpolated;
      else ++unavailable;
    }
    summary->setText(tr("Preview: %4 periods · %1 original · %2 estimated / derived · %3 missing")
                     .arg(existing).arg(interpolated).arg(unavailable).arg(periods.size()));
    bool hasValues = existing + interpolated > 0;
    if(!hasValues) for(const auto &name: result->get_stations_names()) {
      const auto &station = result->get_station_data(name);
      for(unsigned f = 0; f < station.frequencies().size() && !hasValues; ++f)
        for(const auto &type: type_to_column_table) {
          double value, error;
          if(station.scalar_value(type.first, f, value, error)) { hasValues = true; break; }
        }
    }
    apply->setEnabled(hasValues); report->setEnabled(true);
  } catch(const std::exception &error) { summary->setText(QString::fromUtf8(error.what())); }
  updatePlots();
}

void PeriodLayoutWindow::drawLayout(QCustomPlot *plot, bool target)
{
  plot->clearGraphs();
  if(!survey) { plot->replot(); return; }
  const auto names = survey->get_stations_names();
  const unsigned c = component->currentData().toUInt();
  auto ticker = QSharedPointer<QCPAxisTickerText>(new QCPAxisTickerText);
  const unsigned tickStep = std::max(1u, unsigned((names.size() + 29) / 30));
  for(unsigned row = 0; row < names.size(); row += tickStep) ticker->addTick(row, QString::fromStdString(names[row]));
  plot->yAxis->setTicker(ticker); plot->yAxis->setRange(-.7, std::max(.7, double(names.size()) - .3));
  plot->yAxis->setRangeReversed(true); plot->yAxis->setLabel(tr("Station"));
  std::map<Status, std::pair<QVector<double>, QVector<double>>> series;
  double lo = std::numeric_limits<double>::infinity(), hi = 0.;
  auto add = [&](double p, unsigned row, Status status) {
    lo = std::min(lo, p); hi = std::max(hi, p);
    series[status].first.push_back(p); series[status].second.push_back(row);
  };
  if(target && result) {
    std::map<std::string, unsigned> rows;
    for(unsigned i = 0; i < names.size(); ++i) rows[names[i]] = i;
    for(const auto &r: result->resampling_info().records) if(r.component == c)
      add(r.period, rows.at(r.station), survey->is_active(r.station) ? r.status : Status::Masked);
  } else if(!target) {
    for(unsigned row = 0; row < names.size(); ++row)
      for(const auto &sample: PeriodResampling::layout(survey->get_station_data(names[row]), c)) add(sample.period, row, sample.status);
  }
  for(const auto &s: series) {
    auto *graph = plot->addGraph(); graph->setName(statusLabel(s.first)); graph->setData(s.second.first, s.second.second);
    graph->setLineStyle(QCPGraph::lsNone);
    graph->setScatterStyle(QCPScatterStyle(available(s.first) ? QCPScatterStyle::ssDisc : QCPScatterStyle::ssCross, color(s.first), 6));
    graph->setSelectable(QCP::stNone);
  }
  if(hi > 0.) plot->xAxis->setRange(lo / 1.15, hi * 1.15);
  else plot->xAxis->setRange(.1, 10.);
  plot->replot();
}

void PeriodLayoutWindow::updatePlots()
{
  if(!component->currentData().isValid()) return;
  drawLayout(sourcePlot, false); drawLayout(targetPlot, true);
  counts->setRowCount(0); coverage->clearGraphs();
  if(!survey) { coverage->replot(); return; }
  const auto names = survey->get_stations_names(); const unsigned c = component->currentData().toUInt();
  std::map<double, unsigned> sourceCoverage, existingCoverage, addedCoverage;
  std::map<std::string, std::array<unsigned, 3>> targets;
  if(result) for(const auto &r: result->resampling_info().records) if(r.component == c) {
    auto &row = targets[r.station];
    ++row[r.status == Status::Existing ? 0 : (available(r.status) ? 1 : 2)];
    existingCoverage.emplace(r.period, 0); addedCoverage.emplace(r.period, 0);
    if(survey->is_active(r.station)) {
      if(r.status == Status::Existing) ++existingCoverage[r.period];
      else if(available(r.status)) ++addedCoverage[r.period];
    }
  }
  unsigned total = 0, largeGaps = 0;
  counts->setRowCount(names.size());
  for(unsigned row = 0; row < names.size(); ++row) {
    unsigned active = 0, gaps = 0; double previous = 0.;
    for(const auto &sample: PeriodResampling::layout(survey->get_station_data(names[row]), c)) {
      sourceCoverage.emplace(sample.period, 0);
      if(sample.status != Status::Existing) continue;
      ++active; ++sourceCoverage[sample.period];
      if(previous > 0. && std::log(sample.period) - std::log(previous) > std::log(gap->value()) + 1e-14) ++gaps;
      previous = sample.period;
    }
    total += active; largeGaps += gaps;
    const auto count = targets[names[row]];
    const QStringList cells{QString::fromStdString(names[row]), QString::number(active), QString::number(gaps),
      result ? QString::number(count[0]) : "—", result ? QString::number(count[1]) : "—", result ? QString::number(count[2]) : "—",
      survey->is_active(names[row]) ? tr("Yes") : tr("No")};
    for(int col = 0; col < cells.size(); ++col) counts->setItem(row, col, new QTableWidgetItem(cells[col]));
  }
  sourceSummary->setText(tr("%1: %2 source samples · %3 stations · %4 large gaps (ratio > %5)")
                         .arg(component->currentText()).arg(total).arg(names.size()).arg(largeGaps).arg(gap->value()));
  auto graph = [&](const QString &name, const std::map<double, unsigned> &data, const QColor &color) {
    QVector<double> x, y;
    for(const auto &p: data) { x.push_back(p.first); y.push_back(p.second); }
    auto *g = coverage->addGraph(); g->setName(name); g->setData(x, y); g->setPen(QPen(color, 1.5));
    g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, color, 4));
  };
  graph(tr("Source"), sourceCoverage, color(Status::Existing));
  if(result) {
    graph(tr("Original on grid"), existingCoverage, QColor("#7c3aed"));
    coverage->graph(1)->setLineStyle(QCPGraph::lsNone);
    for(auto &p: addedCoverage) p.second += existingCoverage[p.first];
    graph(tr("Usable on grid"), addedCoverage, color(Status::Interpolated));
  }
  coverage->yAxis->setLabel(tr("Stations with this component")); coverage->rescaleAxes();
  coverage->yAxis->setRange(0., std::max(1., double(names.size()))); coverage->replot();
}

void PeriodLayoutWindow::exportReport(const Info &info)
{
  auto path = QFileDialog::getSaveFileName(this, tr("Export resampling audit"), "resampling.csv", tr("CSV (*.csv)"));
  if(path.isEmpty()) return;
  if(!path.endsWith(".csv", Qt::CaseInsensitive)) path += ".csv";
  QSaveFile file(path);
  if(!file.open(QIODevice::WriteOnly | QIODevice::Text)) { QMessageBox::warning(this, tr("Export audit"), file.errorString()); return; }
  QTextStream out(&file); out.setLocale(QLocale::c()); out.setRealNumberPrecision(17);
  out << "# " << QString::fromStdString(info.method) << '\n';
  out << "# maximum_bounding_period_ratio=" << info.maximumRatio << "; masked_barriers=" << info.maskedBarriers << '\n';
  out << "station,component,target_period_s,status,source_lower_period_s,source_upper_period_s\n";
  for(const auto &r: info.records)
    out << csv(QString::fromStdString(r.station)) << ',' << component_name(r.component) << ',' << r.period << ','
        << csv(status_name(r.status)) << ',' << r.lower << ',' << r.upper << '\n';
  out.flush();
  if(out.status() != QTextStream::Ok || !file.commit()) QMessageBox::warning(this, tr("Export audit"), tr("Could not write the complete audit file."));
}
