#include "PolarityWindow.h"
#include "include/MTSurveyData.h"
#include "include/SurveyHistory.h"
#include "qcustomplot.h"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>

PolarityWindow::PolarityWindow(QWidget *parent, SurveyHistory &history, std::function<void()> changed)
  : QDialog(parent, Qt::Window), history(history), changed(std::move(changed)) {
  setObjectName("polarityWindow");
  setWindowTitle(tr("Polarity review"));
  resize(1100, 720);
  auto *layout = new QVBoxLayout(this);
  auto *controls = new QHBoxLayout;
  layout->addLayout(controls);
  controls->addWidget(new QLabel(tr("Neighbours"), this));
  neighbours = new QSpinBox(this); neighbours->setRange(3, 20); neighbours->setValue(5);
  neighbours->setKeyboardTracking(false);
  neighbours->setObjectName("polarityNeighbours"); controls->addWidget(neighbours);
  neighbours->setToolTip(tr("Compare this many nearest enabled stations within the distance limit. Coordinates must use the same x/y convention. A survey-wide reversal cannot be identified from neighbours."));
  controls->addWidget(new QLabel(tr("Max distance"), this));
  radius = new QDoubleSpinBox(this); radius->setRange(.1, 20000.); radius->setValue(100.); radius->setSuffix(tr(" km"));
  radius->setKeyboardTracking(false);
  radius->setObjectName("polarityRadius"); controls->addWidget(radius);
  radius->setToolTip(tr("Maximum WGS84 geodesic distance. Stations without positions can be corrected manually but cannot be compared automatically."));
  controls->addWidget(new QLabel(tr("Tolerance"), this));
  tolerance = new QDoubleSpinBox(this); tolerance->setRange(1., 60.); tolerance->setValue(30.); tolerance->setSuffix(tr("°"));
  tolerance->setKeyboardTracking(false);
  tolerance->setObjectName("polarityTolerance"); controls->addWidget(tolerance);
  tolerance->setToolTip(tr("Maximum circular phase difference from a neighbour consensus, after a possible 180° reversal. At least 3 neighbours and 75% agreement are required per period. A candidate needs at least 3 matched periods, 60% coverage of usable periods and 80% reversed matches. SNR must be at least 3; periods must match within 0.1%, without interpolation."));
  controls->addStretch();
  auto *notice = new QLabel(tr("Candidates need review: real structure and coordinate conventions can also reverse phases."), this);
  notice->setWordWrap(true); layout->addWidget(notice);
  auto *splitter = new QSplitter(this); layout->addWidget(splitter, 1);
  stations = new QTableWidget(splitter);
  stations->setObjectName("polarityStations");
  stations->setColumnCount(3); stations->setHorizontalHeaderLabels({tr("Station"), "XY", "YX"});
  stations->setSelectionBehavior(QAbstractItemView::SelectRows);
  stations->setSelectionMode(QAbstractItemView::SingleSelection);
  stations->setEditTriggers(QAbstractItemView::NoEditTriggers);
  stations->verticalHeader()->hide(); stations->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  auto *right = new QWidget(splitter); auto *plotLayout = new QVBoxLayout(right);
  plotLayout->setContentsMargins(0, 0, 0, 0);
  for(unsigned c = 0; c < 2; ++c) {
    auto *plot = plots[c] = new QCustomPlot(right);
    plot->setObjectName(c == 0 ? "polarityXY" : "polarityYX");
    plotLayout->addWidget(plot, 1);
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    plot->xAxis->setLabel(tr("Period (s)"));
    plot->yAxis->setLabel(c == 0 ? tr("XY phase (°)") : tr("YX phase (°)"));
    auto ticks = QSharedPointer<QCPAxisTickerFixed>(new QCPAxisTickerFixed);
    ticks->setTickStep(90.); plot->yAxis->setTicker(ticks);
    plot->yAxis->setRange(-190., 190.);
    plot->legend->setVisible(true);
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot->setToolTip(tr("Raw phases in [-180°, 180°], without the main plot's 180° phase wrapping. Grey crosses show the local consensus. Preview changes only the selected field-channel signs."));
  }
  splitter->setSizes({360, 740}); splitter->setStretchFactor(1, 1);
  auto *corrections = new QHBoxLayout; layout->addLayout(corrections);
  corrections->addWidget(new QLabel(tr("Reverse"), this));
  const QStringList names{"Ex", "Ey", "Hx", "Hy", "Hz"};
  const QStringList help{
    tr("Reverse Zxx and Zxy at all periods. Tippers and phase tensor stay unchanged."),
    tr("Reverse Zyx and Zyy at all periods. Tippers and phase tensor stay unchanged."),
    tr("Reverse Zxx, Zyx and Tzx at all periods, and the off-diagonal phase-tensor entries."),
    tr("Reverse Zxy, Zyy and Tzy at all periods, and the off-diagonal phase-tensor entries."),
    tr("Reverse both tippers at all periods. Impedance and phase tensor stay unchanged.")};
  for(unsigned c = 0; c < channels.size(); ++c) {
    channels[c] = new QCheckBox(names[c], this); channels[c]->setObjectName("polarity" + names[c]);
    channels[c]->setToolTip(help[c]); corrections->addWidget(channels[c]);
    connect(channels[c], &QCheckBox::toggled, this, [this] { preview(); });
  }
  suggestion = new QPushButton(tr("Suggest Ex / Ey"), this);
  suggestion->setObjectName("polaritySuggest"); corrections->addWidget(suggestion);
  suggestion->setToolTip(tr("Preview an electric-channel correction for flagged rows. Off-diagonal phases alone cannot distinguish an electric reversal from a magnetic reversal. Confirm the wiring/convention before applying; Hz is not detected automatically."));
  corrections->addStretch();
  applyButton = new QPushButton(tr("Apply to station"), this);
  applyButton->setObjectName("polarityApply"); corrections->addWidget(applyButton);
  applyButton->setToolTip(tr("Apply these channel reversals to every stored period, including masked values. Errors, error floors and masks are preserved. Original input files and calculated responses are untouched. Use the main Undo command to restore the station."));
  status = new QLabel(this); status->setObjectName("polarityStatus"); status->setWordWrap(true); layout->addWidget(status);
  connect(stations, &QTableWidget::currentCellChanged, this, [this] { select(); });
  connect(neighbours, qOverload<int>(&QSpinBox::valueChanged), this, [this] { analyse(); });
  for(auto *box : {radius, tolerance}) connect(box, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this] { analyse(); });
  connect(suggestion, &QPushButton::clicked, this, [this] {
    const int row = stations->currentRow(); if(row < 0) return;
    for(unsigned c = 0; c < channels.size(); ++c) {
      const QSignalBlocker block(channels[c]);
      channels[c]->setChecked(c < 2 && results[row].components[c].candidate());
    }
    preview();
  });
  connect(applyButton, &QPushButton::clicked, this, [this] { apply(); });
  select();
}
void PolarityWindow::setData(std::shared_ptr<MTSurveyData> data, const QString &selected) {
  survey = std::move(data); analyse();
  for(unsigned i = 0; i < results.size(); ++i)
    if(QString::fromStdString(results[i].name) == selected) { stations->setCurrentCell(i, 0); break; }
}
void PolarityWindow::analyse() {
  const auto *item = stations->item(stations->currentRow(), 0);
  const QString previous = item ? item->text() : QString();
  results = survey ? PolarityAnalysis::analyse(*survey, {unsigned(neighbours->value()), radius->value(), tolerance->value()})
                   : std::vector<PolarityAnalysis::Station>{};
  const auto candidate = [](const PolarityAnalysis::Station &s) { return s.components[0].candidate() || s.components[1].candidate(); };
  std::stable_sort(results.begin(), results.end(), [&](const auto &a, const auto &b) { return candidate(a) > candidate(b); });
  {
    const QSignalBlocker block(stations);
    stations->setRowCount(int(results.size()));
    int selected = 0;
    for(unsigned i = 0; i < results.size(); ++i) {
      const auto name = QString::fromStdString(results[i].name);
      stations->setItem(i, 0, new QTableWidgetItem(name));
      if(name == previous) selected = i;
      for(unsigned c = 0; c < 2; ++c) {
        const auto &result = results[i].components[c];
        const auto label = result.candidate() ? tr("Candidate") : result.supported() ? tr("No flag") : tr("Insufficient");
        auto *cell = new QTableWidgetItem(tr("%1 (%2/%3)").arg(label).arg(result.reversed).arg(result.points.size()));
        cell->setToolTip(tr("%1 reversed / %2 matched periods; %3 usable target periods. A flag is evidence relative to neighbours, not a confirmed wiring error.")
          .arg(result.reversed).arg(result.points.size()).arg(result.usable));
        if(result.candidate()) { auto font = cell->font(); font.setBold(true); cell->setFont(font); }
        stations->setItem(i, c + 1, cell);
      }
    }
    if(!results.empty()) stations->setCurrentCell(selected, 0);
  }
  select();
}
std::array<bool, 5> PolarityWindow::selection() const {
  std::array<bool, 5> result;
  for(unsigned c = 0; c < channels.size(); ++c) result[c] = channels[c]->isChecked();
  return result;
}
void PolarityWindow::select() {
  for(auto *channel : channels) { const QSignalBlocker block(channel); channel->setChecked(false); }
  preview();
}
void PolarityWindow::preview() {
  const int row = stations->currentRow();
  const bool valid = survey && row >= 0 && row < int(results.size());
  const auto signs = selection();
  const bool changes = std::any_of(signs.begin(), signs.end(), [](bool on) { return on; }) &&
                       !std::all_of(signs.begin(), signs.end(), [](bool on) { return on; });
  applyButton->setEnabled(valid && changes);
  suggestion->setEnabled(valid && (results[row].components[0].candidate() || results[row].components[1].candidate()));
  for(auto *channel : channels) channel->setEnabled(valid);
  for(auto *plot : plots) plot->clearGraphs();
  status->setText(valid ? tr("%1 — preview only. Apply reverses the selected channels at every period; Undo restores the edit.").arg(QString::fromStdString(results[row].name)) : tr("Load a survey to review polarity."));
  if(valid) {
    const auto &station = survey->get_station_data(results[row].name);
    auto corrected = station; corrected.reverse_polarity(signs);
    for(unsigned c = 0; c < 2; ++c) {
      auto *plot = plots[c];
      const QColor color = c == 0 ? QColor(0, 90, 180) : QColor(200, 85, 0);
      const auto graph = [&](const QString &name, QColor ink, QCPScatterStyle::ScatterShape shape) {
        auto *line = plot->addGraph(); line->setName(name); line->setPen(QPen(ink));
        line->setLineStyle(QCPGraph::lsNone); line->setScatterStyle(QCPScatterStyle(shape, 6)); return line;
      };
      auto *observed = graph(tr("Observed"), color, QCPScatterStyle::ssCircle);
      auto *reference = graph(tr("Neighbours"), Qt::gray, QCPScatterStyle::ssCross);
      auto *preview = changes ? graph(tr("Preview"), color.darker(), QCPScatterStyle::ssPlus) : nullptr;
      for(const auto &point : results[row].components[c].points) reference->addData(point.period, point.referencePhase);
      for(unsigned i = 0; i < station.frequencies().size(); ++i) {
        const double f = station.frequencies()[i]; if(!std::isfinite(f) || f <= 0.) continue;
        double value, error;
        const auto type = c == 0 ? PhsZxy : PhsZyx;
        if(station.scalar_value(type, i, value, error)) observed->addData(1. / f, std::remainder(value, 360.));
        if(preview && corrected.scalar_value(type, i, value, error)) preview->addData(1. / f, std::remainder(value, 360.));
      }
      plot->rescaleAxes(); plot->yAxis->setRange(-190., 190.);
    }
  }
  for(auto *plot : plots) plot->replot();
}
void PolarityWindow::apply() {
  const int row = stations->currentRow();
  if(!survey || row < 0 || !applyButton->isEnabled()) return;
  const auto name = results[row].name;
  QStringList reversed;
  for(auto *channel : channels) if(channel->isChecked()) reversed << channel->text();
  survey->get_station_data(name).reverse_polarity(selection());
  const bool recorded = history.record(tr("Reverse %1: %2").arg(QString::fromStdString(name), reversed.join(", ")));
  analyse();
  if(recorded && changed) changed();
  if(!recorded) status->setText(tr("No stored values changed."));
}
