#include "include/SurveyReport.h"
#include "include/ApparentResistivityPlot.h"
#include "include/PhasePlot.h"
#include "include/TipperPlot.h"
#include "include/PhaseTensorPlot.h"
#include "include/FitStatistics.h"
#include "include/HelpButton.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QPdfWriter>
#include <QProgressDialog>
#include <QPushButton>
#include <QSaveFile>
#include <algorithm>
#include <cmath>

namespace {
using namespace SurveyReport;
QString tr(const char *text) { return QObject::tr(text); }
QFont font(int pixels, bool bold = false) { QFont f("Sans Serif"); f.setPixelSize(pixels); f.setBold(bold); return f; }
QString number(double value, int precision = 6) { return std::isfinite(value) ? QString::number(value, 'g', precision) : tr("N/A"); }
void text(QCPPainter &p, const QRectF &rect, const QString &value, int pixels = 16, bool bold = false) {
  p.setPen(QColor("#172033")); p.setFont(font(pixels, bold));
  p.drawText(rect, Qt::AlignLeft | Qt::AlignTop | Qt::TextWordWrap, value);
}
void line(QCPPainter &p, double y, double width) { p.setPen(QPen(QColor("#cbd5e1"), 1)); p.drawLine(QPointF(0, y), QPointF(width, y)); }
void plotAt(QCPPainter &p, QCustomPlot &plot, const QRect &rect) {
  p.save(); p.translate(rect.topLeft()); plot.toPainter(&p, rect.width(), rect.height()); p.restore();
}
void style(QCustomPlot &plot) {
  plot.setBackground(Qt::white);
  for(auto *axis: {plot.xAxis, plot.yAxis}) {
    axis->setTickLabelFont(font(14)); axis->setLabelFont(font(15)); axis->setNumberPrecision(5);
  }
  plot.legend->setFont(font(13)); plot.legend->setBorderPen(Qt::NoPen); plot.legend->setBrush(Qt::white);
  plot.axisRect()->insetLayout()->take(plot.legend);
  plot.plotLayout()->insertRow(0); plot.plotLayout()->addElement(0, 0, plot.legend);
  plot.legend->setFillOrder(QCPLayoutGrid::foColumnsFirst); plot.legend->setWrap(4);
  plot.plotLayout()->setRowStretchFactor(0, .001);
}
void title(QCustomPlot &plot, const QString &name) {
  plot.plotLayout()->insertRow(0); plot.plotLayout()->addElement(0, 0, new QCPTextElement(&plot, name, font(17, true)));
  plot.plotLayout()->setRowStretchFactor(0, .001);
}
const std::array<std::vector<RealDataType>, 3> familyTypes{{
  {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy},
  {RealTzx, ImagTzx, RealTzy, ImagTzy}, {PTxx, PTxy, PTyx, PTyy}}};
class LocationMap {
public:
  QCustomPlot plot;
  std::map<std::string, QPointF> positions;
  QString description;
  QCPRange x{-1., 1.}, y{-1., 1.};
  LocationMap(const MTSurveyData &survey) {
    style(plot); title(plot, tr("Station map"));
    std::vector<std::string> names;
    std::vector<std::array<double, 3>> locations;
    for(const auto &name: survey.get_stations_names()) {
      auto p = survey.get_station_data(name).position();
      if(!std::isfinite(p[0]) || !std::isfinite(p[1]) || p[0] < -80. || p[0] > 84. || std::abs(p[1]) > 180.) continue;
      if(!std::isfinite(p[2])) p[2] = 0.; // Elevation is not used by the horizontal map.
      names.push_back(name); locations.push_back(p);
    }
    try {
      if(!locations.empty()) {
        auto coordinates = survey.coordinates();
        if(!coordinates.utm) {
          const auto suggestion = SurveyCoordinates::suggested(locations);
          coordinates = SurveyCoordinates::calculate(locations, suggestion.zone, suggestion.north, true);
        }
        description = tr("WGS84 / UTM %1%2").arg(coordinates.zone).arg(coordinates.north ? "N" : "S");
        plot.xAxis->setLabel(coordinates.centered ? tr("East offset (km)") : tr("Easting (km)"));
        plot.yAxis->setLabel(coordinates.centered ? tr("North offset (km)") : tr("Northing (km)"));
        for(unsigned i = 0; i < locations.size(); ++i) {
          try {
            const auto projected = coordinates.transform({locations[i]}).front();
            positions[names[i]] = QPointF(projected[1] / 1000., projected[0] / 1000.);
          } catch(const std::exception &) { /* Keep other valid locations in the report. */ }
        }
      }
    } catch(const std::exception &) { description = tr("Map projection unavailable"); }
    for(int group = 0; group < 3; ++group) {
      auto *g = plot.addGraph(); g->setLineStyle(QCPGraph::lsNone);
      const QColor color(group == 2 ? "#dc2626" : (group == 0 ? "#64748b" : "#cbd5e1"));
      g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, color, group == 2 ? 9 : 5));
      g->setName(group == 2 ? tr("This station") : (group == 0 ? tr("Enabled") : tr("Disabled")));
      QVector<double> px, py;
      if(group < 2) for(const auto &p: positions) if(survey.is_active(p.first) == (group == 0)) { px.push_back(p.second.x()); py.push_back(p.second.y()); }
      g->setData(px, py);
    }
    if(!positions.empty()) {
      x = QCPRange(positions.begin()->second.x(), positions.begin()->second.x());
      y = QCPRange(positions.begin()->second.y(), positions.begin()->second.y());
      for(const auto &p: positions) { x.lower = std::min(x.lower, p.second.x()); x.upper = std::max(x.upper, p.second.x()); y.lower = std::min(y.lower, p.second.y()); y.upper = std::max(y.upper, p.second.y()); }
      const double margin = .1 * std::max({x.size(), y.size(), 1.});
      x = QCPRange(x.lower - margin, x.upper + margin); y = QCPRange(y.lower - margin, y.upper + margin);
    }
    QObject::connect(&plot, &QCustomPlot::afterLayout, &plot, [this] {
      const double w = plot.axisRect()->width(), h = plot.axisRect()->height();
      if(w <= 0. || h <= 0.) return;
      if(plot.xAxis->range().size() / w > plot.yAxis->range().size() / h) plot.yAxis->setScaleRatio(plot.xAxis);
      else plot.xAxis->setScaleRatio(plot.yAxis);
    });
  }
  void render(QCPPainter &p, const QRect &rect, const std::string &station = {}) {
    plot.clearItems(); plot.graph(2)->data()->clear();
    const auto found = positions.find(station);
    plot.legend->setVisible(true);
    plot.legend->item(2)->setVisible(found != positions.end());
    if(found != positions.end()) plot.graph(2)->addData(found->second.x(), found->second.y());
    if(positions.empty()) {
      auto *label = new QCPItemText(&plot); label->position->setType(QCPItemPosition::ptAxisRectRatio);
      label->position->setCoords(.5, .5); label->setText(tr("No usable coordinates")); label->setFont(font(15));
    }
    plot.xAxis->setRange(x); plot.yAxis->setRange(y); plotAt(p, plot, rect);
  }
};
void countTable(QCPPainter &p, double x, double y, double width, const std::array<std::array<unsigned, 4>, 3> &counts, bool aggregate = false) {
  const double first = width * .38, step = width * .155;
  text(p, QRectF(x, y, first, 24), aggregate ? tr("Station-periods") : tr("Periods"), 13, true);
  for(int i = 0; i < 4; ++i) text(p, QRectF(x + first + i * step, y, step, 24), QStringList{tr("Full"), tr("Part."), tr("Mask"), tr("Miss.")}[i], 13, true);
  for(unsigned row = 0; row < 3; ++row) {
    const double top = y + 30 + row * 27;
    text(p, QRectF(x, top, first, 24), QStringList{tr("Impedance"), tr("Tipper"), tr("Phase tensor")}[row], 13);
    for(int i = 0; i < 4; ++i) text(p, QRectF(x + first + i * step, top, step, 24), QString::number(counts[row][i]), 13);
  }
}
}

namespace SurveyReport {
PeriodSummary summarize_periods(const MTStationData &station)
{
  PeriodSummary info;
  auto raw = station; raw.set_active(true);
  for(double f: raw.frequencies()) for(const auto &types: familyTypes) for(auto type: types) raw.set_data_mask(type, f, true);
  for(unsigned f = 0; f < station.frequencies().size(); ++f) {
    const double period = 1. / station.frequencies()[f];
    if(!std::isfinite(period) || period <= 0.) continue;
    info.minimum = std::min(info.minimum, period); info.maximum = std::max(info.maximum, period);
    auto &states = info.periods[period];
    for(unsigned family = 0; family < familyTypes.size(); ++family) {
      unsigned usable = 0, stored = 0;
      double value, error;
      for(auto type: familyTypes[family]) {
        usable += station.scalar_value(type, f, value, error);
        stored += raw.scalar_value(type, f, value, error);
      }
      const auto state = usable == familyTypes[family].size() ? PeriodState::Full :
        (usable ? PeriodState::Partial : (stored ? PeriodState::Masked : PeriodState::Missing));
      states[family] = state; ++info.counts[family][static_cast<unsigned>(state)];
    }
  }
  return info;
}

bool write_pdf(const QString &path, const MTSurveyData &survey, const MTSurveyData *response,
               const Options &options, const std::function<bool(unsigned, unsigned)> &progress)
{
  std::vector<std::string> included;
  for(const auto &name: survey.get_stations_names()) if(options.includeDisabled || survey.is_active(name)) included.push_back(name);
  if(included.empty()) throw std::runtime_error("No stations match the report selection.");
  const unsigned pageCount = included.size() + 1;
  if(progress && !progress(0, pageCount)) return false;
  std::map<std::string, PeriodSummary> information;
  std::map<double, std::array<unsigned, 3>> coverage;
  std::array<std::array<unsigned, 4>, 3> totals{};
  unsigned storedPeriods = 0;
  unsigned enabled = 0;
  double minimum = std::numeric_limits<double>::infinity(), maximum = 0.;
  for(const auto &name: included) {
    const auto &s = survey.get_station_data(name); auto &info = information[name]; info = summarize_periods(s);
    enabled += s.active(); minimum = std::min(minimum, info.minimum); maximum = std::max(maximum, info.maximum);
    storedPeriods += info.periods.size();
    for(const auto &period: info.periods) {
      auto &counts = coverage[period.first];
      for(unsigned c = 0; c < 3; ++c) if(period.second[c] == PeriodState::Full || period.second[c] == PeriodState::Partial) ++counts[c];
    }
    for(unsigned c = 0; c < 3; ++c) for(unsigned state = 0; state < 4; ++state) totals[c][state] += info.counts[c][state];
  }
  FitStatistics::Result fit;
  if(response) fit = FitStatistics::compare(survey, *response, FitStatistics::ErrorSource::Observed);
  LocationMap map(survey);
  QSaveFile file(path);
  if(!file.open(QIODevice::WriteOnly)) throw std::runtime_error(file.errorString().toStdString());
  bool canceled = false;
  {
    QPdfWriter pdf(&file); pdf.setTitle(options.title); pdf.setCreator("EDITools");
    pdf.setResolution(144);
    pdf.setPageLayout(QPageLayout(QPageSize(QPageSize::A4), QPageLayout::Landscape, QMarginsF(10, 10, 10, 10)));
    QCPPainter painter;
    if(!painter.begin(&pdf)) throw std::runtime_error("Cannot initialize the PDF writer.");
    const int w = 1400, h = qRound(double(w) * pdf.height() / pdf.width());
    auto beginPage = [&](const QString &heading, const QString &subheading, unsigned page) {
      painter.resetTransform(); painter.scale(double(pdf.width()) / w, double(pdf.width()) / w);
      painter.fillRect(QRect(0, 0, w, h), Qt::white);
      painter.setFont(font(27, true));
      text(painter, QRectF(0, 0, w, 40), painter.fontMetrics().elidedText(heading, Qt::ElideRight, w), 27, true);
      text(painter, QRectF(0, 45, w, 40), subheading, 15); line(painter, 88, w);
      line(painter, h - 33, w);
      painter.setFont(font(12));
      text(painter, QRectF(0, h - 24, w - 180, 22), painter.fontMetrics().elidedText(options.title, Qt::ElideRight, w - 180), 12);
      text(painter, QRectF(w - 145, h - 24, 145, 22), tr("Page %1 of %2").arg(page).arg(pageCount), 12);
    };
    beginPage(options.title, tr("Survey overview · %1").arg(QDateTime::currentDateTime().toString("yyyy-MM-dd HH:mm")), 1);
    const int left = 690;
    text(painter, QRectF(0, 106, left, 36), tr("%1 stations · %2 enabled · %3 disabled").arg(included.size()).arg(enabled).arg(included.size() - enabled), 21, true);
    text(painter, QRectF(0, 150, left, 66), tr("Period range: %1 – %2 s\nScope: %3 of %4 survey stations · %5 mapped locations")
         .arg(number(minimum)).arg(maximum > 0. ? number(maximum) : tr("N/A")).arg(included.size()).arg(survey.n_stations()).arg(map.positions.size()), 17);
    text(painter, QRectF(0, 217, left, 22), tr("%1 distinct periods · %2 stored station-periods").arg(coverage.size()).arg(storedPeriods), 14);
    countTable(painter, 0, 249, left, totals, true);
    text(painter, QRectF(0, 365, left, 78), tr("Each stored period counts once per data type. Full: all components usable; Part.: some. Mask: finite values exist but none are enabled; Miss.: no finite values. Disabled data count as masked. Counts ignore error quality and plot visibility."), 14);
    int noteY = 450;
    if(response) {
      painter.setFont(font(16));
      text(painter, QRectF(0, noteY, left, 55), tr("Response: %1\nnRMS: %2 · %3 matched scalars (observed errors)")
           .arg(painter.fontMetrics().elidedText(options.responseName, Qt::ElideMiddle, left - 100))
           .arg(fit.total.count ? number(fit.total.rms()) : tr("N/A")).arg(fit.total.count), 16);
      text(painter, QRectF(0, noteY + 48, left, 24), tr("Fit uses all enabled components, independent of plot visibility."), 13);
      noteY += 80;
    }
    if(!survey.resampling_info().method.empty()) {
      text(painter, QRectF(0, noteY, left, 67), tr("Resampled survey: %1 target periods; maximum gap ratio %2. Estimates share source information; coverage does not represent new measurements.")
           .arg(survey.resampling_info().periods.size()).arg(survey.resampling_info().maximumRatio), 15);
      noteY += 74;
    }
    QCustomPlot coveragePlot; style(coveragePlot); title(coveragePlot, tr("Period coverage · full + partial")); coveragePlot.legend->setVisible(true);
    coveragePlot.xAxis->setScaleType(QCPAxis::stLogarithmic); coveragePlot.xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    coveragePlot.xAxis->setLabel(tr("Period (s)")); coveragePlot.yAxis->setLabel(tr("Enabled stations with data"));
    const QStringList coverageNames{tr("Impedance"), tr("Tipper"), tr("Phase tensor")};
    const QColor colors[] = {QColor("#2563eb"), QColor("#d97706"), QColor("#9333ea")};
    for(unsigned c = 0; c < 3; ++c) {
      auto *curve = coveragePlot.addGraph(); QVector<double> cx, cy;
      for(const auto &point: coverage) { cx.push_back(point.first); cy.push_back(point.second[c]); }
      curve->setName(coverageNames[c]); curve->setData(cx, cy);
      curve->setPen(QPen(colors[c], 1.5, c == 2 ? Qt::DashLine : Qt::SolidLine));
      curve->setScatterStyle(QCPScatterStyle(c == 2 ? QCPScatterStyle::ssCircle : QCPScatterStyle::ssDisc, colors[c], 4));
    }
    if(!coverage.empty()) coveragePlot.xAxis->setRange(coverage.begin()->first / 1.1, coverage.rbegin()->first * 1.1);
    else coveragePlot.xAxis->setRange(.1, 10.);
    coveragePlot.yAxis->setRange(0., std::max(1., enabled * 1.05));
    plotAt(painter, coveragePlot, QRect(0, noteY, left, h - noteY - 55));
    map.render(painter, QRect(735, 108, w - 735, h - 295));
    text(painter, QRectF(755, h - 175, w - 755, 60), map.description + tr("\nAll survey locations shown; missing coordinates are omitted."), 15);
    painter.setFont(font(14));
    text(painter, QRectF(755, h - 100, w - 755, 45), tr("Source: %1").arg(painter.fontMetrics().elidedText(options.source.isEmpty() ? tr("Unsaved survey") : options.source, Qt::ElideMiddle, w - 805)), 14);
    if(progress && !progress(1, pageCount)) canceled = true;

    std::array<std::unique_ptr<QCustomPlot>, 4> plots;
    std::array<std::unique_ptr<MTDataPlot>, 4> handlers;
    const QStringList plotNames{tr("Apparent resistivity"), tr("Phase"), tr("Tipper"), tr("Phase tensor")};
    for(unsigned i = 0; i < 4; ++i) {
      plots[i].reset(new QCustomPlot);
      plots[i]->setObjectName(QString("surveyReportPlot%1").arg(i));
      if(i == 0) handlers[i].reset(new ApparentResistivityPlot(plots[i].get()));
      if(i == 1) { auto *phase = new PhasePlot(plots[i].get()); phase->set_phase_wrap(options.phaseWrap); handlers[i].reset(phase); }
      if(i == 2) handlers[i].reset(new TipperPlot(plots[i].get()));
      if(i == 3) handlers[i].reset(new PhaseTensorPlot(plots[i].get()));
      handlers[i]->set_error_bars_visible(options.errorBars); handlers[i]->set_component_visibility(options.components[i]);
      handlers[i]->set_y_axis_range(options.axes[i].lower, options.axes[i].upper);
      handlers[i]->set_y_axis_autoscale(options.axes[i].autoscale);
      // Static paper legends omit interactive controls while using the same data/curve rendering.
      plots[i]->legend->clearItems();
      for(unsigned c = 0; c < 4; ++c) if(options.components[i][c]) plots[i]->legend->addItem(new QCPPlottableLegendItem(plots[i]->legend, plots[i]->graph(c)));
      style(*plots[i]); title(*plots[i], plotNames[i]);
    }
    for(unsigned page = 0; page < included.size() && !canceled; ++page) {
      if(!pdf.newPage()) throw std::runtime_error("Cannot append a PDF page.");
      const auto &name = included[page]; auto station = survey.get_station_data(name);
      const auto &info = information.at(name);
      beginPage(QString::fromStdString(name), tr("%1 · %2 periods · %3 – %4 s")
                .arg(station.active() ? tr("Enabled") : tr("Disabled")).arg(station.frequencies().size()).arg(number(info.minimum)).arg(info.maximum > 0. ? number(info.maximum) : tr("N/A")), page + 2);
      const MTStationData *prediction = response && response->is_station_present(name) ? &response->get_station_data(name) : nullptr;
      const int plotWidth = 510, plotHeight = (h - 155) / 2;
      for(unsigned i = 0; i < 4; ++i) {
        auto &plot = *plots[i]; plot.clearItems(); handlers[i]->clear_predicted_data();
        // Reset empty panels so a previous station's scales cannot leak into this page.
        plot.xAxis->setRange(std::isfinite(info.minimum) ? info.minimum / 1.5 : .1, info.maximum > 0. ? info.maximum * 1.5 : 10.);
        plot.yAxis->setRange(i == 0 ? QCPRange(.1, 1000.) : (i == 1 ? QCPRange(-180., 180.) : QCPRange(-1., 1.)));
        handlers[i]->set_observed_data(station);
        if(prediction) handlers[i]->set_predicted_data(*prediction, true);
        bool hasData = false;
        for(int g = 0; g < plot.graphCount(); ++g) if(plot.graph(g)->visible())
          for(const auto &value: *plot.graph(g)->data()) if(std::isfinite(value.value) && (i != 0 || value.value > 0.)) { hasData = true; break; }
        if(!hasData) {
          auto *label = new QCPItemText(&plot); label->position->setType(QCPItemPosition::ptAxisRectRatio); label->position->setCoords(.5, .5);
          label->setText(tr("No displayed data")); label->setFont(font(16)); label->setColor(QColor("#64748b"));
        }
        plotAt(painter, plot, QRect((i % 2) * (plotWidth + 12), 100 + (i / 2) * (plotHeight + 8), plotWidth, plotHeight));
      }
      const int sideX = 1060, sideW = w - sideX;
      map.render(painter, QRect(sideX, 100, sideW, 330), name);
      text(painter, QRectF(sideX, 440, sideW, 50), map.positions.count(name) ? map.description : tr("Station location unavailable"), 14);
      const auto pos = station.position();
      text(painter, QRectF(sideX, 492, sideW, 86), tr("Latitude: %1°\nLongitude: %2°\nElevation: %3 m")
           .arg(number(pos[0], 9)).arg(number(pos[1], 9)).arg(number(pos[2], 7)), 16);
      countTable(painter, sideX, 596, sideW, info.counts);
      text(painter, QRectF(sideX, 706, sideW, 34), tr("Full / partial / masked / missing periods\nEach row sums to this station's stored periods."), 11);
      if(response) {
        const auto match = fit.stations.find(name);
        text(painter, QRectF(sideX, 757, sideW, 62), prediction ? tr("Response overlay\nnRMS: %1 · N: %2")
             .arg(match != fit.stations.end() && match->second.count ? number(match->second.rms()) : tr("N/A"))
             .arg(match != fit.stations.end() ? match->second.count : 0) : tr("No response for this station"), 16);
      }
      const bool fixed = std::any_of(options.axes.begin(), options.axes.end(), [](const Options::Axis &axis) { return !axis.autoscale; });
      text(painter, QRectF(sideX, h - 126, sideW, 80), tr("Points: observations\nGrey points: masked / disabled\nLines: computed response\n") +
           (fixed ? tr("Y: GUI fixed ranges where set\nOther axes: auto per station") : tr("Axes auto-scaled per station")), 13);
      if(progress && !progress(page + 2, pageCount)) canceled = true;
    }
    if(!painter.end()) throw std::runtime_error("Cannot finish the PDF report.");
  }
  if(canceled) { file.cancelWriting(); return false; }
  if(!file.commit()) throw std::runtime_error(file.errorString().toStdString());
  return true;
}

void show_dialog(QWidget *parent, const MTSurveyData &survey, const std::map<std::string, MTSurveyData> &responses,
                 const QString &selectedResponse, const Options &defaults, QString &lastDirectory)
{
  if(!survey.n_stations()) { QMessageBox::information(parent, tr("Survey report"), tr("Load a survey first.")); return; }
  QDialog dialog(parent); dialog.setObjectName("surveyReportDialog"); dialog.setWindowTitle(tr("Export survey report"));
  auto *layout = new QVBoxLayout(&dialog); auto *form = new QFormLayout;
  auto *title = new QLineEdit(defaults.title, &dialog); title->setObjectName("reportTitle"); form->addRow(tr("Title:"), title);
  auto *scope = new QComboBox(&dialog); scope->setObjectName("reportStations"); scope->addItems({tr("All stations"), tr("Enabled only")}); form->addRow(tr("Stations:"), scope);
  auto *response = new QComboBox(&dialog); response->setObjectName("reportResponse"); response->addItem(tr("None"), QString());
  for(const auto &entry: responses) { const QString path = QString::fromStdString(entry.first); response->addItem(QFileInfo(path).fileName(), path); response->setItemData(response->count() - 1, path, Qt::ToolTipRole); }
  response->setCurrentIndex(std::max(0, response->findData(selectedResponse))); form->addRow(tr("Response:"), response);
  auto *errors = new QCheckBox(tr("Error bars"), &dialog); errors->setObjectName("reportErrorBars"); errors->setChecked(defaults.errorBars); form->addRow(errors);
  auto *visible = new QCheckBox(tr("Visible components only"), &dialog); visible->setObjectName("reportVisibleOnly"); form->addRow(visible);
  layout->addLayout(form);
  layout->addWidget(UiHelp::label(&dialog, tr("Report"), tr("A landscape PDF with a survey overview and one page per included station. Each station has resistivity, phase, tipper and phase-tensor plots beside a location map. All components are included unless Visible components only is checked. Current phase wrapping and fixed GUI Y ranges are retained. Y axes set to autoscale and period axes scale to each station. Optional response curves use observed errors for nRMS. Maps show all valid survey locations and highlight the station; missing locations are noted. Period tables count each stored period once per data type: Full = all components usable; Part. = some usable; Mask = finite data but none enabled; Miss. = no finite data. Impedance needs four complex components, tipper two, and phase tensor four entries for Full. Counts ignore error quality and plot visibility. Coverage shows stations with full or partial data at each exact period. Export does not change your survey or current plots."), "reportHelp"));
  auto *pages = new QLabel(&dialog); pages->setObjectName("reportPageCount"); layout->addWidget(pages);
  auto updatePages = [&] {
    unsigned count = 0; for(const auto &name: survey.get_stations_names()) if(scope->currentIndex() == 0 || survey.is_active(name)) ++count;
    pages->setText(tr("%1 station pages + 1 overview").arg(count));
  };
  QObject::connect(scope, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, updatePages); updatePages();
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel, &dialog); layout->addWidget(buttons);
  buttons->button(QDialogButtonBox::Save)->setText(tr("Export…"));
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if(dialog.exec() != QDialog::Accepted) return;
  QFileDialog destination(parent, tr("Export survey report"), lastDirectory + "/survey_report.pdf", tr("PDF (*.pdf)"));
  destination.setAcceptMode(QFileDialog::AcceptSave);
  destination.setDefaultSuffix("pdf"); // Resolve the suffix before Qt checks whether the destination exists.
  if(destination.exec() != QDialog::Accepted || destination.selectedFiles().isEmpty()) return;
  const auto path = destination.selectedFiles().first();
  Options options = defaults; options.title = title->text().trimmed(); if(options.title.isEmpty()) options.title = tr("Survey data report");
  options.includeDisabled = scope->currentIndex() == 0; options.errorBars = errors->isChecked();
  if(!visible->isChecked()) for(auto &group: options.components) group.fill(true);
  const auto found = responses.find(response->currentData().toString().toStdString());
  options.responseName = found == responses.end() ? QString() : response->currentText();
  const auto snapshot = survey;
  const auto prediction = found == responses.end() ? std::shared_ptr<MTSurveyData>() : std::make_shared<MTSurveyData>(found->second);
  QProgressDialog progress(tr("Writing survey report…"), tr("Cancel"), 0, 0, parent); progress.setWindowModality(Qt::ApplicationModal); progress.setMinimumDuration(0);
  progress.setAutoClose(false); progress.setAutoReset(false);
  try {
    const bool saved = write_pdf(path, snapshot, prediction.get(), options, [&](unsigned done, unsigned total) {
      progress.setRange(0, total); progress.setValue(done); QApplication::processEvents(); return !progress.wasCanceled();
    });
    progress.close();
    if(saved) lastDirectory = QFileInfo(path).absolutePath();
  } catch(const std::exception &error) { progress.close(); QMessageBox::warning(parent, tr("Survey report"), QString::fromUtf8(error.what())); }
}
}
