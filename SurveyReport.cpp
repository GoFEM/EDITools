#include "include/MapBackground.h"
#include "include/FileLabels.h"
#include "include/SurveyReport.h"
#include "PeriodMapWindow.h"
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
#include <QGridLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QPdfWriter>
#include <QProgressDialog>
#include <QPushButton>
#include <QSaveFile>
#include <QRegularExpression>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace {
using namespace SurveyReport;
QString tr(const char *text) { return QObject::tr(text); }
QFont font(int pixels, bool bold = false) { QFont f("Sans Serif"); f.setPixelSize(pixels); f.setBold(bold); return f; }
QString number(double value, int precision = 6) { return std::isfinite(value) ? QString::number(value, 'g', precision) : tr("N/A"); }
std::vector<double> validatedMapPeriods(std::vector<double> periods) {
  if(periods.empty()) throw std::invalid_argument("Choose at least one map period.");
  for(double p: periods) if(!std::isfinite(p) || p <= 0.) throw std::invalid_argument("Map periods must be positive numbers in seconds.");
  std::sort(periods.begin(), periods.end()); periods.erase(std::unique(periods.begin(), periods.end()), periods.end()); return periods;
}
std::vector<double> parseMapPeriods(const QString &text) {
  std::vector<double> result;
  for(const auto &token: text.split(QRegularExpression("[,;\\s]+"), Qt::SkipEmptyParts)) {
    bool ok = false; const double value = token.toDouble(&ok);
    if(!ok) throw std::invalid_argument("Map periods must be positive numbers in seconds.");
    result.push_back(value);
  }
  return validatedMapPeriods(result);
}
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
  LocationMap(const MTSurveyData &survey, int layers) {
    auto *background = new MapBackground(&plot); background->setLayers(layers);
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
        background->setCoordinates(coordinates, .001);
        description = tr("WGS84 / UTM %1%2").arg(coordinates.zone).arg(coordinates.north ? "N" : "S");
        description += tr("\nUTM origin E: %1 m\nUTM origin N: %2 m")
          .arg(coordinates.origin_easting, 0, 'f', 2).arg(coordinates.origin_northing, 0, 'f', 2);
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
  void render(QCPPainter &p, const QRect &rect, const std::string &station = {}, bool showNames = false) {
    plot.clearItems(); plot.graph(2)->data()->clear();
    const auto found = positions.find(station);
    plot.legend->setVisible(true);
    plot.legend->item(2)->setVisible(found != positions.end());
    if(found != positions.end()) plot.graph(2)->addData(found->second.x(), found->second.y());
    if(showNames) for(const auto &entry: positions) {
      auto *label = new QCPItemText(&plot);
      label->position->setType(QCPItemPosition::ptPlotCoords);
      label->position->setCoords(entry.second);
      label->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom);
      label->setPadding(QMargins(5, 2, 2, 3));
      label->setText(QString::fromStdString(entry.first));
      label->setFont(font(12)); label->setColor(QColor("#334155"));
      label->setBrush(QColor(255, 255, 255, 210));
    }
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
  const bool maps = options.phaseTensorMaps || options.inductionMaps;
  const auto mapPeriods = maps ? validatedMapPeriods(options.mapPeriods) : std::vector<double>{};
  if(maps && options.mapData != Options::MapData::Observed && !response) throw std::invalid_argument("Select a response for response maps.");
  if(options.inductionMaps && !options.mapSettings.real && !options.mapSettings.imaginary) throw std::invalid_argument("Choose real or imaginary induction vectors.");
  const unsigned pageCount = unsigned(options.overview) + (options.stationPages ? included.size() : 0) +
    mapPeriods.size() * (options.mapData == Options::MapData::Both ? 2 : 1);
  if(!pageCount) throw std::invalid_argument("Select at least one report section.");
  if(progress && !progress(0, pageCount)) return false;
  // Response labels may come from projects saved on another platform. Never
  // print directories, drive letters or network host names in the report.
  auto responseName = FileLabels::fileName(options.responseName);
  if(responseName.isEmpty()) responseName = tr("Computed response");
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
  LocationMap map(survey, options.mapLayers);
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
    unsigned writtenPages = 0;
    auto beginPage = [&](const QString &heading, const QString &subheading) {
      if(writtenPages && !pdf.newPage()) throw std::runtime_error("Cannot append a PDF page.");
      ++writtenPages;
      painter.resetTransform(); painter.scale(double(pdf.width()) / w, double(pdf.width()) / w);
      painter.fillRect(QRect(0, 0, w, h), Qt::white);
      painter.setFont(font(27, true));
      text(painter, QRectF(0, 0, w, 40), painter.fontMetrics().elidedText(heading, Qt::ElideRight, w), 27, true);
      text(painter, QRectF(0, 45, w, 40), subheading, 15); line(painter, 88, w);
      line(painter, h - 33, w);
      painter.setFont(font(12));
      text(painter, QRectF(0, h - 24, w - 180, 22), painter.fontMetrics().elidedText(options.title, Qt::ElideRight, w - 180), 12);
      text(painter, QRectF(w - 145, h - 24, 145, 22), tr("Page %1 of %2").arg(writtenPages).arg(pageCount), 12);
    };
    if(options.overview) {
    beginPage(options.title, tr("Survey overview · %1").arg(QDateTime::currentDateTime().toString("yyyy-MM-dd HH:mm")));
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
           .arg(painter.fontMetrics().elidedText(responseName, Qt::ElideMiddle, left - 100))
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
    map.render(painter, QRect(735, 108, w - 735, h - 315), {}, true);
    text(painter, QRectF(755, h - 190, w - 755, 85), map.description + tr("\nAll survey locations shown; missing coordinates are omitted."), 15);
    if(progress && !progress(writtenPages, pageCount)) canceled = true;
    }

    if(options.stationPages && !canceled) {
    std::array<std::unique_ptr<QCustomPlot>, 4> plots;
    std::array<std::unique_ptr<MTDataPlot>, 4> handlers;
    const QStringList plotNames{tr("Apparent resistivity"), tr("Phase"), tr("Tipper"), tr("Phase tensor")};
    for(unsigned i = 0; i < 4; ++i) {
      plots[i].reset(new QCustomPlot);
      plots[i]->setObjectName(QString("surveyReportPlot%1").arg(i));
      if(i == 0) handlers[i].reset(new ApparentResistivityPlot(plots[i].get()));
      if(i == 1) { auto *phase = new PhasePlot(plots[i].get()); phase->set_phase_wrap(options.phaseWrap); handlers[i].reset(phase); }
      if(i == 2) { auto *tipper = new TipperPlot(plots[i].get()); tipper->set_arrow_mode(options.tipperArrows); handlers[i].reset(tipper); }
      if(i == 3) handlers[i].reset(new PhaseTensorPlot(plots[i].get()));
      handlers[i]->set_error_bars_visible(options.errorBars); handlers[i]->set_component_visibility(options.components[i]);
      handlers[i]->set_y_axis_range(options.axes[i].lower, options.axes[i].upper);
      handlers[i]->set_y_axis_autoscale(options.axes[i].autoscale);
      // Static paper legends omit interactive controls while using the same data/curve rendering.
      plots[i]->legend->clearItems();
      for(unsigned c = 0; c < 4; ++c) {
        const bool shown = i == 2 && options.tipperArrows ? (c % 2 == 0 && (options.components[i][c] || options.components[i][c + 1])) : options.components[i][c];
        if(shown) plots[i]->legend->addItem(new QCPPlottableLegendItem(plots[i]->legend, plots[i]->graph(c)));
      }
      style(*plots[i]); title(*plots[i], i == 2 && options.tipperArrows ? tr("Tipper arrows") : plotNames[i]);
    }
    for(unsigned page = 0; page < included.size() && !canceled; ++page) {
      const auto &name = included[page]; auto station = survey.get_station_data(name);
      const auto &info = information.at(name);
      beginPage(QString::fromStdString(name), tr("%1 · %2 periods · %3 – %4 s")
                .arg(station.active() ? tr("Enabled") : tr("Disabled")).arg(station.frequencies().size()).arg(number(info.minimum)).arg(info.maximum > 0. ? number(info.maximum) : tr("N/A")));
      const MTStationData *prediction = response && response->is_station_present(name) ? &response->get_station_data(name) : nullptr;
      const int plotWidth = 510, plotHeight = (h - 155) / 2;
      for(unsigned i = 0; i < 4; ++i) {
        auto &plot = *plots[i];
        for(int item = plot.itemCount() - 1; item >= 0; --item) if(plot.item(item)->objectName() == "reportNoData") plot.removeItem(item);
        handlers[i]->clear_predicted_data();
        // Reset empty panels so a previous station's scales cannot leak into this page.
        plot.xAxis->setRange(std::isfinite(info.minimum) ? info.minimum / 1.5 : .1, info.maximum > 0. ? info.maximum * 1.5 : 10.);
        plot.yAxis->setRange(i == 0 ? QCPRange(.1, 1000.) : (i == 1 ? QCPRange(-180., 180.) : QCPRange(-1., 1.)));
        handlers[i]->set_observed_data(station);
        if(prediction) handlers[i]->set_predicted_data(*prediction, true);
        if(i == 2 && !options.tipperArrows) {
          const auto mask = station.tipper_mask();
          std::vector<dvector> values, errors; station.get_tipper(values, errors);
          for(unsigned c = 0; c < 4; ++c) {
            // Include gaps in both the line and its error vector, so inserting
            // a masked period cannot shift subsequent uncertainty bars.
            std::map<double, std::pair<double, double>> samples;
            for(unsigned f = 0; f < station.frequencies().size(); ++f)
              samples[1. / station.frequencies()[f]] = mask[c][f] ? std::make_pair(values[c][f], errors[c][f]) : std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
            QVector<double> x, y, uncertainty;
            for(const auto &sample: samples) { x << sample.first; y << sample.second.first; uncertainty << sample.second.second; }
            plot.graph(c)->setData(x, y, true); plot.graph(c)->setLineStyle(QCPGraph::lsLine);
            for(int index = 0; index < plot.plottableCount(); ++index) {
              auto *bars = dynamic_cast<QCPErrorBars *>(plot.plottable(index));
              if(bars && bars->dataPlottable() == plot.graph(c)) bars->setData(uncertainty);
            }
          }
        }
        bool hasData = false;
        for(int g = 0; g < plot.graphCount(); ++g) if(plot.graph(g)->visible())
          for(const auto &value: *plot.graph(g)->data()) if(std::isfinite(value.value) && (i != 0 || value.value > 0.)) { hasData = true; break; }
        if(i == 2 && options.tipperArrows) {
          hasData = false;
          for(const MTStationData *data: std::array<const MTStationData *, 2>{{&station, prediction}}) if(data) {
            std::vector<dvector> values, errors; data->get_tipper(values, errors);
            for(unsigned offset: {0u, 2u}) if(options.components[i][offset] || options.components[i][offset + 1])
              for(unsigned f = 0; f < data->frequencies().size(); ++f)
                hasData |= (!options.components[i][offset] || std::isfinite(values[offset][f])) &&
                           (!options.components[i][offset + 1] || std::isfinite(values[offset + 1][f]));
          }
        }
        if(!hasData) {
          auto *label = new QCPItemText(&plot); label->setObjectName("reportNoData"); label->position->setType(QCPItemPosition::ptAxisRectRatio); label->position->setCoords(.5, .5);
          label->setText(tr("No displayed data")); label->setFont(font(16)); label->setColor(QColor("#64748b"));
        }
        QMetaObject::Connection arrowLayout;
        if(i == 2 && options.tipperArrows) {
          for(int item = 0; item < plot.itemCount(); ++item) {
            if(auto *label = dynamic_cast<QCPItemText *>(plot.item(item))) if(label->text().startsWith("|T|")) label->setFont(font(13));
            if(auto *arrow = dynamic_cast<QCPItemLine *>(plot.item(item)))
              if(arrow->start->type() == QCPItemPosition::ptAxisRectRatio) arrow->start->setCoords(.6, .12);
          }
          arrowLayout = QObject::connect(&plot, &QCustomPlot::afterLayout, &plot, [&plot] {
            // Scale all vectors and their reference together to fit the printed
            // panel. PDF panels can be much narrower than the interactive plot.
            double scale = 1.; const auto bounds = plot.axisRect()->rect().adjusted(5, 5, -5, -5);
            for(int item = 0; item < plot.itemCount(); ++item) if(auto *arrow = dynamic_cast<QCPItemLine *>(plot.item(item))) {
              if(arrow->start->type() != QCPItemPosition::ptPlotCoords) continue;
              const auto start = arrow->start->pixelPosition(), delta = arrow->end->coords();
              if(delta.x() > 0.) scale = std::min(scale, (bounds.right() - start.x()) / delta.x());
              if(delta.x() < 0.) scale = std::min(scale, (bounds.left() - start.x()) / delta.x());
              if(delta.y() > 0.) scale = std::min(scale, (bounds.bottom() - start.y()) / delta.y());
              if(delta.y() < 0.) scale = std::min(scale, (bounds.top() - start.y()) / delta.y());
            }
            if(scale > 0. && scale < 1.) for(int item = 0; item < plot.itemCount(); ++item)
              if(auto *arrow = dynamic_cast<QCPItemLine *>(plot.item(item))) arrow->end->setCoords(arrow->end->coords() * scale);
          });
        }
        plotAt(painter, plot, QRect((i % 2) * (plotWidth + 12), 100 + (i / 2) * (plotHeight + 8), plotWidth, plotHeight));
        if(arrowLayout) QObject::disconnect(arrowLayout);
      }
      const int sideX = 1060, sideW = w - sideX;
      map.render(painter, QRect(sideX, 100, sideW, 330), name);
      text(painter, QRectF(sideX, 435, sideW, 66), map.description + (map.positions.count(name) ? QString() : tr("\nStation location unavailable")), 14);
      const auto pos = station.position();
      text(painter, QRectF(sideX, 510, sideW, 72), tr("Latitude: %1°\nLongitude: %2°\nElevation: %3 m")
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
      text(painter, QRectF(sideX, h - 126, sideW, 80), (options.tipperArrows ? tr("Points / solid arrows: observations\nGrey: masked / disabled\nCurves / dashed arrows: response\n") : tr("Points / tipper lines: observations\nGrey points: masked / disabled\nCurves: computed response\n")) +
           (fixed ? tr("Y: GUI fixed ranges where set\nOther axes: auto per station") : tr("Axes auto-scaled per station")), 13);
      if(progress && !progress(writtenPages, pageCount)) canceled = true;
    }
    }
    if(maps && !canceled) {
      auto mapSurvey = std::make_shared<MTSurveyData>(survey);
      if(!options.includeDisabled) for(const auto &name: survey.get_stations_names()) if(!survey.is_active(name)) mapSurvey->remove_station(name);
      const auto responseKey = responseName;
      std::map<std::string, MTSurveyData> mapResponses;
      if(response) mapResponses.emplace(responseKey.toStdString(), *response);
      PeriodMapWindow mapWindow(nullptr);
      mapWindow.setData(mapSurvey, mapResponses);
      auto settings = options.useMapSettings ? options.mapSettings : mapWindow.displaySettings();
      settings.mapLayers = options.mapLayers;
      settings.real = options.mapSettings.real; settings.imaginary = options.mapSettings.imaginary;
      mapWindow.setDisplaySettings(settings);
      QStringList datasets;
      if(options.mapData != Options::MapData::Response) datasets.push_back(QString());
      if(options.mapData != Options::MapData::Observed) datasets.push_back(responseKey);
      for(double seconds: mapPeriods) {
        if(canceled) break;
        for(const auto &key: datasets) {
          if(canceled) break;
          const auto heading = options.phaseTensorMaps && options.inductionMaps ? tr("Phase tensors and induction vectors") :
                               (options.phaseTensorMaps ? tr("Phase-tensor ellipses") : tr("Induction vectors"));
          beginPage(heading, tr("T = %1 s · %2 · Tolerance: %3%")
                    .arg(number(seconds, 12)).arg(key.isEmpty() ? tr("Observed data") : key).arg(settings.tolerancePercent));
          const auto summary = mapWindow.renderReport(painter, QRect(0, 100, w, h - 215), seconds, options.phaseTensorMaps, options.inductionMaps, key);
          text(painter, QRectF(0, h - 105, w, 55), summary + tr("\nMasked, disabled, missing and unmatched values are omitted. Map extent is fitted to this page."), 14);
          if(progress && !progress(writtenPages, pageCount)) canceled = true;
        }
      }
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
  for(const auto &entry: responses) { const QString path = QString::fromStdString(entry.first); response->addItem(FileLabels::fileName(path), path); response->setItemData(response->count() - 1, path, Qt::ToolTipRole); }
  response->setCurrentIndex(std::max(0, response->findData(selectedResponse))); form->addRow(tr("Response:"), response);
  auto *errors = new QCheckBox(tr("Error bars"), &dialog); errors->setObjectName("reportErrorBars"); errors->setChecked(defaults.errorBars); form->addRow(errors);
  auto *visible = new QCheckBox(tr("Visible components only"), &dialog); visible->setObjectName("reportVisibleOnly"); form->addRow(visible);
  auto checkbox = [&](const QString &label, const char *name, bool checked) {
    auto *box = new QCheckBox(label, &dialog); box->setObjectName(name); box->setChecked(checked); return box;
  };
  auto *overview = checkbox(tr("Survey overview"), "reportOverview", defaults.overview);
  auto *stations = checkbox(tr("Station plots"), "reportStationPages", defaults.stationPages);
  auto *ptMaps = checkbox(tr("PT ellipse maps"), "reportPTMaps", defaults.phaseTensorMaps);
  auto *induction = checkbox(tr("Induction maps"), "reportInductionMaps", defaults.inductionMaps);
  auto *sections = new QGridLayout;
  sections->addWidget(overview, 0, 0); sections->addWidget(stations, 0, 1);
  sections->addWidget(ptMaps, 1, 0); sections->addWidget(induction, 1, 1);
  form->addRow(tr("Include:"), sections);
  auto *tipper = new QComboBox(&dialog); tipper->setObjectName("reportTipperStyle");
  tipper->addItems({tr("Lines"), tr("Arrows")}); tipper->setCurrentIndex(defaults.tipperArrows ? 1 : 0);
  form->addRow(tr("Station tippers:"), tipper);
  auto *mapData = new QComboBox(&dialog); mapData->setObjectName("reportMapData");
  mapData->addItems({tr("Observed"), tr("Response"), tr("Both")}); mapData->setCurrentIndex(int(defaults.mapData));
  form->addRow(tr("Map data:"), mapData);
  std::set<double> surveyPeriods;
  for(const auto &name: survey.get_stations_names()) for(double f: survey.get_station_data(name).frequencies())
    if(std::isfinite(f) && f > 0.) surveyPeriods.insert(1. / f);
  QStringList allPeriods; for(double p: surveyPeriods) allPeriods << QString::number(p, 'g', 17);
  double initialPeriod = defaults.mapSettings.period;
  if(!defaults.useMapSettings && !surveyPeriods.empty()) initialPeriod = *std::min_element(surveyPeriods.begin(), surveyPeriods.end(), [](double a, double b) { return std::abs(std::log(a)) < std::abs(std::log(b)); });
  QStringList initialPeriods; for(double p: defaults.mapPeriods) initialPeriods << QString::number(p, 'g', 17);
  auto *periods = new QLineEdit(initialPeriods.isEmpty() ? QString::number(initialPeriod, 'g', 17) : initialPeriods.join(", "), &dialog);
  periods->setObjectName("reportMapPeriods");
  auto *all = new QPushButton(tr("All"), &dialog); all->setObjectName("reportAllPeriods");
  auto *periodRow = new QHBoxLayout; periodRow->addWidget(periods); periodRow->addWidget(all);
  form->addRow(UiHelp::label(&dialog, tr("Periods (s):"), tr("Comma-separated map periods. One page per period and dataset, with selected phase-tensor ellipses and induction vectors together on the same map. Each station uses its nearest period within the map window's tolerance; no interpolation. All selects the survey's exact periods."), "reportPeriodsHelp"), periodRow);
  QObject::connect(all, &QPushButton::clicked, &dialog, [&] { periods->setText(allPeriods.join(", ")); });
  auto *real = checkbox(tr("Real"), "reportRealArrows", defaults.mapSettings.real);
  auto *imaginary = checkbox(tr("Imaginary"), "reportImagArrows", defaults.mapSettings.imaginary);
  auto *parts = new QHBoxLayout; parts->addWidget(real); parts->addWidget(imaginary); parts->addStretch();
  form->addRow(tr("Map vectors:"), parts);
  int mapLayers = defaults.mapLayers;
  form->addRow(tr("Geography:"), MapBackground::button(&dialog, mapLayers, [&](int mask) { mapLayers = mask; }, true));
  layout->addLayout(form);
  layout->addWidget(UiHelp::label(&dialog, tr("Report"), tr("Choose a survey overview, one page per station, and/or maps at selected periods. Station plots retain current phase wrapping and fixed GUI Y ranges. All components are included unless Visible components only is checked. Tipper lines connect usable observations; arrows use the station plot convention. Map pages inherit the period map window's sizes, colors, tolerance, station labels and vector convention. Maps omit masked, disabled, missing and unmatched values. UTM origin coordinates describe the offsets on each map. Full / Part. / Mask / Miss. count complete, partial, masked and missing periods per data type, independent of plot visibility. Export leaves your survey and current plots unchanged."), "reportHelp"));
  auto *pages = new QLabel(&dialog); pages->setObjectName("reportPageCount"); layout->addWidget(pages);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel, &dialog); layout->addWidget(buttons);
  buttons->button(QDialogButtonBox::Save)->setText(tr("Export…"));
  auto updatePages = [&] {
    const bool maps = ptMaps->isChecked() || induction->isChecked();
    for(auto *control: QList<QWidget *>{mapData, periods, all}) control->setEnabled(maps);
    real->setEnabled(induction->isChecked()); imaginary->setEnabled(induction->isChecked());
    tipper->setEnabled(stations->isChecked()); errors->setEnabled(stations->isChecked()); visible->setEnabled(stations->isChecked());
    unsigned count = 0; for(const auto &name: survey.get_stations_names()) if(scope->currentIndex() == 0 || survey.is_active(name)) ++count;
    try {
      if(!count) throw std::invalid_argument("No stations match the selection.");
      if(maps && mapData->currentIndex() != 0 && response->currentData().toString().isEmpty()) throw std::invalid_argument("Select a response for response maps.");
      if(induction->isChecked() && !real->isChecked() && !imaginary->isChecked()) throw std::invalid_argument("Choose real or imaginary vectors.");
      unsigned total = unsigned(overview->isChecked()) + (stations->isChecked() ? count : 0);
      if(maps) total += parseMapPeriods(periods->text()).size() * (mapData->currentIndex() == 2 ? 2 : 1);
      if(!total) throw std::invalid_argument("Select at least one report section.");
      pages->setText(tr("%1 pages").arg(total)); buttons->button(QDialogButtonBox::Save)->setEnabled(true);
    } catch(const std::exception &e) { pages->setText(QString::fromUtf8(e.what())); buttons->button(QDialogButtonBox::Save)->setEnabled(false); }
  };
  for(auto *box: {overview, stations, ptMaps, induction, real, imaginary}) QObject::connect(box, &QCheckBox::toggled, &dialog, updatePages);
  for(auto *combo: {scope, mapData, response}) QObject::connect(combo, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, updatePages);
  QObject::connect(periods, &QLineEdit::textChanged, &dialog, updatePages); updatePages();
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if(dialog.exec() != QDialog::Accepted) return;
  QFileDialog destination(parent, tr("Export survey report"), lastDirectory + "/survey_report.pdf", tr("PDF (*.pdf)"));
  destination.setAcceptMode(QFileDialog::AcceptSave);
  destination.setDefaultSuffix("pdf"); // Resolve the suffix before Qt checks whether the destination exists.
  if(destination.exec() != QDialog::Accepted || destination.selectedFiles().isEmpty()) return;
  const auto path = destination.selectedFiles().first();
  Options options = defaults; options.title = title->text().trimmed(); if(options.title.isEmpty()) options.title = tr("Survey data report");
  options.mapLayers = mapLayers;
  options.overview = overview->isChecked(); options.stationPages = stations->isChecked(); options.tipperArrows = tipper->currentIndex() == 1;
  options.phaseTensorMaps = ptMaps->isChecked(); options.inductionMaps = induction->isChecked();
  options.mapData = Options::MapData(mapData->currentIndex());
  if(options.phaseTensorMaps || options.inductionMaps) options.mapPeriods = parseMapPeriods(periods->text());
  options.mapSettings.real = real->isChecked(); options.mapSettings.imaginary = imaginary->isChecked();
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
