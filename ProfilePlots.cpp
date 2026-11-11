#include "ProfilePlots.h"
#include <QToolTip>
#include <algorithm>
#include <cmath>

class ProfileSectionPlot::Canvas : public QCPLayerable {
public:
  explicit Canvas(QCustomPlot *plot) : QCPLayerable(plot) { setLayer("main"); }
  std::vector<ProfileGlyph> samples;
  std::vector<ProfileAnalysis::Offset> stations;
  std::vector<std::array<double, 2>> bounds;
  QString selected, up = tr("geographic north");
  bool frequency = false, ellipses = false, cells = true, logarithmic = false;
  double glyphSize = 12.;
  QString label;
  QCPColorGradient colors;
  QCPRange range;
  QPointF center(const ProfileGlyph &s) const {
    return {mParentPlot->xAxis->coordToPixel(s.along),
            mParentPlot->yAxis->coordToPixel(frequency ? 1. / s.period : s.period)};
  }
  QRectF cell(const ProfileGlyph &s) const {
    const double a = std::pow(10., frequency ? -s.logBottom : s.logBottom),
                 b = std::pow(10., frequency ? -s.logTop : s.logTop);
    return QRectF(QPointF(mParentPlot->xAxis->coordToPixel(s.left), mParentPlot->yAxis->coordToPixel(a)),
                  QPointF(mParentPlot->xAxis->coordToPixel(s.right), mParentPlot->yAxis->coordToPixel(b)))
        .normalized();
  }

protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(ellipses || !cells); }
  void draw(QCPPainter *p) override {
    for (const auto &s : samples) {
      if (!std::isfinite(s.value) || (logarithmic && s.value <= 0.))
        continue;
      const auto point = center(s);
      const QColor color = QColor::fromRgb(
          colors.color(std::max(range.lower, std::min(range.upper, s.value)), range, logarithmic));
      if (ellipses) {
        if (!clipRect().adjusted(-glyphSize, -glyphSize, glyphSize, glyphSize).contains(point.toPoint()))
          continue;
        p->save();
        p->translate(point);
        p->rotate(s.tensor.azimuth);
        p->setPen(QPen(QColor("#444444"), .5));
        p->setBrush(color);
        const double radius = .5 * glyphSize;
        p->drawEllipse(
            QRectF(-radius * s.tensor.axisRatio, -radius, 2 * radius * s.tensor.axisRatio, 2 * radius));
        p->restore();
      } else if (cells) {
        const auto rect = cell(s);
        if (rect.intersects(clipRect()))
          p->fillRect(rect, color);
      } else {
        p->setPen(Qt::NoPen);
        p->setBrush(color);
        p->drawEllipse(point, 3.5, 3.5);
      }
    }
    for (unsigned i = 0; i < stations.size(); ++i)
      if (QString::fromStdString(stations[i].name) == selected) {
        const double left = mParentPlot->xAxis->coordToPixel(bounds[i][0]);
        const double right = mParentPlot->xAxis->coordToPixel(bounds[i][1]);
        p->setBrush(Qt::NoBrush);
        p->setPen(QPen(QColor("#d81b60"), 1.5));
        p->drawRect(
            QRectF(QPointF(left, clipRect().top()), QPointF(right, clipRect().bottom())).normalized());
        break;
      }
  }
};

ProfileSectionPlot::ProfileSectionPlot(QWidget *parent) : QCustomPlot(parent), canvas(new Canvas(this)) {
  setMinimumSize(250, 200);
  setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  xAxis->setLabel(tr("Distance along profile (km)"));
  xAxis2->setVisible(true);
  xAxis2->setTickLabelFont(QFont(font().family(), 7));
  xAxis2->setTickLabelRotation(90);
  xAxis2->setSubTicks(false);
  xAxis2->setTickLength(0, 3);
  connect(xAxis, QOverload<const QCPRange &>::of(&QCPAxis::rangeChanged), xAxis2,
          QOverload<const QCPRange &>::of(&QCPAxis::setRange));
  connect(this, &QCustomPlot::mousePress, this, [this](QMouseEvent *event) {
    pressPosition = event->position().toPoint();
    dragged = false;
  });
  connect(this, &QCustomPlot::mouseRelease, this, [this](QMouseEvent *event) {
    if (event->button() != Qt::LeftButton || dragged ||
        (event->position().toPoint() - pressPosition).manhattanLength() > 4 ||
        !axisRect()->rect().contains(event->position().toPoint()))
      return;
    const double along = xAxis->pixelToCoord(event->position().x());
    double best = INFINITY;
    QString name;
    for (unsigned i = 0; i < canvas->stations.size(); ++i) {
      if (along < canvas->bounds[i][0] || along > canvas->bounds[i][1])
        continue;
      const double distance = std::abs(along - canvas->stations[i].along);
      if (distance <= best) {
        best = distance;
        name = QString::fromStdString(canvas->stations[i].name);
      }
    }
    if (!name.isEmpty())
      emit stationClicked(name);
  });
  yAxis->setScaleType(QCPAxis::stLogarithmic);
  yAxis->setTicker(QSharedPointer<QCPAxisTickerLog>::create());
  scale = new QCPColorScale(this);
  scale->setObjectName("profileColorScale");
  plotLayout()->addElement(0, 1, scale);
  scale->setType(QCPAxis::atRight);
  connect(this, &QCustomPlot::mouseMove, this, [this](QMouseEvent *event) {
    if ((event->position().toPoint() - pressPosition).manhattanLength() > 4)
      dragged = true;
    const ProfileGlyph *found = nullptr;
    double best = std::max(10., canvas->glyphSize * .6);
    for (const auto &s : canvas->samples) {
      const auto delta = canvas->center(s) - event->position();
      const double d = std::hypot(delta.x(), delta.y());
      if ((!canvas->ellipses && canvas->cells && canvas->cell(s).contains(event->position())) || d < best) {
        found = &s;
        best = d;
      }
    }
    if (!found || !axisRect()->rect().contains(event->position().toPoint())) {
      setToolTip({});
      return;
    }
    QString text = tr("%1\nDistance: %2 km\nPeriod: %3 s\n%4: %5")
                       .arg(found->station)
                       .arg(found->along, 0, 'g', 7)
                       .arg(found->period, 0, 'g', 8)
                       .arg(canvas->label)
                       .arg(found->value, 0, 'g', 7);
    if (canvas->ellipses)
      text += tr("\nAzimuth in plot frame: %1°; skew: %2°\nEllipse up = %3; right = +90°.")
                  .arg(found->tensor.azimuth, 0, 'f', 2)
                  .arg(found->tensor.skew, 0, 'f', 2)
                  .arg(canvas->up);
    setToolTip(text);
  });
}

void ProfileSectionPlot::setStations(const std::vector<ProfileAnalysis::Offset> &stations, double cellWidth) {
  canvas->stations = stations;
  std::vector<double> positions;
  QMap<double, QString> labels;
  QMap<double, unsigned> counts;
  for (const auto &station : stations) {
    positions.push_back(station.along);
    labels[station.along] = QString::fromStdString(station.name);
    ++counts[station.along];
  }
  const QFontMetrics font(xAxis2->tickLabelFont());
  for (auto it = labels.begin(); it != labels.end(); ++it) {
    const auto extra = counts[it.key()] - 1;
    const auto suffix = extra ? QString(" (+%1)").arg(extra) : QString();
    it.value() = font.elidedText(it.value(), Qt::ElideRight, 140 - font.horizontalAdvance(suffix)) + suffix;
  }
  canvas->bounds = ProfileAnalysis::cells(positions, cellWidth);
  auto ticker = QSharedPointer<QCPAxisTickerText>::create();
  ticker->addTicks(labels);
  xAxis2->setTicker(ticker);
  xAxis2->setRange(xAxis->range());
}
void ProfileSectionPlot::setSelectedStation(const QString &name) {
  canvas->selected = name;
  replot(QCustomPlot::rpQueuedReplot);
}

void ProfileSectionPlot::setReference(const QString &up) { canvas->up = up; }

void ProfileSectionPlot::setSamples(std::vector<ProfileGlyph> samples, bool frequency, bool ellipses,
                                    bool cells, double glyphSize, const QString &label,
                                    QCPColorGradient colors, QCPRange range, bool logarithmic) {
  canvas->samples = std::move(samples);
  canvas->frequency = frequency;
  canvas->ellipses = ellipses;
  canvas->cells = cells;
  canvas->glyphSize = glyphSize;
  canvas->label = label;
  canvas->colors = colors;
  canvas->range = range;
  canvas->logarithmic = logarithmic;
  setProperty("profileSampleCount", unsigned(canvas->samples.size()));
  yAxis->setLabel(frequency ? tr("Frequency (Hz)") : tr("Period (s)"));
  yAxis->setRangeReversed(!frequency);
  scale->setDataScaleType(logarithmic ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
  if (logarithmic)
    scale->axis()->setTicker(QSharedPointer<QCPAxisTickerLog>::create());
  else
    scale->axis()->setTicker(QSharedPointer<QCPAxisTicker>::create());
  scale->setDataRange(range);
  scale->setGradient(colors);
  scale->axis()->setLabel(label);
  replot();
}

class ProfileRosePlot::Canvas : public QCPLayerable {
public:
  explicit Canvas(QCustomPlot *plot) : QCPLayerable(plot) { setLayer("main"); }
  std::vector<Panel> panels;
  QString title;
  bool strikeFrame = false;

protected:
  QRect clipRect() const override { return mParentPlot->viewport(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(true); }
  void draw(QCPPainter *p) override {
    const QRect area = mParentPlot->viewport().adjusted(15, 15, -15, -15);
    p->setPen(Qt::black);
    p->drawText(QRect(area.left(), area.top(), area.width(), 30), Qt::AlignCenter, title);
    if (panels.empty())
      return;
    const double width = double(area.width()) / panels.size();
    const double radius = std::max(15., std::min(width * .38, (area.height() - 170.) * .5));
    double maximum = 0.;
    for (const auto &panel : panels)
      for (double value : panel.histogram.percent)
        maximum = std::max(maximum, value);
    const double radians = std::acos(-1.) / 180.;
    for (unsigned i = 0; i < panels.size(); ++i) {
      const auto &panel = panels[i];
      const auto &hist = panel.histogram;
      const QPointF center(area.left() + width * (i + .5), area.top() + 70. + radius);
      p->setPen(Qt::black);
      p->drawText(QRectF(area.left() + width * i, area.top() + 35., width, 25.), Qt::AlignCenter,
                  panel.title);
      p->setPen(QPen(QColor("#d1d5db"), .7));
      p->setBrush(Qt::NoBrush);
      for (double fraction : {.25, .5, .75, 1.})
        p->drawEllipse(center, radius * fraction, radius * fraction);
      p->drawLine(center - QPointF(radius, 0), center + QPointF(radius, 0));
      p->drawLine(center - QPointF(0, radius), center + QPointF(0, radius));
      p->setPen(QPen(panel.color, .6));
      QColor fill = panel.color;
      fill.setAlpha(140);
      p->setBrush(fill);
      const double bin = hist.percent.empty() ? 180. : 180. / hist.percent.size();
      for (unsigned j = 0; j < hist.percent.size(); ++j)
        for (double opposite : {0., 180.}) {
          const double r = maximum > 0. ? radius * hist.percent[j] / maximum : 0.;
          if (r == 0.)
            continue;
          QPainterPath wedge(center);
          for (int step = 0; step <= 8; ++step) {
            const double angle = (j * bin + bin * step / 8. + opposite) * radians;
            wedge.lineTo(center + QPointF(r * std::sin(angle), -r * std::cos(angle)));
          }
          wedge.closeSubpath();
          p->drawPath(wedge);
        }
      p->setPen(Qt::black);
      p->drawText(QRectF(center.x() - 20, center.y() - radius - 22, 40, 20), Qt::AlignCenter,
                  strikeFrame ? "x" : "N");
      p->drawText(QRectF(center.x() - 20, center.y() + radius + 2, 40, 20), Qt::AlignCenter,
                  strikeFrame ? "−x" : "S");
      p->drawText(QRectF(center.x() + radius + 2, center.y() - 10, 20, 20), Qt::AlignCenter,
                  strikeFrame ? "y" : "E");
      p->drawText(QRectF(center.x() - radius - 22, center.y() - 10, 20, 20), Qt::AlignCenter,
                  strikeFrame ? "−y" : "W");
      const QString mean =
          std::isfinite(hist.mean) ? QString::number(hist.mean, 'f', 1) + "°" : tr("undefined");
      p->drawText(QRectF(area.left() + width * i, center.y() + radius + 30, width, 85),
                  Qt::AlignHCenter | Qt::TextWordWrap,
                  tr("%1 samples · %2 stations\nAxial mean %3 · R = %4\nOuter ring %5% (shared scale)")
                      .arg(hist.samples)
                      .arg(hist.stations)
                      .arg(mean)
                      .arg(hist.concentration, 0, 'f', 2)
                      .arg(maximum, 0, 'f', 1));
    }
  }
};
ProfileRosePlot::ProfileRosePlot(QWidget *parent) : QCustomPlot(parent), canvas(new Canvas(this)) {
  setMinimumSize(400, 240);
  for (auto *axis : {xAxis, yAxis, xAxis2, yAxis2})
    axis->setVisible(false);
}
void ProfileRosePlot::setPanels(std::vector<Panel> panels, const QString &title, bool strikeFrame) {
  unsigned n = 0;
  for (const auto &p : panels)
    n += p.histogram.samples;
  setProperty("profileDirectionCount", n);
  canvas->panels = std::move(panels);
  canvas->title = title;
  canvas->strikeFrame = strikeFrame;
  replot();
}
