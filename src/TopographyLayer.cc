#include "include/TopographyLayer.h"
#include "include/PlotColorMaps.h"
#include <cmath>
#include <QPagedPaintDevice>

namespace {
QCPColorGradient palette(const std::string &name) {
  if (name != "Terrain")
    return PlotColorMaps::gradient(QString::fromStdString(name));
  QCPColorGradient g;
  g.setColorStopAt(0., QColor("#47795a"));
  g.setColorStopAt(.35, QColor("#b1bc80"));
  g.setColorStopAt(.65, QColor("#c6aa7c"));
  g.setColorStopAt(.85, QColor("#938275"));
  g.setColorStopAt(1., QColor("#f5f3ec"));
  return g;
}
} // namespace
TopographyLayer::TopographyLayer(QCustomPlot *plot) : QCPLayerable(plot) {
  setObjectName("topographyLayer");
  if (!plot->layer("terrain"))
    plot->addLayer("terrain", plot->layer("grid"), QCustomPlot::limBelow);
  setLayer("terrain");
}
void TopographyLayer::draw(QCPPainter *p) {
  const auto &s = terrain.style;
  if (!terrain.grid || !s.visible || scale <= 0. || !std::isfinite(s.lower) || !std::isfinite(s.upper) ||
      s.upper <= s.lower)
    return;
  const auto &g = *terrain.grid;
  const auto rect = clipRect();
  if (rect.width() < 2 || rect.height() < 2)
    return;
  const auto x = mParentPlot->xAxis->range(), y = mParentPlot->yAxis->range();
  QString key = QString::fromStdString(s.palette);
  for (double v : {x.lower, x.upper, y.lower, y.upper, scale, coordinates.origin_easting,
                   coordinates.origin_northing, s.lower, s.upper, s.interval, double(s.colors),
                   double(s.shade), double(s.contours), double(coordinates.utm),
                   double(Topography::utmEpsg(coordinates)), double(rect.width()), double(rect.height())})
    key += '|' + QString::number(v, 'g', 17);
  if (key != signature || cachedGrid.get() != terrain.grid.get()) {
    const int w = std::min(1000, rect.width()), h = std::min(1000, rect.height());
    image = QImage(w, h, QImage::Format_ARGB32);
    image.fill(Qt::transparent);
    std::vector<double> e(size_t(w) * h), n(e.size()), z(e.size(), NAN);
    for (int r = 0; r < h; ++r)
      for (int c = 0; c < w; ++c) {
        const auto i = size_t(r) * w + c;
        e[i] = x.lower + (c + .5) / w * x.size();
        n[i] = y.upper - (r + .5) / h * y.size();
        if (coordinates.utm) {
          e[i] = e[i] / scale + coordinates.origin_easting;
          n[i] = n[i] / scale + coordinates.origin_northing;
        }
      }
    try {
      Topography::transform(coordinates.utm ? Topography::utmEpsg(coordinates) : 4326, g.epsg, e, n);
    } catch (const std::exception &) {
      return;
    }
    auto gradient = palette(s.palette);
    const QCPRange range(s.lower, s.upper);
    for (size_t i = 0; i < z.size(); ++i)
      z[i] = g.sample(e[i], n[i]);
    for (int r = 0; r < h; ++r)
      for (int c = 0; c < w; ++c) {
        const size_t i = size_t(r) * w + c;
        if (!std::isfinite(z[i]))
          continue;
        QColor color = s.colors ? QColor(gradient.color(z[i], range)) : QColor("#eeeeee");
        double light = 1.;
        if (s.shade) {
          const double l = g.sample(e[i] - g.spacing, n[i]), rr = g.sample(e[i] + g.spacing, n[i]),
                       b = g.sample(e[i], n[i] - g.spacing), t = g.sample(e[i], n[i] + g.spacing);
          if (std::isfinite(l) && std::isfinite(rr) && std::isfinite(b) && std::isfinite(t)) {
            const double dx = (rr - l) / (2 * g.spacing), dy = (t - b) / (2 * g.spacing);
            light =
                .45 + .65 * std::max(0., (.5 * dx - .5 * dy + .70710678) / std::sqrt(1 + dx * dx + dy * dy));
          }
        }
        if (s.contours && s.interval > 0. && std::isfinite(s.interval)) {
          const double level = std::floor(z[i] / s.interval);
          if ((c && std::isfinite(z[i - 1]) && std::floor(z[i - 1] / s.interval) != level) ||
              (r && std::isfinite(z[i - w]) && std::floor(z[i - w] / s.interval) != level))
            light *= .5;
        }
        image.setPixelColor(c, r,
                            QColor(std::min(255, int(color.red() * light)),
                                   std::min(255, int(color.green() * light)),
                                   std::min(255, int(color.blue() * light))));
      }
    signature = key;
    cachedGrid = terrain.grid;
  }
  p->save();
  p->setOpacity(std::max(0., std::min(1., s.opacity)));
  p->drawImage(rect, image);
  p->restore();
  p->save();
  QFont f = p->font();
  f.setPixelSize(10);
  p->setFont(f);
  p->setPen(QColor("#27313b"));
  const QString label = QString("Elevation %1–%2 m")
                            .arg(s.lower, 0, 'f', 0)
                            .arg(s.upper, 0, 'f', 0);
  const QRect area = rect.adjusted(5, 5, -5, -5);
  const auto box = p->boundingRect(area, Qt::AlignRight | Qt::AlignTop | Qt::TextWordWrap, label);
  p->fillRect(box.adjusted(-3, -2, 3, 2), QColor(255, 255, 255, 210));
  p->drawText(area, Qt::AlignRight | Qt::AlignTop | Qt::TextWordWrap, label);
  p->restore();
  if (s.colors) {
    QImage bar(120, 7, QImage::Format_RGB32);
    auto gradient = palette(s.palette);
    for (int c = 0; c < bar.width(); ++c)
      for (int r = 0; r < bar.height(); ++r)
        bar.setPixel(c, r, gradient.color(c, QCPRange(0, bar.width() - 1)));
    p->drawImage(QPoint(rect.right() - bar.width() - 5, box.bottom() + 3), bar);
  }
  QString attribution = QString::fromStdString(g.source);
  if (!g.notice.empty() && dynamic_cast<QPagedPaintDevice *>(p->device()))
    attribution += '\n' + QString::fromStdString(g.notice);
  if (!attribution.isEmpty()) {
    p->save();
    QFont credit = p->font();
    credit.setPixelSize(9);
    p->setFont(credit);
    p->setPen(Qt::black);
    // Leave the bottom row available for the other geographic layers' credits.
    const auto area = rect.adjusted(5, 5, -5, -20);
    const auto box = p->boundingRect(area, Qt::AlignLeft | Qt::AlignBottom | Qt::TextWordWrap, attribution);
    p->fillRect(box.adjusted(-2, -2, 2, 2), QColor(255, 255, 255, 230));
    p->drawText(area, Qt::AlignLeft | Qt::AlignBottom | Qt::TextWordWrap, attribution);
    p->restore();
  }
}
