#include "include/OsmBasemap.h"
#include <QApplication>
#include <QEvent>
#include <QDesktopServices>
#include <QHelpEvent>
#include <QMouseEvent>
#include <QToolTip>
#include <proj.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace Osm {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double mercatorLimit = 85.0511287798066;
class Credit : public QCPLayerable {
public:
  Credit(QCustomPlot *plot, Basemap *owner): QCPLayerable(plot, "overlay"), owner(owner) { setObjectName("osmAttribution"); }
protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *) const override {}
  void draw(QCPPainter *painter) override { if(owner) owner->drawCredit(painter); }
private:
  QPointer<Basemap> owner;
};
bool finite(const QPointF &p) { return std::isfinite(p.x()) && std::isfinite(p.y()); }
}
QString View::signature() const {
  QString result;
  for(double value: {left, right, top, bottom, unitsPerMetre, coordinates.origin_easting, coordinates.origin_northing}) result += QString::number(value, 'g', 17) + ",";
  return result + QString("%1,%2,%3,%4,%5").arg(coordinates.utm).arg(coordinates.zone).arg(coordinates.north).arg(size.width()).arg(size.height());
}
Frame render(const View &view, TileStore &store, QObject *viewer) {
  Frame result;
  if(view.size.width() < 2 || view.size.height() < 2 || view.unitsPerMetre <= 0. ||
     !std::isfinite(view.left) || !std::isfinite(view.right) || !std::isfinite(view.top) || !std::isfinite(view.bottom)) return result;
  // Limit work at very large export resolutions; tiles retain screen resolution.
  const auto size = view.size.scaled(1600, 1600, Qt::KeepAspectRatio);
  const int width = std::min(view.size.width(), size.width()), height = std::min(view.size.height(), size.height());
  const int columns = (width + 15) / 16, rows = (height + 15) / 16;
  std::vector<QPointF> grid((columns + 1) * (rows + 1));
  std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)> context(proj_context_create(), proj_context_destroy);
  const auto definition = "+proj=utm +zone=" + std::to_string(view.coordinates.zone) + " +ellps=WGS84" + (view.coordinates.north ? "" : " +south");
  std::unique_ptr<PJ, decltype(&proj_destroy)> projection(view.coordinates.utm ? proj_create(context.get(), definition.c_str()) : nullptr, proj_destroy);
  if(view.coordinates.utm && !projection) return result;
  const double reference = view.coordinates.utm ? view.coordinates.zone * 6. - 183. : .5 * (view.left + view.right);
  double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY;
  for(int row = 0; row <= rows; ++row) for(int column = 0; column <= columns; ++column) {
    double lon = view.left + (view.right - view.left) * column / columns;
    double lat = view.top + (view.bottom - view.top) * row / rows;
    if(projection) {
      proj_errno_reset(projection.get());
      const auto point = proj_trans(projection.get(), PJ_INV, proj_coord(lon / view.unitsPerMetre + view.coordinates.origin_easting,
                                     lat / view.unitsPerMetre + view.coordinates.origin_northing, 0., 0.));
      lon = point.lp.lam * 180. / pi; lat = point.lp.phi * 180. / pi;
      if(proj_errno(projection.get()) || lat < -80. || lat > 84.) lon = NAN;
      lon = reference + std::remainder(lon - reference, 360.);
    }
    QPointF point(NAN, NAN);
    if(std::isfinite(lon) && std::isfinite(lat) && std::abs(lat) <= mercatorLimit)
      point = QPointF((lon + 180.) / 360., .5 - std::asinh(std::tan(lat * pi / 180.)) / (2. * pi));
    grid[row * (columns + 1) + column] = point;
    if(finite(point)) { minX = std::min(minX, point.x()); maxX = std::max(maxX, point.x()); minY = std::min(minY, point.y()); maxY = std::max(maxY, point.y()); }
  }
  if(!std::isfinite(minX) || maxX - minX > 2.) return result;
  const double density = std::min(width / std::max(1e-12, maxX - minX), height / std::max(1e-12, maxY - minY));
  int zoom = std::max(0, std::min(19, int(std::floor(std::log2(density / 256.)))));
  // A single bounded zoom level, without adjacent-area or higher-zoom prefetch.
  while(zoom > 0 && (std::ceil((maxX - minX) * (1 << zoom)) + 2) * (std::ceil((maxY - minY) * (1 << zoom)) + 2) > 48) --zoom;
  const int count = 1 << zoom;
  struct Sample { QImage image; int zoom = 0; };
  std::map<Tile, Sample> samples;
  result.image = QImage(width, height, QImage::Format_ARGB32); result.image.fill(Qt::transparent);
  for(int y = 0; y < height; ++y) {
    auto *pixels = reinterpret_cast<QRgb *>(result.image.scanLine(y));
    const double gy = double(y) * rows / (height - 1); const int row = std::min(rows - 1, int(gy)); const double fy = gy - row;
    for(int x = 0; x < width; ++x) {
      const double gx = double(x) * columns / (width - 1); const int column = std::min(columns - 1, int(gx)); const double fx = gx - column;
      const int index = row * (columns + 1) + column;
      const auto a = grid[index], b = grid[index + 1], c = grid[index + columns + 1], d = grid[index + columns + 2];
      if(!finite(a) || !finite(b) || !finite(c) || !finite(d)) continue;
      const auto point = (a * (1. - fx) + b * fx) * (1. - fy) + (c * (1. - fx) + d * fx) * fy;
      const double wx = point.x() - std::floor(point.x()), wy = point.y();
      if(wy < 0. || wy >= 1.) continue;
      const Tile tile{zoom, std::min(count - 1, int(wx * count)), std::min(count - 1, int(wy * count))};
      auto found = samples.find(tile);
      if(found == samples.end()) {
        Sample sample;
        for(int z = zoom; z >= 0; --z) {
          sample.image = store.image({z, tile.x >> (zoom - z), tile.y >> (zoom - z)}, viewer); sample.zoom = z;
          if(!sample.image.isNull()) break;
        }
        found = samples.emplace(tile, sample).first;
      }
      const auto &sample = found->second;
      if(sample.image.isNull()) continue;
      const double tx = wx * (1 << sample.zoom) * 256., ty = wy * (1 << sample.zoom) * 256.;
      pixels[x] = sample.image.pixel(int(tx) % 256, int(ty) % 256);
    }
  }
  for(const auto &sample: samples) { result.tiles.push_back(sample.first); result.loaded += !sample.second.image.isNull(); }
  return result;
}
Basemap::Basemap(QCustomPlot *plot, TileStore *tiles): QCPLayerable(plot), store(tiles ? tiles : TileStore::shared()) {
  setObjectName("osmBasemap");
  if(!plot->layer("osm")) plot->addLayer("osm", plot->layer("grid"), QCustomPlot::limBelow);
  setLayer("osm"); new Credit(plot, this);
  debounce.setSingleShot(true); debounce.setInterval(350);
  connect(&debounce, &QTimer::timeout, this, [this] {
    if(!canRequest() || !haveScreenView) { stopRequests(); return; }
    const auto view = screenView;
    const auto visible = render(view, *store, this);
    requestedSignature = view.signature(); store->replaceView(this, visible.tiles);
    mParentPlot->setProperty("osmRequestedTiles", unsigned(visible.tiles.size()));
  });
  connect(plot, &QCustomPlot::afterReplot, this, [this] {
    // Export rendering never emits afterReplot. Remember only an actual widget
    // frame so temporary PDF sizes/ranges cannot authorize extra downloads.
    if(mParentPlot->isVisible() && mParentPlot->viewport() == mParentPlot->rect()) { screenView = currentView(); haveScreenView = true; }
    scheduleView();
  });
  connect(store, &TileStore::changed, this, [this] {
    frameSignature.clear(); if(enabled) mParentPlot->replot(QCustomPlot::rpQueuedReplot);
  });
  plot->installEventFilter(this); plot->window()->installEventFilter(this);
  connect(qApp, &QGuiApplication::applicationStateChanged, this, [this](Qt::ApplicationState) { scheduleView(); });
}
Basemap::~Basemap() { store->forgetView(this); }
void Basemap::setEnabled(bool on) {
  enabled = on; frameSignature.clear(); requestedSignature.clear();
  creditRect = {}; creditPressed = false;
  if(!on) { stopRequests(); frame = {}; } else scheduleView();
}
void Basemap::setCoordinates(const SurveyCoordinates &value, double unitsPerMetre) {
  coordinates = value; units = unitsPerMetre; haveScreenView = false; frameSignature.clear(); requestedSignature.clear();
  stopRequests(); scheduleView();
}
QRect Basemap::clipRect() const { return mParentPlot->axisRect()->rect(); }
View Basemap::currentView() const {
  const auto rect = clipRect(); View view; view.coordinates = coordinates; view.unitsPerMetre = units; view.size = rect.size();
  view.left = mParentPlot->xAxis->pixelToCoord(rect.left()); view.right = mParentPlot->xAxis->pixelToCoord(rect.right());
  view.top = mParentPlot->yAxis->pixelToCoord(rect.top()); view.bottom = mParentPlot->yAxis->pixelToCoord(rect.bottom()); return view;
}
bool Basemap::canRequest() const {
  return enabled && TileStore::opacity() > 0 && mParentPlot->isVisible() && !mParentPlot->visibleRegion().isEmpty() &&
    !mParentPlot->window()->isMinimized() && mParentPlot->window()->isActiveWindow() &&
    QApplication::applicationState() == Qt::ApplicationActive && mParentPlot->viewport() == mParentPlot->rect();
}
void Basemap::stopRequests() { debounce.stop(); store->forgetView(this); requestedSignature.clear(); }
void Basemap::scheduleView() {
  if(!canRequest()) { stopRequests(); return; }
  if(haveScreenView && screenView.signature() != requestedSignature) {
    // Invalidate the old queue immediately; the debounce delays only new work.
    if(!requestedSignature.isEmpty()) { store->forgetView(this); requestedSignature.clear(); }
    debounce.start();
  }
}
bool Basemap::eventFilter(QObject *object, QEvent *event) {
  if(object == mParentPlot && enabled && TileStore::opacity() > 0) {
    if(event->type() == QEvent::ToolTip) {
      auto *help = static_cast<QHelpEvent *>(event);
      if(creditRect.contains(help->pos())) {
        QToolTip::showText(help->globalPos(), creditTooltip, mParentPlot, creditRect);
        return true;
      }
    }
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      auto *mouse = static_cast<QMouseEvent *>(event);
      if(mouse->button() == Qt::LeftButton) {
        const bool inside = creditRect.contains(mouse->pos());
        if(event->type() == QEvent::MouseButtonPress) {
          creditPressed = inside;
          if(inside) return true;
        } else if(creditPressed) {
          creditPressed = false;
          if(inside) QDesktopServices::openUrl(QUrl("https://www.openstreetmap.org/copyright"));
          return true;
        }
      }
    }
  }
  switch(event->type()) {
    case QEvent::Hide: case QEvent::WindowDeactivate: stopRequests(); break;
    case QEvent::Show: case QEvent::WindowActivate: case QEvent::Resize: case QEvent::Wheel: case QEvent::MouseButtonRelease:
      QTimer::singleShot(0, this, [this] { scheduleView(); }); break;
    default: break;
  }
  return false;
}
void Basemap::draw(QCPPainter *painter) {
  if(!enabled || TileStore::opacity() == 0) return;
  const auto view = currentView(); const auto signature = view.signature();
  if(signature != frameSignature) { frame = render(view, *store, this); frameSignature = signature; }
  painter->save(); painter->setOpacity(TileStore::opacity() / 100.); painter->setRenderHint(QPainter::SmoothPixmapTransform);
  painter->drawImage(clipRect(), frame.image); painter->restore();
  mParentPlot->setProperty("osmLoadedTiles", frame.loaded);
  mParentPlot->setProperty("osmVisibleTiles", unsigned(frame.tiles.size()));
}
void Basemap::drawCredit(QCPPainter *painter) {
  if(!enabled || TileStore::opacity() == 0) return;
  const bool exporting = painter->modes().testFlag(QCPPainter::pmNoCaching);
  QString text = QString::fromUtf8("© OpenStreetMap contributors");
  if(exporting) {
    text += "\nhttps://www.openstreetmap.org/copyright";
    if(frame.loaded < frame.tiles.size() || frame.tiles.empty()) text += "\n" + tr("Partial basemap");
  }
  if(!store->config().extraAttribution.isEmpty()) text += "\n" + store->config().extraAttribution;
  mParentPlot->setProperty("osmAttribution", text);
  painter->save(); QFont font = mParentPlot->font(); font.setPixelSize(11); font.setUnderline(!exporting); painter->setFont(font);
  const auto bounds = clipRect().adjusted(4, 4, -4, -4);
  if(bounds.width() <= 6) { painter->restore(); return; }
  const auto box = painter->fontMetrics().boundingRect(QRect(0, 0, bounds.width() - 6, 300), Qt::TextWordWrap, text);
  const int width = std::min(bounds.width(), box.width() + 6);
  const QRect panel(bounds.right() - width + 1, bounds.bottom() - box.height() - 5, width, box.height() + 6);
  painter->fillRect(panel, Qt::white); painter->setPen(QColor(exporting ? "#26323c" : "#205782"));
  painter->drawText(panel.adjusted(3, 3, -3, -3), Qt::TextWordWrap | Qt::AlignLeft, text); painter->restore();
  if(!exporting) {
    creditRect = panel;
    creditTooltip = tr("OpenStreetMap attribution and licence — click to open.\n%1/%2 tiles available.").arg(frame.loaded).arg(frame.tiles.size());
    if(!store->status().isEmpty()) creditTooltip += "\n" + store->status();
    if(frame.loaded < frame.tiles.size() || frame.tiles.empty()) creditTooltip += "\n" + tr("View this area online to load missing tiles before exporting.");
  }
}
}
