#include "include/MapBackground.h"
#include "include/OsmBasemap.h"
#include <QDesktopServices>
#include <QHBoxLayout>
#include <QLabel>
#include <QSlider>
#include <QWidgetAction>
#include <QDataStream>
#include <QFile>
#include <QDir>
#include <QStandardPaths>
#include <QtEndian>
#include <QMenu>
#include <QToolButton>
#include <proj.h>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace {
const QStringList files{"coast", "countries", "states", "rivers", "lakes"};
const std::vector<QPolygonF> &geometry(int dataset) {
  // External files are immutable during a session; restart after installing data.
  static std::map<QString, std::vector<QPolygonF>> cache;
  static const std::vector<QPolygonF> empty;
  const auto directory = MapBackground::dataDirectory();
  if(directory.isEmpty()) return empty;
  const auto path = QDir(directory).absoluteFilePath(files[dataset] + ".bin");
  const auto found = cache.find(path);
  if(found != cache.end()) return found->second;
  auto &result = cache[path];
  QFile file(path);
  if(!file.open(QIODevice::ReadOnly)) return result;
  auto invalid = [&]() -> const std::vector<QPolygonF> & {
    result.clear(); qWarning() << "Invalid map geometry:" << path; return result;
  };
  if(file.size() < 4 || file.size() > 64 * 1024 * 1024) return invalid();
  const auto packed = file.readAll();
  const auto size = qFromBigEndian<quint32>(reinterpret_cast<const uchar *>(packed.constData()));
  if(size < 4 || size > 128 * 1024 * 1024) return invalid();
  const auto data = qUncompress(packed);
  if(quint32(data.size()) != size) return invalid();
  QDataStream input(data); input.setByteOrder(QDataStream::LittleEndian); input.setFloatingPointPrecision(QDataStream::SinglePrecision);
  quint32 count = 0; input >> count;
  if(count > quint32(data.size() / 4)) return invalid();
  for(quint32 i = 0; i < count; ++i) {
    quint32 points = 0; input >> points;
    if(points < 2 || points > quint64(input.device()->bytesAvailable()) / 8) return invalid();
    QPolygonF line; line.reserve(points);
    for(quint32 j = 0; j < points; ++j) {
      float lon, lat; input >> lon >> lat;
      if(!std::isfinite(lon) || !std::isfinite(lat) || std::abs(lon) > 180. || std::abs(lat) > 90.) return invalid();
      line << QPointF(lon, lat);
    }
    result.push_back(std::move(line));
  }
  if(input.status() != QDataStream::Ok || !input.atEnd()) return invalid();
  return result;
}
}
QString MapBackground::dataDirectory() {
  const auto override = qEnvironmentVariable("EDITOOLS_MAP_DATA_DIR");
  if(!override.isEmpty()) return override;
  return QDir(QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation))
    .filePath("EDITools/natural-earth");
}
int MapBackground::availableLayers() {
  int mask = 0;
  for(int i = 0; i < 3; ++i) if(!geometry(i).empty()) mask |= 1 << i;
  if(!geometry(3).empty() || !geometry(4).empty()) mask |= Water;
  return mask;
}
MapBackground::MapBackground(QCustomPlot *plot): QCPLayerable(plot) {
  setObjectName("mapBackground");
  osm = new Osm::Basemap(plot);
  if(!plot->layer("geography")) plot->addLayer("geography", plot->layer("grid"), QCustomPlot::limBelow);
  setLayer("geography");
}
void MapBackground::setCoordinates(const SurveyCoordinates &value, double units) {
  if(coordinates.utm == value.utm && coordinates.zone == value.zone && coordinates.north == value.north &&
     coordinates.origin_easting == value.origin_easting && coordinates.origin_northing == value.origin_northing && scale == units) return;
  coordinates = value; scale = units;
  osm->setCoordinates(value, units);
  ready.fill(false); for(auto &paths: projected) paths.clear();
}
void MapBackground::setLayers(int layers) {
  const int previous = enabled; enabled = layers ? layers & (availableLayers() | OpenStreetMap) : 0;
  if((previous & OpenStreetMap) != (enabled & OpenStreetMap)) osm->setEnabled(enabled & OpenStreetMap);
}
QToolButton *MapBackground::button(QWidget *parent, int initial, const std::function<void(int)> &changed, bool reportOnly) {
  auto *button = new QToolButton(parent); button->setObjectName("mapLayers"); button->setText(QObject::tr("Map layers…"));
  button->setPopupMode(QToolButton::InstantPopup);
  button->setToolTip(QObject::tr("Geographic backgrounds for maps and exports. Natural Earth uses optional local files; OpenStreetMap loads online tiles for the visible area. See README for setup."));
  auto *menu = new QMenu(button); button->setMenu(menu);
  const QStringList labels{QObject::tr("Coastlines"), QObject::tr("Country boundaries"), QObject::tr("State / province boundaries"), QObject::tr("Rivers and lakes")};
  const int available = availableLayers();
  for(int i = 0; i < labels.size(); ++i) {
    auto *action = menu->addAction(labels[i]); action->setObjectName(QString("mapLayer%1").arg(1 << i));
    action->setCheckable(true); action->setChecked(initial & available & (1 << i)); action->setEnabled(available & (1 << i));
    if(!action->isEnabled()) action->setToolTip(QObject::tr("Optional map data unavailable. See README: Geographic map layers.")); action->setData(1 << i);
    QObject::connect(action, &QAction::toggled, button, [menu, changed] {
      int mask = 0; for(auto *item: menu->actions()) if(item->isChecked()) mask |= item->data().toInt(); changed(mask);
    });
  }
  if(!available) {
    auto *hint = menu->addAction(QObject::tr("Natural Earth: install optional data (README)")); hint->setEnabled(false);
  }
  menu->addSeparator();
  auto *online = menu->addAction(reportOnly ? QObject::tr("OpenStreetMap (viewed tiles)") : QObject::tr("OpenStreetMap (online)"));
  online->setObjectName("mapLayer16"); online->setData(OpenStreetMap); online->setCheckable(true); online->setChecked(initial & OpenStreetMap);
  online->setToolTip(QObject::tr("Loads tiles for the visible map area from OpenStreetMap. Exports reuse available tiles without downloading. Pan/zoom to view the area before exporting. Requires internet for new areas; no offline downloads."));
  QObject::connect(online, &QAction::toggled, button, [menu, changed] {
    int mask = 0; for(auto *item: menu->actions()) if(item->isChecked()) mask |= item->data().toInt(); changed(mask);
  });
  auto *opacityAction = new QWidgetAction(menu); auto *opacityRow = new QWidget(menu); auto *layout = new QHBoxLayout(opacityRow);
  layout->setContentsMargins(12, 4, 12, 4); layout->addWidget(new QLabel(QObject::tr("OSM opacity"), opacityRow));
  auto *opacity = new QSlider(Qt::Horizontal, opacityRow); opacity->setObjectName("osmOpacity"); opacity->setRange(0, 100); opacity->setValue(Osm::TileStore::opacity());
  opacity->setToolTip(QObject::tr("Background opacity for all OSM maps, including exports.")); layout->addWidget(opacity); opacityAction->setDefaultWidget(opacityRow); menu->addAction(opacityAction);
  QObject::connect(opacity, &QSlider::valueChanged, button, [](int value) { Osm::TileStore::setOpacity(value); });
  QObject::connect(menu, &QMenu::aboutToShow, opacity, [opacity] { const QSignalBlocker block(opacity); opacity->setValue(Osm::TileStore::opacity()); });
  auto *license = menu->addAction(QObject::tr("OSM attribution / licence…"));
  QObject::connect(license, &QAction::triggered, button, [] { QDesktopServices::openUrl(QUrl("https://www.openstreetmap.org/copyright")); });
  menu->addSeparator(); auto *credit = menu->addAction(QObject::tr("Natural Earth · public domain")); credit->setEnabled(false);
  return button;
}
QRect MapBackground::clipRect() const { return mParentPlot->axisRect()->rect(); }
void MapBackground::applyDefaultAntialiasingHint(QCPPainter *painter) const { painter->setAntialiasing(true); }
void MapBackground::prepare(int dataset) {
  if(ready[dataset]) return;
  auto &paths = projected[dataset]; paths.clear();
  std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)> context(proj_context_create(), proj_context_destroy);
  const auto definition = "+proj=utm +zone=" + std::to_string(coordinates.zone) + " +ellps=WGS84" + (coordinates.north ? "" : " +south");
  std::unique_ptr<PJ, decltype(&proj_destroy)> projection(coordinates.utm ? proj_create(context.get(), definition.c_str()) : nullptr, proj_destroy);
  if(coordinates.utm && !projection) return;
  const double radians = std::acos(-1.) / 180.;
  for(const auto &line: geometry(dataset)) {
    QPolygonF chunk;
    auto flush = [&] { if(chunk.size() > 1) paths.push_back({chunk, chunk.boundingRect()}); chunk.clear(); };
    double previousLongitude = NAN;
    for(const auto &point: line) {
      const double lon = point.x(), lat = point.y();
      if(std::isfinite(previousLongitude) && std::abs(lon - previousLongitude) > 180.) flush();
      previousLongitude = lon;
      QPointF target = point;
      if(coordinates.utm) {
        const double delta = std::remainder(lon - (coordinates.zone * 6. - 183.), 360.);
        if(lat < -80. || lat > 84. || std::abs(delta) > 30.) { flush(); continue; }
        proj_errno_reset(projection.get());
        const auto p = proj_trans(projection.get(), PJ_FWD, proj_coord(lon * radians, lat * radians, 0., 0.));
        if(proj_errno(projection.get()) || !std::isfinite(p.xy.x) || !std::isfinite(p.xy.y)) { flush(); continue; }
        target = QPointF((p.xy.x - coordinates.origin_easting) * scale, (p.xy.y - coordinates.origin_northing) * scale);
      }
      chunk << target;
      if(chunk.size() == 256) { flush(); chunk << target; }
    }
    flush();
  }
  ready[dataset] = true;
}
void MapBackground::draw(QCPPainter *painter) {
  if(!(enabled & 15)) return;
  auto *x = mParentPlot->xAxis, *y = mParentPlot->yAxis;
  const auto xr = x->range(), yr = y->range();
  const QColor colors[]{QColor("#7998a6"), QColor("#89909c"), QColor("#a0a6af"), QColor("#9bbfcf"), QColor("#8db3c6")};
  const int bits[]{Coast, Countries, States, Water, Water};
  for(int dataset = 0; dataset < 5; ++dataset) if(enabled & bits[dataset]) {
    prepare(dataset); painter->setPen(QPen(colors[dataset], dataset == 1 ? 1. : .8, dataset == 2 ? Qt::DashLine : Qt::SolidLine));
    for(const auto &path: projected[dataset]) {
      const auto &b = path.bounds;
      if(b.right() < xr.lower || b.left() > xr.upper || b.bottom() < yr.lower || b.top() > yr.upper) continue;
      QPolygonF pixels; pixels.reserve(path.points.size());
      for(const auto &p: path.points) pixels << QPointF(x->coordToPixel(p.x()), y->coordToPixel(p.y()));
      painter->drawPolyline(pixels);
    }
  }
  painter->setPen(QColor("#697586")); QFont credit = painter->font(); credit.setPixelSize(10); painter->setFont(credit);
  const auto rect = clipRect().adjusted(5, 5, -5, -3);
  painter->drawText(rect, (enabled & OpenStreetMap ? Qt::AlignLeft : Qt::AlignRight) | Qt::AlignBottom, "Natural Earth");
}
