#include "include/OsmTiles.h"
#include "include/OsmBasemap.h"
#include <QApplication>
#include <QBuffer>
#include <QElapsedTimer>
#include <QDir>
#include <QFile>
#include <QPdfWriter>
#include <QSettings>
#include <QTemporaryDir>
#include <QThread>
#include <QTimer>
#include <deque>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition, const char *message) { if(!condition) throw std::runtime_error(message); }
void events(int milliseconds = 30) {
  QElapsedTimer timer; timer.start();
  do { QApplication::processEvents(); QThread::msleep(1); } while(timer.elapsed() < milliseconds);
}
QByteArray png() {
  QImage image(256, 256, QImage::Format_ARGB32);
  for(int y = 0; y < 256; ++y) for(int x = 0; x < 256; ++x) image.setPixel(x, y, qRgb(x, y, 70));
  QByteArray data; QBuffer buffer(&data); buffer.open(QIODevice::WriteOnly); image.save(&buffer, "PNG"); return data;
}
struct Response {
  int status = 200, delay = 0;
  QNetworkReply::NetworkError error = QNetworkReply::NoError;
  QByteArray body = png();
  std::map<QByteArray, QByteArray> headers{{"Cache-Control", "max-age=604800"}, {"ETag", "\"tile-v1\""}, {"Last-Modified", "Fri, 25 Sep 2026 12:00:00 GMT"}};
};
class Reply : public QNetworkReply {
public:
  Reply(const QNetworkRequest &request, const Response &response, QObject *parent): QNetworkReply(parent), body(response.body) {
    setRequest(request); setUrl(request.url()); setOperation(QNetworkAccessManager::GetOperation);
    setAttribute(QNetworkRequest::HttpStatusCodeAttribute, response.status);
    for(const auto &header: response.headers) setRawHeader(header.first, header.second);
    if(response.error != QNetworkReply::NoError) setError(response.error, "Simulated transport error");
    open(QIODevice::ReadOnly);
    QTimer::singleShot(response.delay, this, [this] { complete(); });
  }
  void abort() override { if(isFinished()) return; canceled = true; setError(OperationCanceledError, "Canceled"); complete(); }
  bool canceled = false;
  qint64 bytesAvailable() const override { return body.size() - offset + QNetworkReply::bytesAvailable(); }
protected:
  qint64 readData(char *data, qint64 length) override {
    const auto size = std::min<qint64>(length, body.size() - offset);
    if(!size) return -1;
    memcpy(data, body.constData() + offset, size); offset += size; return size;
  }
private:
  void complete() { if(isFinished()) return; setFinished(true); if(!canceled) emit readyRead(); emit finished(); }
  QByteArray body;
  qint64 offset = 0;
};
class Network : public QNetworkAccessManager {
public:
  std::vector<QNetworkRequest> requests;
  std::vector<QPointer<Reply>> replies;
  std::deque<Response> responses;
  Response fallback;
protected:
  QNetworkReply *createRequest(Operation operation, const QNetworkRequest &request, QIODevice *) override {
    check(operation == GetOperation, "Unexpected network operation"); requests.push_back(request);
    const auto response = responses.empty() ? fallback : responses.front(); if(!responses.empty()) responses.pop_front();
    auto *reply = new Reply(request, response, this); replies.push_back(reply); return reply;
  }
};
Osm::Config configuration(const QString &directory, const QString &name) {
  Osm::Config config; config.cacheDirectory = directory + "/" + name; return config;
}
void policy() {
  const auto now = QDateTime::fromString("2026-09-25T12:00:00Z", Qt::ISODate);
  check(Osm::cachePolicy({}, {}, {}, {}, now).expires == now.addDays(7), "Missing headers need seven-day caching");
  check(Osm::cachePolicy("public, max-age=3600", {}, "Fri, 25 Sep 2026 11:50:00 GMT", "900", now).expires == now.addSecs(2700), "Date/Age not subtracted");
  check(Osm::cachePolicy({}, "Fri, 25 Sep 2026 13:00:00 GMT", "Fri, 25 Sep 2026 11:00:00 GMT", {}, now).expires == now.addSecs(3600), "Expires ignored");
  check(Osm::cachePolicy("max-age=3600, no-cache", {}, {}, {}, now).expires == now, "no-cache must revalidate");
  check(Osm::cachePolicy({}, "0", {}, {}, now).expires == now, "Invalid Expires must count as expired");
  check(Osm::cachePolicy("no-store", {}, {}, {}, now).noStore, "no-store ignored");
  check(Osm::cachePolicy("max-age=\"120\"", {}, {}, {}, now).expires == now.addSecs(120), "Quoted max-age ignored");
}
void cacheAndHeaders(const QString &directory) {
  Network network; auto config = configuration(directory, "cache");
  auto now = QDateTime::currentDateTimeUtc(); config.clock = [&] { return now; };
  Osm::TileStore store(config, nullptr, &network); QObject viewer, secondViewer;
  const Osm::Tile tile{3, 2, 3};
  check(store.image(tile).isNull() && network.requests.empty(), "Cache lookup triggered network");
  store.replaceView(&viewer, {tile}); store.replaceView(&secondViewer, {tile}); events();
  check(network.requests.size() == 1, "Duplicate visible windows downloaded the same tile");
  const auto &request = network.requests.front();
  check(request.url() == QUrl("https://tile.openstreetmap.org/3/2/3.png"), "Wrong public tile URL");
  check(request.rawHeader("User-Agent").startsWith("EDITools") && request.rawHeader("User-Agent").contains("github.com/GoFEM/EDITools"), "Missing application identification");
  check(!request.hasRawHeader("Cache-Control") && !request.hasRawHeader("Pragma"), "Request bypasses HTTP caches");
  check(!request.hasRawHeader("Referer") && !request.hasRawHeader("Authorization"), "Unexpected native-app headers");
  check(request.attribute(QNetworkRequest::RedirectPolicyAttribute).toInt() == QNetworkRequest::SameOriginRedirectPolicy, "Redirects may leave provider origin");
  check(!store.image(tile).isNull() && QFile::exists(store.cachePath(tile)), "Tile was not cached");
  store.replaceView(&viewer, {tile}); events(); check(network.requests.size() == 1, "Fresh tile requested again");
  Network otherNetwork; Osm::TileStore reopened(config, nullptr, &otherNetwork);
  reopened.replaceView(&viewer, {tile}); events(); check(otherNetwork.requests.empty() && !reopened.image(tile).isNull(), "Persistent cache was not reused");
  now = now.addDays(8);
  check(reopened.image(tile).isNull(), "Expired cache entry reused without revalidation");
  Response unchanged; unchanged.status = 304; unchanged.body.clear(); unchanged.headers = {{"Cache-Control", "max-age=3600"}, {"ETag", "\"tile-v1\""}};
  otherNetwork.responses.push_back(unchanged);
  reopened.replaceView(&viewer, {tile}); events();
  check(otherNetwork.requests.size() == 1, "Expired tile not revalidated");
  check(otherNetwork.requests.back().rawHeader("If-None-Match") == "\"tile-v1\"" &&
        otherNetwork.requests.back().rawHeader("If-Modified-Since").startsWith("Fri,"), "Conditional cache validators missing");
  check(!reopened.image(tile).isNull(), "304 discarded cached image");
  reopened.replaceView(&viewer, {tile}); events(); check(otherNetwork.requests.size() == 1, "304 did not refresh cache lifetime");
  Network fresh; Osm::TileStore third(config, nullptr, &fresh); third.replaceView(&viewer, {tile}); events();
  check(fresh.requests.empty(), "304 metadata was not persisted");

  Network differentProvider; auto alternative = config; alternative.url = "https://tiles.example.invalid/{z}/{x}/{y}.png";
  Osm::TileStore separateCache(alternative, nullptr, &differentProvider);
  check(separateCache.image(tile).isNull(), "Provider cache entries were mixed");

  Network noStoreNetwork; noStoreNetwork.fallback.headers = {{"Cache-Control", "no-store"}};
  Osm::TileStore noStore(configuration(directory, "no-store"), nullptr, &noStoreNetwork);
  noStore.replaceView(&viewer, {tile}); events();
  check(!noStore.image(tile, &viewer).isNull() && !QFile::exists(noStore.cachePath(tile)), "no-store response persisted");
  noStore.forgetView(&viewer); check(noStore.image(tile).isNull(), "no-store response retained after view closed");
  noStore.replaceView(&viewer, {tile}); events(); check(noStoreNetwork.requests.size() == 2, "no-store revisit skipped network");
}
void limiting(const QString &directory) {
  Network network; network.fallback.delay = 1000;
  Osm::TileStore store(configuration(directory, "limited"), nullptr, &network); QObject viewer;
  store.replaceView(&viewer, {{3, 0, 0}, {3, 1, 0}, {3, 2, 0}, {3, 3, 0}});
  check(network.requests.size() == 2, "Too many concurrent requests");
  const auto first = network.replies[0], second = network.replies[1];
  store.replaceView(&viewer, {{3, 7, 7}});
  check(first->canceled && second->canceled, "Pan did not cancel offscreen requests");
  events(); check(network.requests.size() == 3, "Old viewport queue survived pan");
  const auto last = network.replies.back(); store.forgetView(&viewer); check(last->canceled, "Hidden map continued downloading");
  events(); check(network.requests.size() == 3, "Hidden map drained pending queue");

  for(int status: {403, 429, 503}) {
    Network denied; denied.fallback.status = status; denied.fallback.headers = {{"Retry-After", "120"}};
    auto config = configuration(directory, QString("denied-%1").arg(status));
    auto now = QDateTime::currentDateTimeUtc(); config.clock = [&] { return now; };
    Osm::TileStore blocked(config, nullptr, &denied);
    blocked.replaceView(&viewer, {{3, 0, 0}, {3, 1, 0}, {3, 2, 0}}); events();
    check(denied.requests.size() == 2, "Rejection did not stop the queue");
    blocked.replaceView(&viewer, {{3, 3, 0}}); events(); check(denied.requests.size() == 2, "Rejection/cooldown was bypassed");
    now = now.addSecs(121); blocked.replaceView(&viewer, {{3, 4, 0}}); events();
    check(denied.requests.size() == (status == 403 ? 2 : 3), "Retry-After or access-denied handling incorrect");
  }
  Network tls; tls.fallback.status = 0; tls.fallback.error = QNetworkReply::SslHandshakeFailedError;
  Osm::TileStore tlsFailure(configuration(directory, "tls"), nullptr, &tls);
  tlsFailure.replaceView(&viewer, {{3, 0, 0}, {3, 1, 0}, {3, 2, 0}}); events();
  check(tls.requests.size() == 2 && !tlsFailure.status().isEmpty(), "TLS failure did not stop queued requests");
  for(const auto &request: tls.requests) check(request.url().scheme() == "https", "TLS failure fell back to HTTP");

  Network invalid; auto config = configuration(directory, "http"); config.url = "http://tile.openstreetmap.org/{z}/{x}/{y}.png";
  Osm::TileStore badUrl(config, nullptr, &invalid); badUrl.replaceView(&viewer, {{0, 0, 0}}); events();
  check(invalid.requests.empty() && !badUrl.status().isEmpty(), "Public tiles allowed over HTTP");
}
void projectionsAndExports(const QString &directory) {
  Network network; Osm::TileStore store(configuration(directory, "render"), nullptr, &network); QObject viewer;
  Osm::View view; view.left = -113.1; view.right = -112.9; view.bottom = 37.9; view.top = 38.1; view.size = {400, 400};
  const auto empty = Osm::render(view, store);
  check(network.requests.empty() && !empty.tiles.empty() && empty.tiles.size() <= 48, "Renderer requested or over-selected tiles");
  store.replaceView(&viewer, empty.tiles); events(100);
  const auto loaded = Osm::render(view, store);
  check(loaded.loaded == loaded.tiles.size() && loaded.loaded > 0, "Raster tiles not rendered");
  check(qAlpha(loaded.image.pixel(200, 200)) == 255, "OSM background missing");
  const auto lon = .5 * (view.left + view.right), lat = .5 * (view.top + view.bottom);
  const double pi = std::acos(-1.);
  // Independent Web Mercator image coordinates for the center pixel.
  const auto expected = [&](double longitude, double latitude, int zoom) {
    const double tx = (longitude + 180.) / 360. * (1 << zoom) * 256.;
    const double ty = (.5 - std::asinh(std::tan(latitude * pi / 180.)) / (2. * pi)) * (1 << zoom) * 256.;
    return QPoint(int(tx) % 256, int(ty) % 256);
  };
  const auto sample = expected(view.left + (view.right-view.left) * 200 / 399., view.top + (view.bottom-view.top) * 200 / 399., loaded.tiles.front().z);
  check(std::abs(qRed(loaded.image.pixel(200, 200)) - sample.x()) <= 1 && std::abs(qGreen(loaded.image.pixel(200, 200)) - sample.y()) <= 1, "Geographic pixels misregistered");
  for(bool centered: {false, true}) for(bool north: {false, true}) {
    auto projected = view;
    projected.coordinates = SurveyCoordinates::calculate({{{lat, lon, 0.}}}, 12, north, centered);
    projected.unitsPerMetre = centered ? .001 : 1.;
    const auto center = projected.coordinates.transform({{{lat, lon, 0.}}}).front();
    projected.left = (center[1] - 10000.) * projected.unitsPerMetre; projected.right = (center[1] + 10000.) * projected.unitsPerMetre;
    projected.top = (center[0] + 10000.) * projected.unitsPerMetre; projected.bottom = (center[0] - 10000.) * projected.unitsPerMetre;
    const auto frame = Osm::render(projected, store); store.replaceView(&viewer, frame.tiles); events(100);
    const auto map = Osm::render(projected, store);
    check(map.loaded == map.tiles.size() && map.loaded > 0 && qAlpha(map.image.pixel(200, 200)) == 255, "UTM/origin/units lost basemap");
    const auto centerSample = expected(lon, lat, map.tiles.front().z);
    const auto color = map.image.pixel(200, 200);
    check(std::abs(qRed(color) - centerSample.x()) <= 3 && std::abs(qGreen(color) - centerSample.y()) <= 3,
          "UTM basemap is displaced from the known station location");
  }
  Osm::View dateline = view; dateline.left = 179.99; dateline.right = 180.01; dateline.top = .01; dateline.bottom = -.01;
  const auto wrapped = Osm::render(dateline, store);
  check(!wrapped.tiles.empty() && wrapped.tiles.size() <= 48, "Date line generated unbounded requests");
  for(const auto &tile: wrapped.tiles) check(tile.x == 0 || tile.x == (1 << tile.z) - 1, "Date line tiles unwrapped incorrectly");

  QCustomPlot plot; plot.resize(600, 450); plot.xAxis->setRange(view.left, view.right); plot.yAxis->setRange(view.bottom, view.top);
  auto *background = new Osm::Basemap(&plot, &store); background->setEnabled(true);
  const auto requestsBefore = network.requests.size();
  check(plot.savePdf(directory + "/osm-export.pdf", 1200, 900), "PDF export failed");
  check(plot.toPixmap(1000, 800).save(directory + "/osm-export.png"), "Raster export failed");
  check(plot.property("osmAttribution").toString().contains("© OpenStreetMap contributors") &&
        plot.property("osmAttribution").toString().contains("https://www.openstreetmap.org/copyright"), "Export attribution missing");
  events(400);
  check(network.requests.size() == requestsBefore, "Hidden/export plots downloaded tiles");
  background->setEnabled(false); plot.replot(); events(400); check(network.requests.size() == requestsBefore, "Disabled layer requested tiles");

  // All following network calls use the injected fake transport, never OSM.
  background->setEnabled(true); plot.show(); plot.activateWindow(); plot.raise(); plot.replot(); events(700);
  check(network.requests.size() > requestsBefore || plot.property("osmRequestedTiles").toUInt() > 0, "Visible map did not authorize its viewport");
  check(plot.property("osmAttribution").toString() == QString::fromUtf8("© OpenStreetMap contributors"), "Interactive attribution should be one compact line");
  check(plot.grab().save(directory + "/osm-screen.png"), "Screen preview failed");
  const auto visibleCount = network.requests.size();
  plot.savePdf(directory + "/osm-visible-export.pdf", 2400, 1600); events(500);
  check(network.requests.size() == visibleCount, "Visible PDF export requested higher-resolution tiles");
  plot.hide(); events(100);
  plot.xAxis->setRange(-90., -89.); plot.replot(); events(500);
  check(network.requests.size() == visibleCount, "Hidden window fetched a new area");
}
}
int main(int argc, char **argv) {
  QApplication app(argc, argv); QCoreApplication::setOrganizationName("EDIToolsTests"); QCoreApplication::setApplicationName("OsmTests");
  QTemporaryDir directory; QSettings::setDefaultFormat(QSettings::IniFormat); QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
  try {
    check(directory.isValid(), "Temporary directory unavailable"); policy(); cacheAndHeaders(directory.path()); limiting(directory.path()); projectionsAndExports(directory.path());
    if(argc > 1) {
      const QString destination = QString::fromLocal8Bit(argv[1]); check(QDir().mkpath(destination), "Cannot create artifact directory");
      for(const auto *name: {"osm-export.png", "osm-export.pdf", "osm-visible-export.pdf", "osm-screen.png"})
        check(QFile::copy(directory.path() + "/" + name, destination + "/" + name), "Cannot retain visual-check artifact");
    }
    std::cout << "OSM cache, HTTP policy, projection, visibility and export checks passed (mock transport only).\n";
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
