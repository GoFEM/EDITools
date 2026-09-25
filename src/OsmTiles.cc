#include "include/OsmTiles.h"
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDataStream>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QImageReader>
#include <QBuffer>
#include <QSaveFile>
#include <QSettings>
#include <QStandardPaths>
#include <QTemporaryFile>
#include <QTimer>
#include <algorithm>
#include <tuple>

namespace Osm {
namespace {
constexpr qint64 week = 7 * 24 * 60 * 60;
constexpr quint32 magic = 0x4544494f;
QDateTime httpDate(const QByteArray &value) { return QDateTime::fromString(QString::fromLatin1(value), Qt::RFC2822Date).toUTC(); }
bool validTile(const Tile &t) { return t.z >= 0 && t.z <= 19 && t.x >= 0 && t.y >= 0 && t.x < (1 << t.z) && t.y < (1 << t.z); }
QImage decode(const QByteArray &png) {
  if(png.size() > 1024 * 1024) return {};
  QBuffer buffer; buffer.setData(png); buffer.open(QIODevice::ReadOnly);
  QImageReader reader(&buffer, "PNG");
  if(reader.size() != QSize(256, 256)) return {};
  return reader.read().convertToFormat(QImage::Format_ARGB32);
}
}
QString Tile::key() const { return QString("%1/%2/%3").arg(z).arg(x).arg(y); }
bool Tile::operator<(const Tile &o) const { return std::tie(z, x, y) < std::tie(o.z, o.x, o.y); }
CachePolicy cachePolicy(const QByteArray &control, const QByteArray &expires, const QByteArray &date,
                        const QByteArray &age, const QDateTime &now) {
  CachePolicy result; bool revalidate = false; qint64 lifetime = -1;
  for(auto part: control.toLower().split(',')) {
    part = part.trimmed();
    const auto directive = part.split('=').first().trimmed();
    if(directive == "no-store") result.noStore = true;
    if(directive == "no-cache") revalidate = true;
    if(directive == "max-age") {
      bool ok; const auto seconds = part.mid(part.indexOf('=') + 1).trimmed().replace("\"", "").toLongLong(&ok);
      if(ok && seconds >= 0) lifetime = std::min<qint64>(seconds, 10LL * 365 * 24 * 3600);
    }
  }
  const auto serverDate = httpDate(date), expiryDate = httpDate(expires);
  if(lifetime < 0 && !expires.isEmpty()) lifetime = expiryDate.isValid() ? std::max<qint64>(0, (serverDate.isValid() ? serverDate : now).secsTo(expiryDate)) : 0;
  if(lifetime < 0) { result.expires = now.addSecs(revalidate || result.noStore ? 0 : week); return result; }
  bool ageOk; const auto headerAge = age.toLongLong(&ageOk);
  const auto currentAge = std::max<qint64>(ageOk ? std::max<qint64>(0, headerAge) : 0,
                                         serverDate.isValid() ? std::max<qint64>(0, serverDate.secsTo(now)) : 0);
  result.expires = now.addSecs(revalidate || result.noStore ? 0 : std::max<qint64>(0, lifetime - currentAge));
  return result;
}
TileStore::TileStore(const Config &config, QObject *parent, QNetworkAccessManager *transport)
  : QObject(parent), settings(config), network(transport ? transport : new QNetworkAccessManager(this)) {
  if(settings.cacheDirectory.isEmpty()) settings.cacheDirectory = QDir(QStandardPaths::writableLocation(QStandardPaths::GenericCacheLocation)).filePath("EDITools/osm");
  const auto endpoint = url({0, 0, 0});
  const bool local = endpoint.host() == "127.0.0.1" || endpoint.host() == "localhost" || endpoint.host() == "::1";
  valid = endpoint.isValid() && !endpoint.host().isEmpty() && endpoint.userInfo().isEmpty() && endpoint.fragment().isEmpty() &&
    (endpoint.scheme() == "https" || (local && endpoint.scheme() == "http")) &&
    settings.url.contains("{z}") && settings.url.contains("{x}") && settings.url.contains("{y}");
  if(!valid) message = tr("OSM: invalid tile URL");
}
TileStore *TileStore::shared() {
  static QPointer<TileStore> store;
  if(!store) {
    Config config;
    const auto endpoint = qEnvironmentVariable("EDITOOLS_OSM_TILE_URL"); if(!endpoint.isEmpty()) config.url = endpoint;
    config.extraAttribution = qEnvironmentVariable("EDITOOLS_OSM_ATTRIBUTION");
    store = new TileStore(config, QCoreApplication::instance());
  }
  return store;
}
int TileStore::opacity() { return std::max(0, std::min(100, QSettings().value("maps/osmOpacity", 65).toInt())); }
void TileStore::setOpacity(int percent) { QSettings().setValue("maps/osmOpacity", std::max(0, std::min(100, percent))); emit shared()->changed(); }
QUrl TileStore::url(const Tile &tile) const {
  auto text = settings.url; text.replace("{z}", QString::number(tile.z)).replace("{x}", QString::number(tile.x)).replace("{y}", QString::number(tile.y));
  return QUrl(text);
}
QString TileStore::cachePath(const Tile &tile) const {
  return QDir(settings.cacheDirectory).filePath(QString::fromLatin1(QCryptographicHash::hash(url(tile).toEncoded(), QCryptographicHash::Sha256).toHex()) + ".tile");
}
std::shared_ptr<TileStore::Entry> TileStore::entry(const Tile &tile) {
  const auto key = tile.key(); auto found = entries.find(key);
  if(found != entries.end()) { found->second->touched = ++tick; return found->second; }
  auto value = std::make_shared<Entry>(); value->touched = ++tick;
  QFile file(cachePath(tile));
  if(file.size() <= 2 * 1024 * 1024 && file.open(QIODevice::ReadOnly)) {
    QDataStream input(&file); input.setVersion(QDataStream::Qt_5_0);
    quint32 signature; qint64 expiry;
    input >> signature >> expiry >> value->etag >> value->modified >> value->control >> value->expiresHeader >> value->png;
    if(input.status() == QDataStream::Ok && signature == magic) {
      value->image = decode(value->png); value->expires = QDateTime::fromMSecsSinceEpoch(expiry, Qt::UTC);
    }
    if(value->image.isNull()) *value = Entry{};
  }
  entries[key] = value; trimMemory(); return value;
}
QImage TileStore::image(const Tile &tile, QObject *viewer) {
  if(!valid || !validTile(tile)) return {};
  const auto value = entry(tile);
  const auto view = views.find(viewer);
  const bool validated = view != views.end() && view->second.validated.contains(tile.key());
  return validated || value->expires > settings.clock() ? value->image : QImage();
}
bool TileStore::needed(const QString &key) const { for(const auto &view: views) if(view.second.tiles.count(key)) return true; return false; }
void TileStore::replaceView(QObject *owner, const std::vector<Tile> &tiles) {
  if(!owner || !valid || tiles.size() > 64) return;
  if(!owners.contains(owner)) {
    owners.insert(owner); connect(owner, &QObject::destroyed, this, [this, owner] { owners.remove(owner); forgetView(owner); });
  }
  View view;
  for(const auto &tile: tiles) if(validTile(tile)) view.tiles[tile.key()] = tile;
  for(const auto &flight: flights) if(view.tiles.count(flight.first)) view.attempted.insert(flight.first);
  views[owner] = std::move(view);
  for(auto it = entries.begin(); it != entries.end();) {
    if(it->second->noStore && !needed(it->first)) it = entries.erase(it); else ++it;
  }
  // Stop transfers that no visible map needs after a pan, hide or layer toggle.
  std::vector<QPointer<QNetworkReply>> obsolete;
  for(const auto &flight: flights) if(!needed(flight.first)) obsolete.push_back(flight.second.reply);
  for(const auto &reply: obsolete) if(reply) reply->abort();
  pump();
}
void TileStore::forgetView(QObject *owner) {
  views.erase(owner);
  for(auto it = entries.begin(); it != entries.end();) {
    if(it->second->noStore && !needed(it->first)) it = entries.erase(it); else ++it;
  }
  std::vector<QPointer<QNetworkReply>> obsolete;
  for(const auto &flight: flights) if(!needed(flight.first)) obsolete.push_back(flight.second.reply);
  for(const auto &reply: obsolete) if(reply) reply->abort();
}
void TileStore::pump() {
  if(pumping || !valid || blocked || cooldown > settings.clock()) return;
  pumping = true;
  while(flights.size() < 2) {
    Tile next; QString key;
    for(auto &view: views) {
      for(const auto &candidate: view.second.tiles) {
        if(view.second.attempted.contains(candidate.first) || flights.count(candidate.first)) continue;
        const auto cached = entry(candidate.second);
        if(!cached->image.isNull() && cached->expires > settings.clock()) {
          view.second.attempted.insert(candidate.first); view.second.validated.insert(candidate.first); continue;
        }
        if(failedUntil[candidate.first] > settings.clock()) { view.second.attempted.insert(candidate.first); continue; }
        next = candidate.second; key = candidate.first; break;
      }
      if(!key.isEmpty()) break;
    }
    if(key.isEmpty()) break;
    // A persistent cache is required before contacting the public service.
    QDir().mkpath(settings.cacheDirectory); QTemporaryFile probe(QDir(settings.cacheDirectory).filePath("write-test-XXXXXX"));
    if(!probe.open()) { blocked = true; message = tr("OSM: tile cache is not writable"); emit changed(); break; }
    for(auto &view: views) if(view.second.tiles.count(key)) view.second.attempted.insert(key);
    const auto previous = entry(next);
    QNetworkRequest request(url(next));
    request.setRawHeader("User-Agent", "EDITools (+https://github.com/GoFEM/EDITools)");
    request.setRawHeader("Accept", "image/png");
    request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::SameOriginRedirectPolicy);
    request.setMaximumRedirectsAllowed(3);
    request.setAttribute(QNetworkRequest::Http2AllowedAttribute, true);
    request.setAttribute(QNetworkRequest::CookieLoadControlAttribute, QNetworkRequest::Manual);
    request.setAttribute(QNetworkRequest::CookieSaveControlAttribute, QNetworkRequest::Manual);
    request.setAttribute(QNetworkRequest::AuthenticationReuseAttribute, QNetworkRequest::Manual);
    if(!previous->image.isNull()) {
      if(!previous->etag.isEmpty()) request.setRawHeader("If-None-Match", previous->etag);
      if(!previous->modified.isEmpty()) request.setRawHeader("If-Modified-Since", previous->modified);
    }
    auto *reply = network->get(request); flights[key] = {next, reply};
    auto *timeout = new QTimer(reply); timeout->setSingleShot(true); timeout->start(20000);
    connect(timeout, &QTimer::timeout, reply, &QNetworkReply::abort);
    connect(reply, &QNetworkReply::downloadProgress, reply, [reply](qint64 bytes, qint64) { if(bytes > 1024 * 1024) reply->abort(); });
    connect(reply, &QNetworkReply::finished, this, [this, key, reply, previous] { finish(key, reply, previous); });
  }
  pumping = false;
}
void TileStore::finish(const QString &key, QNetworkReply *reply, const std::shared_ptr<Entry> &previous) {
  const auto found = flights.find(key); if(found == flights.end()) { reply->deleteLater(); return; }
  const auto tile = found->second.tile; flights.erase(found);
  const int status = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
  const auto now = settings.clock();
  bool success = false;
  if(reply->error() == QNetworkReply::NoError && (status == 200 || (status == 304 && !previous->image.isNull()))) {
    auto value = std::make_shared<Entry>(*previous);
    if(status == 200) { value->png = reply->readAll(); value->image = decode(value->png); value->control.clear(); value->expiresHeader.clear(); value->etag.clear(); value->modified.clear(); }
    if(!value->image.isNull()) {
      if(reply->hasRawHeader("ETag")) value->etag = reply->rawHeader("ETag");
      if(reply->hasRawHeader("Last-Modified")) value->modified = reply->rawHeader("Last-Modified");
      if(reply->hasRawHeader("Cache-Control")) value->control = reply->rawHeader("Cache-Control");
      if(reply->hasRawHeader("Expires")) value->expiresHeader = reply->rawHeader("Expires");
      const auto policy = cachePolicy(value->control, value->expiresHeader, reply->rawHeader("Date"), reply->rawHeader("Age"), now);
      for(auto &view: views) if(view.second.tiles.count(key)) view.second.validated.insert(key);
      value->noStore = policy.noStore;
      value->expires = policy.expires; value->touched = ++tick; entries[key] = value;
      if(!policy.noStore) save(tile, *value); else QFile::remove(cachePath(tile));
      success = true; failedUntil.erase(key); if(!blocked && cooldown <= now) message.clear();
    }
  }
  if(!success && (reply->error() != QNetworkReply::OperationCanceledError || needed(key))) {
    failedUntil[key] = now.addSecs(60);
    message = status ? tr("OSM unavailable (HTTP %1)").arg(status) : tr("OSM unavailable: network or TLS error");
    if(status == 0 || status == 200 || reply->error() == QNetworkReply::SslHandshakeFailedError) cooldown = now.addSecs(60);
    if(status == 401 || status == 403) blocked = true;
    if(status == 429 || status == 503) {
      bool ok; const auto seconds = reply->rawHeader("Retry-After").toLongLong(&ok);
      cooldown = ok ? now.addSecs(std::max<qint64>(60, seconds)) : httpDate(reply->rawHeader("Retry-After"));
      if(!cooldown.isValid() || cooldown < now.addSecs(60)) cooldown = now.addSecs(60);
    }
  }
  reply->deleteLater(); trimMemory(); emit changed();
  // No automatic retry timer: retry only on a new user view after the cooldown.
  QTimer::singleShot(0, this, [this] { pump(); });
}
void TileStore::save(const Tile &tile, const Entry &value) {
  QSaveFile file(cachePath(tile));
  if(!file.open(QIODevice::WriteOnly)) { blocked = true; message = tr("OSM: cannot save tile cache"); return; }
  QDataStream output(&file); output.setVersion(QDataStream::Qt_5_0);
  output << magic << value.expires.toMSecsSinceEpoch() << value.etag << value.modified << value.control << value.expiresHeader << value.png;
  if(output.status() != QDataStream::Ok || !file.commit()) { blocked = true; message = tr("OSM: cannot save tile cache"); }
  if(++saved % 32 == 0) trimDisk();
}
void TileStore::trimMemory() {
  while(entries.size() > 256) {
    auto oldest = entries.end();
    for(auto it = entries.begin(); it != entries.end(); ++it)
      if(!needed(it->first) && !flights.count(it->first) && (oldest == entries.end() || it->second->touched < oldest->second->touched)) oldest = it;
    if(oldest == entries.end()) break;
    entries.erase(oldest);
  }
}
void TileStore::trimDisk() {
  auto files = QDir(settings.cacheDirectory).entryInfoList({"*.tile"}, QDir::Files, QDir::Time | QDir::Reversed);
  qint64 bytes = 0; for(const auto &file: files) bytes += file.size();
  for(const auto &info: files) {
    if(bytes <= 512LL * 1024 * 1024) break;
    if(info.lastModified().addSecs(week) > settings.clock()) continue;
    QFile file(info.filePath()); if(!file.open(QIODevice::ReadOnly)) continue;
    QDataStream input(&file); input.setVersion(QDataStream::Qt_5_0); quint32 signature; qint64 expires; input >> signature >> expires;
    if(signature != magic || expires > settings.clock().toMSecsSinceEpoch()) continue;
    file.close(); if(QFile::remove(info.filePath())) bytes -= info.size();
  }
}
}
