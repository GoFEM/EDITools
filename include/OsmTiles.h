#ifndef EDITOOLS_OSM_TILES_H
#define EDITOOLS_OSM_TILES_H

#include <QObject>
#include <QDateTime>
#include <QImage>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QPointer>
#include <QSet>
#include <functional>
#include <map>
#include <memory>

namespace Osm {
struct Tile {
  int z = 0, x = 0, y = 0;
  QString key() const;
  bool operator<(const Tile &other) const;
};
struct Config {
  QString url = "https://tile.openstreetmap.org/{z}/{x}/{y}.png";
  QString cacheDirectory, extraAttribution;
  std::function<QDateTime()> clock = [] { return QDateTime::currentDateTimeUtc(); };
};
struct CachePolicy {
  QDateTime expires;
  bool noStore = false;
};
CachePolicy cachePolicy(const QByteArray &control, const QByteArray &expires,
                        const QByteArray &date, const QByteArray &age, const QDateTime &now);

// All network access is initiated by replaceView, never by image/render/export.
// Each owner supplies only tiles intersecting its current, visible viewport.
class TileStore : public QObject {
  Q_OBJECT
public:
  explicit TileStore(const Config &config, QObject *parent = nullptr, QNetworkAccessManager *transport = nullptr);
  static TileStore *shared();
  const Config &config() const { return settings; }
  QImage image(const Tile &tile, QObject *viewer = nullptr); // Cache-only; expired entries need validation for this view.
  void replaceView(QObject *owner, const std::vector<Tile> &tiles);
  void forgetView(QObject *owner);
  QString status() const { return message; }
  static int opacity();
  static void setOpacity(int percent);
  QString cachePath(const Tile &tile) const;
signals:
  void changed();
private:
  struct Entry {
    QImage image;
    QByteArray png, etag, modified, control, expiresHeader;
    QDateTime expires;
    quint64 touched = 0;
    bool noStore = false;
  };
  struct View { std::map<QString, Tile> tiles; QSet<QString> attempted, validated; };
  struct Flight { Tile tile; QPointer<QNetworkReply> reply; };
  std::shared_ptr<Entry> entry(const Tile &tile);
  QUrl url(const Tile &tile) const;
  bool needed(const QString &key) const;
  void pump();
  void finish(const QString &key, QNetworkReply *reply, const std::shared_ptr<Entry> &previous);
  void save(const Tile &tile, const Entry &entry);
  void trimMemory();
  void trimDisk();
  Config settings;
  QNetworkAccessManager *network;
  std::map<QObject *, View> views;
  QSet<QObject *> owners;
  std::map<QString, std::shared_ptr<Entry>> entries;
  std::map<QString, Flight> flights;
  std::map<QString, QDateTime> failedUntil;
  QDateTime cooldown;
  QString message;
  bool valid = true, blocked = false, pumping = false;
  quint64 tick = 0;
  unsigned saved = 0;
};
}
#endif
