#ifndef EDITOOLS_OSM_BASEMAP_H
#define EDITOOLS_OSM_BASEMAP_H
#include "SurveyCoordinates.h"
#include "OsmTiles.h"
#include "qcustomplot.h"
#include <QTimer>

namespace Osm {
struct View {
  SurveyCoordinates coordinates;
  double unitsPerMetre = 1., left = -180., right = 180., top = 85., bottom = -85.;
  QSize size;
  QString signature() const;
};
struct Frame {
  QImage image;
  std::vector<Tile> tiles;
  unsigned loaded = 0;
};
// Inverse-project a geographic/UTM viewport into Web Mercator. No network calls.
Frame render(const View &view, TileStore &store, QObject *viewer = nullptr);

class Basemap : public QCPLayerable {
public:
  explicit Basemap(QCustomPlot *plot, TileStore *store = nullptr);
  ~Basemap() override;
  void setEnabled(bool enabled);
  void setCoordinates(const SurveyCoordinates &value, double unitsPerMetre);
  void drawCredit(QCPPainter *painter);
protected:
  bool eventFilter(QObject *object, QEvent *event) override;
  QRect clipRect() const override;
  void applyDefaultAntialiasingHint(QCPPainter *) const override {}
  void draw(QCPPainter *painter) override;
private:
  View currentView() const;
  bool canRequest() const;
  void scheduleView();
  void stopRequests();
  bool enabled = false;
  SurveyCoordinates coordinates;
  double units = 1.;
  TileStore *store;
  QTimer debounce;
  QString frameSignature, requestedSignature;
  QRect creditRect;
  QString creditTooltip;
  bool creditPressed = false;
  Frame frame;
  View screenView;
  bool haveScreenView = false;
};
}
#endif
