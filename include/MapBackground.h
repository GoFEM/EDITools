#ifndef EDITOOLS_MAP_BACKGROUND_H
#define EDITOOLS_MAP_BACKGROUND_H

#include "SurveyCoordinates.h"
#include "qcustomplot.h"
#include <functional>
#include <QToolButton>

namespace Osm { class Basemap; }
// A nonselectable layer below data/grid. Survives clearItems/clearPlottables and
// does not participate in axis rescaling. Owned by its QCustomPlot.
class MapBackground : public QCPLayerable {
public:
  enum Layer { Coast = 1, Countries = 2, States = 4, Water = 8, OpenStreetMap = 16 };
  explicit MapBackground(QCustomPlot *plot);
  void setCoordinates(const SurveyCoordinates &coordinates, double unitsPerMetre = 1.);
  void setLayers(int layers);
  int layers() const { return enabled; }
  static QString dataDirectory();
  static int availableLayers();
  static QToolButton *button(QWidget *parent, int initial, const std::function<void(int)> &changed, bool reportOnly = false);
protected:
  QRect clipRect() const override;
  void applyDefaultAntialiasingHint(QCPPainter *painter) const override;
  void draw(QCPPainter *painter) override;
private:
  struct Segment { QPolygonF points; QRectF bounds; };
  void prepare(int dataset);
  Osm::Basemap *osm;
  SurveyCoordinates coordinates;
  double scale = 1.;
  int enabled = 0;
  std::array<bool, 5> ready{{false, false, false, false, false}};
  std::array<std::vector<Segment>, 5> projected;
};
#endif
