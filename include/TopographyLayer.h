#ifndef EDITOOLS_TOPOGRAPHY_LAYER_H
#define EDITOOLS_TOPOGRAPHY_LAYER_H
#include "Topography.h"
#include "qcustomplot.h"

// Renders only local grid data; drawing/export never initiates a download.
class TopographyLayer : public QCPLayerable {
public:
  explicit TopographyLayer(QCustomPlot *plot);
  void setData(const Topography::State &state) { terrain = state; }
  void setCoordinates(const SurveyCoordinates &value, double units) {
    coordinates = value;
    scale = units;
  }

protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *) const override {}
  void draw(QCPPainter *painter) override;

private:
  Topography::State terrain;
  SurveyCoordinates coordinates;
  double scale = 1.;
  std::shared_ptr<const Topography::Grid> cachedGrid;
  QString signature;
  QImage image;
};
#endif
