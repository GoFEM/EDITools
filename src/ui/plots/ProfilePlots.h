#ifndef EDITOOLS_PROFILE_PLOTS_H
#define EDITOOLS_PROFILE_PLOTS_H
#include "core/MTMapData.h"
#include "core/ProfileAnalysis.h"
#include "qcustomplot.h"

struct ProfileGlyph {
  QString station;
  double along, period, left, right, logBottom, logTop, value;
  MTMapData::PhaseTensor tensor{};
};

class ProfileSectionPlot : public QCustomPlot {
  Q_OBJECT
public:
  explicit ProfileSectionPlot(QWidget *parent = nullptr);
  void setStations(const std::vector<ProfileAnalysis::Offset> &stations, double cellWidth);
  void setSelectedStation(const QString &name);
  void setReference(const QString &up);
  void setSamples(std::vector<ProfileGlyph> samples, bool frequency, bool ellipses, bool cells,
                  double glyphSize, const QString &label, QCPColorGradient colors, QCPRange range,
                  bool logarithmic);

signals:
  void stationClicked(const QString &name);

private:
  QPoint pressPosition;
  bool dragged = false;
  class Canvas;
  Canvas *canvas;
  QCPColorScale *scale;
};

class ProfileRosePlot : public QCustomPlot {
public:
  struct Panel {
    QString title;
    QColor color;
    ProfileAnalysis::Rose histogram;
  };
  explicit ProfileRosePlot(QWidget *parent = nullptr);
  void setPanels(std::vector<Panel> panels, const QString &title, bool strikeFrame = false);

private:
  class Canvas;
  Canvas *canvas;
};
#endif
