#ifndef EDITOOLS_PROFILE_PLOTS_H
#define EDITOOLS_PROFILE_PLOTS_H
#include "include/MTMapData.h"
#include "include/ProfileAnalysis.h"
#include "qcustomplot.h"

struct ProfileGlyph {
  QString station;
  double along, period, left, right, logBottom, logTop, value;
  MTMapData::PhaseTensor tensor{};
};

class ProfileSectionPlot : public QCustomPlot {
public:
  explicit ProfileSectionPlot(QWidget *parent = nullptr);
  void setSamples(std::vector<ProfileGlyph> samples, bool frequency, bool ellipses, bool cells,
                  double glyphSize, const QString &label, QCPColorGradient colors, QCPRange range,
                  bool logarithmic);

private:
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
  void setPanels(std::vector<Panel> panels, const QString &title);

private:
  class Canvas;
  Canvas *canvas;
};
#endif
