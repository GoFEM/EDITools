#ifndef FIT_STATISTICS_WINDOW_H
#define FIT_STATISTICS_WINDOW_H

#include "include/FitStatistics.h"
#include "qcustomplot.h"
#include <QDialog>
#include <QTreeWidget>
#include <QComboBox>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <functional>
#include <memory>

class FitStatisticsWindow : public QDialog
{
public:
  FitStatisticsWindow(QWidget *parent, std::function<void()> refresh);
  void setData(const std::shared_ptr<MTSurveyData> &survey,
               const std::map<std::string, MTSurveyData> &responses);

private:
  struct Response {
    QString path, name;
    QColor color;
    FitStatistics::Result fit;
  };
  void updatePlots();
  void updateDetails();
  void savePdf();
  void editRanges();
  void applyAxisRanges();
  QCPColorGradient spatialGradient() const;
  std::vector<int> enabled() const;
  QCustomPlot *makePlot(const QString &name, const QString &title, const QString &x, const QString &y);
  void drawMap(QCustomPlot *plot, QCPColorScale *scale, int response, const QCPRange &range);
  void drawHeatmap(QCustomPlot *plot, QCPColorScale *scale, int response, const QCPRange &range);

  struct AxisOptions {
    bool automatic = true;
    QCPRange range;
  };
  struct ColorLimits {
    QCheckBox *automatic;
    QDoubleSpinBox *minimum, *maximum;
  };

  std::shared_ptr<MTSurveyData> survey;
  std::vector<Response> responses;
  QTreeWidget *responseList;
  QComboBox *errorSource, *detailA, *detailB;
  QComboBox *curveColors, *colorMap;
  QCheckBox *reverseColors;
  QCheckBox *stationNames;
  ColorLimits mapLimits, heatLimits;
  std::map<QString, std::array<AxisOptions, 2>> axisOptions;
  QLabel *summary;
  QCustomPlot *overall, *periods, *histogram, *components, *stations;
  QCustomPlot *mapA, *mapB, *heatA, *heatB;
  QCPColorScale *mapScaleA, *mapScaleB, *heatScaleA, *heatScaleB;
  std::vector<std::string> stationOrder;
  std::map<std::string, std::array<double, 3>> positions;
  QString mapXLabel, mapYLabel;
};
#endif
