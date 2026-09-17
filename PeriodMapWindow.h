#ifndef PERIOD_MAP_WINDOW_H
#define PERIOD_MAP_WINDOW_H

#include "include/MTMapData.h"
#include "include/MTSurveyData.h"
#include "qcustomplot.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <memory>
#include <functional>
#include <set>

class QPushButton;

class PeriodMapWindow : public QDialog
{
public:
  explicit PeriodMapWindow(QWidget *parent, std::function<void()> observationsChanged = {});
  void setData(const std::shared_ptr<MTSurveyData> &survey,
               const std::map<std::string, MTSurveyData> &responses);

private:
  struct Site {
    std::string name;
    QPointF position, north, east; // Projected position in km and local unit directions.
    QString details;
  };
  const MTSurveyData *source() const;
  void updatePeriods();
  void updateMap(bool fit = false);
  void fitView();
  void resetColorRange();
  void savePdf();
  void drawEllipse(const Site &site, const MTMapData::PhaseTensor &tensor, double period,
                   QCPColorGradient &colors, const QCPRange &range);
  void drawArrow(const Site &site, const std::array<double, 2> &vector, bool imaginary, double period);
  void extendBounds(const QPointF &point);
  void selectStations(const QRect &rectangle, bool multiple, bool drag);
  void updateSelection(bool replot = true);
  void setSelectedMasks(bool enabled);
  int observedIndex(const std::string &name) const;

  std::shared_ptr<MTSurveyData> survey;
  std::map<QString, const MTSurveyData *> responses;
  std::vector<Site> sites;
  QString coordinatesDescription, locationError;
  QComboBox *dataset, *period, *convention, *colorBy, *colorMap;
  QCheckBox *names, *tensors, *realArrows, *imagArrows, *reverseColors;
  QDoubleSpinBox *tolerance, *ellipseSize, *arrowSize, *colorMin, *colorMax;
  QCustomPlot *plot;
  QCPColorScale *colorScale;
  QCPTextElement *title, *note;
  QLabel *summary;
  QCPRange boundsX, boundsY;
  bool haveBounds = false, geometryChanged = true;
  QCPItemLine *referenceArrow = nullptr;
  std::function<void()> observationsChanged;
  std::set<std::string> selectedStations;
  QCheckBox *selectMode;
  QComboBox *maskGroup;
  QPushButton *maskButton, *unmaskButton, *clearSelection;
  QLabel *selectionSummary;
  QCPGraph *selectionGraph = nullptr;
  QPoint selectionStart;
  bool selectionStarted = false, selectionDragged = false;
};
#endif
