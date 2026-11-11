#ifndef EDITOOLS_STATION_COMPARISON_WINDOW_H
#define EDITOOLS_STATION_COMPARISON_WINDOW_H
#include "include/WorkspaceView.h"
#include "include/SurveyCoordinates.h"
#include "qcustomplot.h"
#include <memory>
#include <vector>
class MTSurveyData;
class QComboBox;
class QSpinBox;
class QListWidget;
class QCheckBox;
class QLabel;
class QCustomPlot;
class MapBackground;
class StationComparisonWindow : public WorkspaceView {
public:
  explicit StationComparisonWindow(QWidget *parent = nullptr, QCustomPlot *sharedMap = nullptr);
  ~StationComparisonWindow();
  void setSharedMapActive(bool active, const SurveyCoordinates &coordinates);
  void setData(std::shared_ptr<MTSurveyData> survey, const QStringList &selected = {});
  void selectReferenceStation(const QString &name);

private:
  void rebuild();
  void selectNearest();
  void draw();
  void updateSharedMap();
  void highlightStation();
  class Overlay;
  QPointer<Overlay> overlay;
  QPointer<QCustomPlot> sharedMap;
  SurveyCoordinates sharedCoordinates;
  bool sharedActive = false;
  QCP::Interactions savedInteractions;
  QPoint sharedPress;
  std::shared_ptr<MTSurveyData> survey;
  QComboBox *anchor, *family, *component;
  QSpinBox *count;
  QListWidget *stations;
  QCheckBox *errors;
  QLabel *status, *referenceLabel;
  QCustomPlot *map, *first, *second;
  MapBackground *background;
};
#endif
