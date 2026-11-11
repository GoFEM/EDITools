#ifndef EDITOOLS_STATION_COMPARISON_WINDOW_H
#define EDITOOLS_STATION_COMPARISON_WINDOW_H
#include "include/WorkspaceView.h"
#include "include/SurveyCoordinates.h"
#include "include/SharedMapInteraction.h"
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
  void setSharedMapActive(bool active, const SurveyCoordinates &coordinates) override;
  void setData(std::shared_ptr<MTSurveyData> survey, const QStringList &selected = {});
  void selectStation(const QString &name) override;

private:
  void refreshData(const std::shared_ptr<MTSurveyData> &survey, const Responses &) override { setData(survey); }
  void rebuild();
  void selectNearest();
  void draw();
  void updateSharedMap();
  void highlightStation();
  class Overlay;
  QPointer<Overlay> overlay;
  QPointer<QCustomPlot> sharedMap;
  SurveyCoordinates sharedCoordinates;
  SharedMapInteraction sharedInteraction;
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
