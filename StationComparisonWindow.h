#ifndef EDITOOLS_STATION_COMPARISON_WINDOW_H
#define EDITOOLS_STATION_COMPARISON_WINDOW_H
#include <QDialog>
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
class StationComparisonWindow : public QDialog {
public:
  explicit StationComparisonWindow(QWidget *parent = nullptr);
  void setData(std::shared_ptr<MTSurveyData> survey, const QStringList &selected = {});

private:
  void rebuild();
  void selectNearest();
  void draw();
  std::shared_ptr<MTSurveyData> survey;
  QComboBox *anchor, *family, *component;
  QSpinBox *count;
  QListWidget *stations;
  QCheckBox *errors;
  QLabel *status;
  QCustomPlot *map, *first, *second;
  MapBackground *background;
};
#endif
