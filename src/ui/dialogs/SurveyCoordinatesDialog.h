#ifndef SURVEY_COORDINATES_DIALOG_H
#define SURVEY_COORDINATES_DIALOG_H
#include <QDialog>
#include "core/SurveyCoordinates.h"
#include "core/UTMAnalysis.h"

class MTSurveyData;
class QSpinBox;
class QComboBox;
class QCheckBox;
class QLineEdit;
class QTableWidget;
class QLabel;
class QPushButton;
class QCustomPlot;
class QDoubleSpinBox;

class SurveyCoordinatesDialog : public QDialog {
public:
  explicit SurveyCoordinatesDialog(MTSurveyData &survey, QWidget *parent = nullptr);

private:
  void updatePreview();
  void updateZoneOverlay();
  void fitSurvey();
  const UTM::Analysis &analysis(int selectedZone);
  MTSurveyData &survey;
  SurveyCoordinates preview;
  SurveyCoordinates suggested;
  std::vector<std::array<double, 3>> geographic;
  UTM::Extent bounds;
  std::map<int, UTM::Analysis> analyses;
  QSpinBox *zone;
  QComboBox *hemisphere;
  QCheckBox *center;
  QLineEdit *easting, *northing;
  QTableWidget *table;
  QTableWidget *candidates;
  QCustomPlot *map;
  QDoubleSpinBox *limit;
  QCheckBox *names;
  QLabel *crossing, *distortion;
  QLabel *status;
  QPushButton *apply, *copy;
};
#endif
