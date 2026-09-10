#ifndef SURVEY_COORDINATES_DIALOG_H
#define SURVEY_COORDINATES_DIALOG_H
#include <QDialog>
#include "include/SurveyCoordinates.h"

class MTSurveyData;
class QSpinBox;
class QComboBox;
class QCheckBox;
class QLineEdit;
class QTableWidget;
class QLabel;
class QPushButton;

class SurveyCoordinatesDialog : public QDialog
{
public:
  explicit SurveyCoordinatesDialog(MTSurveyData &survey, QWidget *parent = nullptr);
private:
  void updatePreview();
  MTSurveyData &survey;
  SurveyCoordinates preview;
  QSpinBox *zone;
  QComboBox *hemisphere;
  QCheckBox *center;
  QLineEdit *easting, *northing;
  QTableWidget *table;
  QLabel *status;
  QPushButton *apply, *copy;
};
#endif
