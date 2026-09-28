#ifndef EDITOOLS_TOPOGRAPHY_DIALOG_H
#define EDITOOLS_TOPOGRAPHY_DIALOG_H
#include "include/MTSurveyData.h"
#include <QDialog>
#include <functional>

class TopographyDialog : public QDialog {
public:
  TopographyDialog(QWidget *parent, MTSurveyData &survey, std::function<void()> changed);
};
#endif
