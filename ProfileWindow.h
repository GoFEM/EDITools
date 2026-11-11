#ifndef EDITOOLS_PROFILE_WINDOW_H
#define EDITOOLS_PROFILE_WINDOW_H
#include "include/MTSurveyData.h"
#include <QDialog>
#include <functional>
#include <memory>

class ProfileWindow : public QDialog {
public:
  explicit ProfileWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> extracted);
  ~ProfileWindow();
  void setData(const std::shared_ptr<MTSurveyData> &survey,
               const std::map<std::string, MTSurveyData> &responses);
  void setMapLayers(int layers);

private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
#endif
