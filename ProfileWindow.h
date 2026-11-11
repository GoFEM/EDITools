#ifndef EDITOOLS_PROFILE_WINDOW_H
#define EDITOOLS_PROFILE_WINDOW_H
#include "include/MTSurveyData.h"
#include "include/WorkspaceView.h"
#include <functional>
#include <memory>

class QCustomPlot;

class ProfileWindow : public WorkspaceView {
public:
  explicit ProfileWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> extracted, QCustomPlot *sharedMap = nullptr);
  ~ProfileWindow();
  void setData(const std::shared_ptr<MTSurveyData> &survey,
               const std::map<std::string, MTSurveyData> &responses);
  void setMapLayers(int layers);
  void setSharedMapActive(bool active, const SurveyCoordinates &coordinates) override;
  void selectStation(const QString &name) override;
  QJsonObject saveViewState() const override;
  void restoreViewState(const QJsonObject &state) override;

private:
  void refreshData(const std::shared_ptr<MTSurveyData> &survey, const Responses &responses) override { setData(survey, responses); }
  struct Impl;
  std::unique_ptr<Impl> impl;
};
#endif
