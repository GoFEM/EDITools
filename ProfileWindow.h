#ifndef EDITOOLS_PROFILE_WINDOW_H
#define EDITOOLS_PROFILE_WINDOW_H
#include "include/MTSurveyData.h"
#include "include/WorkspaceView.h"
#include <functional>
#include <memory>

class ProfileWindow : public WorkspaceView {
public:
  explicit ProfileWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> extracted);
  ~ProfileWindow();
  void setData(const std::shared_ptr<MTSurveyData> &survey,
               const std::map<std::string, MTSurveyData> &responses);
  void setMapLayers(int layers);
  QJsonObject saveViewState() const override;
  void restoreViewState(const QJsonObject &state) override;

private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
#endif
