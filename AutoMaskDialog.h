#pragma once
#include "include/AutoMask.h"
class QWidget;

// Review a single analysis pass, then apply checked candidates. Caller records history.
bool reviewAutoMasks(QWidget *parent, MTSurveyData &survey, const std::set<std::string> &selected,
                     bool linkTensorMasks, const ImpedanceAnalysis::RhoPlusOptions &rhoPlus = {},
                     AutoMask::Options *preferences = nullptr, bool *pairPreference = nullptr);

#include "include/WorkspaceView.h"
#include <memory>
#include <functional>
class AutoMaskReview : public WorkspaceView {
public:
  AutoMaskReview(QWidget *parent, std::shared_ptr<MTSurveyData> survey, std::set<std::string> selected,
                 bool link, const AutoMask::Options &options, bool paired, std::function<void()> changed);
  ~AutoMaskReview() override;
  void setData(std::shared_ptr<MTSurveyData> survey);
  AutoMask::Options options() const;
  bool paired() const;
  QJsonObject saveViewState() const override;
  void restoreViewState(const QJsonObject &state) override;
private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
