#ifndef EDITOOLS_IMPEDANCE_ANALYSIS_WINDOWS_H
#define EDITOOLS_IMPEDANCE_ANALYSIS_WINDOWS_H

#include "ui/workspace/WorkspaceView.h"
#include <memory>
#include <functional>
#include "core/ImpedanceAnalysis.h"
class MTSurveyData;
class SurveyHistory;
class QShowEvent;

class ImpedanceValidityWindow : public WorkspaceView {
public:
  explicit ImpedanceValidityWindow(QWidget *parent = nullptr,
    std::function<void(const ImpedanceAnalysis::RhoPlusOptions &)> optionsChanged = {},
    std::function<void()> masksChanged = {});
  ~ImpedanceValidityWindow() override;
  void setData(const std::shared_ptr<MTSurveyData> &survey, const QString &station = {});
  void setHistory(SurveyHistory *history);
  void setFitOptions(const ImpedanceAnalysis::RhoPlusOptions &options);
  void setLinkTensorMasks(bool linked);
protected:
  void refreshData(const std::shared_ptr<MTSurveyData> &survey, const Responses &) override { setData(survey); }
  bool hasDeferredPlots() const override { return true; }
  void showEvent(QShowEvent *event) override;
private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
class SurveyAverageWindow : public WorkspaceView {
public:
  explicit SurveyAverageWindow(QWidget *parent = nullptr);
  ~SurveyAverageWindow() override;
  void setData(const std::shared_ptr<MTSurveyData> &survey);
protected:
  void refreshData(const std::shared_ptr<MTSurveyData> &survey, const Responses &) override { setData(survey); }
  bool hasDeferredPlots() const override { return true; }
  void showEvent(QShowEvent *event) override;
private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
#endif
