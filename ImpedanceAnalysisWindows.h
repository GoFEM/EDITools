#ifndef EDITOOLS_IMPEDANCE_ANALYSIS_WINDOWS_H
#define EDITOOLS_IMPEDANCE_ANALYSIS_WINDOWS_H

#include "include/WorkspaceView.h"
#include <memory>
#include <functional>
#include "include/ImpedanceAnalysis.h"
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
  void showEvent(QShowEvent *event) override;
private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};
#endif
