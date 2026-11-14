#pragma once
#include "core/AutoMask.h"
#include <QObject>
#include <memory>
class QAction;
class QLabel;
class QMainWindow;
class QTimer;
class MTDataPlot;
class SurveyHistory;
class AutoMaskSettings;
class WorkspaceView;

// Owns auto-mask review state; numerical checks and edit plans live in AutoMask.
class AutoMaskWorkflow : public QObject {
  Q_OBJECT
public:
  AutoMaskWorkflow(QMainWindow &window, SurveyHistory &history,
                   const std::vector<std::shared_ptr<MTDataPlot>> &plots, QAction *linkTensorMasks);
  WorkspaceView *settingsView() const { return settingsWorkspace; }
  AutoMask::Options currentOptions() const { return options; }
  bool pairDiagonals() const;
  void setData(std::shared_ptr<MTSurveyData> survey, const QString &station);
  void setSelection(const std::set<std::string> &stations, const std::vector<std::string> &order);
signals:
  void masksApplied();
  void settingsRequested();
  void reviewRequested();
  void navigateRequested(int direction);
private:
  void refresh();
  void analyseCurrent();
  void updateActions();
  void syncOptions();
  void applyCurrent(bool next);
  void applySelected();
  void review();
  void paint(const AutoMask::Edits &edits);
  QMainWindow &window;
  SurveyHistory &history;
  std::vector<std::shared_ptr<MTDataPlot>> plots;
  std::shared_ptr<MTSurveyData> survey;
  std::string current;
  std::set<std::string> selected;
  std::vector<std::string> stationOrder;
  AutoMask::Options options;
  AutoMaskSettings *settingsEditor = nullptr;
  WorkspaceView *settingsWorkspace = nullptr;
  AutoMask::Edits pending;
  QAction *preview, *settings, *details, *previous, *next, *apply, *applyNext, *batch, *pair, *link;
  std::array<QAction *, AutoMask::criterionCount> criteria;
  std::array<QAction *, 4> components;
  QLabel *status;
  QTimer *timer;
  unsigned revision = 0;
  bool busy = false;
};
