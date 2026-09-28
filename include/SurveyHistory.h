#ifndef EDITOOLS_SURVEY_HISTORY_H
#define EDITOOLS_SURVEY_HISTORY_H
#include "MTSurveyData.h"
#include <QUndoStack>
#include <QByteArray>
#include <functional>

// Session-local survey history. Immutable DEMs and resampling sources remain shared.
// Call reset on loading a survey and record after each completed edit transaction.
class SurveyHistory {
public:
  explicit SurveyHistory(std::function<void()> restored = {});
  void reset(std::shared_ptr<MTSurveyData> survey);
  bool record(const QString &description);
  QUndoStack &stack() { return commands; }
  std::shared_ptr<MTSurveyData> original() const;

private:
  struct State;
  class Command;
  void restore(const std::shared_ptr<const State> &state);
  QUndoStack commands;
  std::shared_ptr<MTSurveyData> survey;
  std::shared_ptr<const State> initial, current;
  std::function<void()> restored;
};
#endif
