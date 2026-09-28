#include "include/SurveyHistory.h"
#include <boost/archive/binary_oarchive.hpp>
#include <QCryptographicHash>
#include <sstream>

namespace {
QByteArray fingerprint(const MTSurveyData &survey) {
  std::ostringstream bytes(std::ios::binary);
  {
    boost::archive::binary_oarchive archive(bytes);
    archive << survey.coordinates() << survey.topography().style;
    for (const auto &name : survey.get_stations_names())
      archive << name << survey.get_station_data(name);
  }
  const auto text = bytes.str();
  return QCryptographicHash::hash(QByteArray(text.data(), int(text.size())), QCryptographicHash::Sha256);
}
} // namespace
struct SurveyHistory::State {
  MTSurveyData data;
  QByteArray fingerprint;
  State(const MTSurveyData &data, QByteArray fingerprint) : data(data), fingerprint(std::move(fingerprint)) {}
};
class SurveyHistory::Command : public QUndoCommand {
public:
  Command(SurveyHistory &owner, std::shared_ptr<const State> before, std::shared_ptr<const State> after,
          const QString &text)
      : QUndoCommand(text), owner(owner), before(std::move(before)), after(std::move(after)) {}
  void undo() override { owner.restore(before); }
  void redo() override {
    if (first)
      first = false;
    else
      owner.restore(after);
  }

private:
  SurveyHistory &owner;
  std::shared_ptr<const State> before, after;
  bool first = true;
};
SurveyHistory::SurveyHistory(std::function<void()> callback) : restored(std::move(callback)) {
  commands.setUndoLimit(20);
}
void SurveyHistory::reset(std::shared_ptr<MTSurveyData> data) {
  commands.clear();
  survey = std::move(data);
  initial = current = survey ? std::make_shared<State>(*survey, fingerprint(*survey)) : nullptr;
}
bool SurveyHistory::record(const QString &description) {
  if (!survey || !current)
    return false;
  auto digest = fingerprint(*survey);
  if (digest == current->fingerprint && survey->topography().grid == current->data.topography().grid)
    return false;
  auto next = std::make_shared<State>(*survey, std::move(digest));
  auto before = current;
  current = next;
  commands.push(new Command(*this, std::move(before), std::move(next), description));
  return true;
}
void SurveyHistory::restore(const std::shared_ptr<const State> &state) {
  survey->restore(state->data);
  current = state;
  if (restored)
    restored();
}
std::shared_ptr<MTSurveyData> SurveyHistory::original() const {
  return initial ? std::make_shared<MTSurveyData>(initial->data) : nullptr;
}
