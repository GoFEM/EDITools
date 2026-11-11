#include "AutoMaskWorkflow.h"
#include "AutoMaskDialog.h"
#include "AutoMaskSettings.h"
#include "include/WorkspaceView.h"
#include "include/MTDataPlot.h"
#include "include/MTSurveyData.h"
#include "include/SurveyHistory.h"
#include <QAction>
#include <QApplication>
#include <QLabel>
#include <QMainWindow>
#include <QProgressDialog>
#include <QSignalBlocker>
#include <QStatusBar>
#include <QTimer>

AutoMaskWorkflow::AutoMaskWorkflow(QMainWindow &window, SurveyHistory &history,
    const std::vector<std::shared_ptr<MTDataPlot>> &plots, QAction *linkTensorMasks)
  : QObject(&window), window(window), history(history), plots(plots), link(linkTensorMasks) {
  auto action = [&](const char *name, const QString &text, const QString &help, bool checkable = false) {
    auto *a = new QAction(text, &window); a->setObjectName(name); a->setToolTip(help); a->setCheckable(checkable); return a;
  };
  preview = action("actionAuto_mask_preview", tr("Preview"), tr("Show auto-mask candidates as red circles on the main plots. Refreshes when the station, data or settings change. Preview does not change masks."), true);
  settings = action("actionAuto_mask_settings", tr("Settings…"), tr("Thresholds, required flags and neighbour limits. Uses the same settings as detailed review."));
  details = action("actionAuto_masking", tr("Detailed review…"), tr("Review scores, check individual candidates and export an audit for selected sites."));
  previous = action("actionAuto_mask_previous", tr("Previous"), tr("Show the previous station without changing the batch selection or applying masks."));
  next = action("actionAuto_mask_next", tr("Next"), tr("Show the next station without changing the batch selection or applying masks."));
  apply = action("actionAuto_mask_apply", tr("Apply"), tr("Apply the current station's previewed masks as one Undo step."));
  applyNext = action("actionAuto_mask_apply_next", tr("Apply + next"), tr("Apply the current preview, then show the next station. Does not change the batch selection."));
  batch = action("actionAuto_mask_selected", tr("Apply to selected"), tr("Analyse all selected, enabled sites using the current checks and apply their candidates in one Undo step. Cancellation applies nothing. Unselected sites provide references only."));
  pair = action("actionAuto_mask_pair", tr("Pair diagonals"), tr("For impedance candidates, also mask Zxx for Zxy and Zyy for Zyx. Does not affect tipper candidates."), true); pair->setChecked(true);
  const auto names = AutoMaskSettings::criterionNames(), help = AutoMaskSettings::criterionHelp();
  for(unsigned k = 0; k < criteria.size(); ++k) {
    criteria[k] = action(("actionAuto_mask_check" + QString::number(k)).toUtf8().constData(), names[k], help[k], true);
    criteria[k]->setChecked(options.criteria[k]);
    connect(criteria[k], &QAction::toggled, this, [this, k](bool on) { options.criteria[k] = on; refresh(); });
  }
  const QStringList labels{"XY", "YX", "Tzx", "Tzy"};
  for(unsigned c = 0; c < components.size(); ++c) {
    components[c] = action(("actionAuto_mask_" + labels[c]).toUtf8().constData(), labels[c], tr("Include %1 in auto-masking. Independent of plot visibility; real and imaginary parts share a mask.").arg(labels[c]), true);
    components[c]->setChecked(options.components[c]);
    connect(components[c], &QAction::toggled, this, [this, c](bool on) { options.components[c] = on; refresh(); });
  }
  status = new QLabel(&window); status->setObjectName("autoMaskStatus"); window.statusBar()->addPermanentWidget(status);
  timer = new QTimer(this); timer->setSingleShot(true); timer->setInterval(120);
  connect(timer, &QTimer::timeout, this, &AutoMaskWorkflow::analyseCurrent);
  connect(preview, &QAction::toggled, this, [this] { refresh(); });
  for(auto *a : {pair, link}) connect(a, &QAction::toggled, this, [this] { refresh(); });
  connect(settings, &QAction::triggered, this, &AutoMaskWorkflow::settingsRequested);
  connect(details, &QAction::triggered, this, &AutoMaskWorkflow::review);
  connect(previous, &QAction::triggered, this, [this] { emit navigateRequested(-1); });
  connect(next, &QAction::triggered, this, [this] { emit navigateRequested(1); });
  connect(apply, &QAction::triggered, this, [this] { applyCurrent(false); });
  connect(applyNext, &QAction::triggered, this, [this] { applyCurrent(true); });
  connect(batch, &QAction::triggered, this, &AutoMaskWorkflow::applySelected);
  settingsWorkspace = new WorkspaceView(&window); settingsWorkspace->hide(); settingsWorkspace->setEmbedded();
  settingsEditor = new AutoMaskSettings(options, settingsWorkspace);
  settingsEditor->moveToRibbon(*settingsWorkspace, tr("Mask settings"));
  connect(settingsEditor, &AutoMaskSettings::changed, this, [this] { options = settingsEditor->options(); syncOptions(); refresh(); });
  settingsWorkspace->setStateRefresh([this] { emit settingsEditor->changed(); });

  refresh();
}
void AutoMaskWorkflow::setData(std::shared_ptr<MTSurveyData> data, const QString &station) {
  survey = std::move(data); current = station.toStdString(); refresh();
}
void AutoMaskWorkflow::setSelection(const std::set<std::string> &stations) { selected = stations; updateActions(); }
void AutoMaskWorkflow::syncOptions() {
  for(unsigned k = 0; k < criteria.size(); ++k) { const QSignalBlocker b(criteria[k]); criteria[k]->setChecked(options.criteria[k]); }
  for(unsigned c = 0; c < components.size(); ++c) { const QSignalBlocker b(components[c]); components[c]->setChecked(options.components[c]); }
}
void AutoMaskWorkflow::paint(const AutoMask::Edits &edits) {
  const auto found = edits.find(current);
  const std::vector<std::pair<RealDataType, double>> empty;
  for(const auto &plot : plots) plot->set_mask_preview(found == edits.end() ? empty : found->second);
}
void AutoMaskWorkflow::refresh() {
  if(settingsEditor) settingsEditor->setOptions(options);
  ++revision; timer->stop(); pending.clear(); paint({}); updateActions();
  status->setText(current.empty() ? QString() : preview->isChecked() ?
    tr("Auto mask · %1: pending").arg(QString::fromStdString(current)) : tr("Station: %1").arg(QString::fromStdString(current)));
  status->setToolTip(QString());
  if(preview->isChecked() && survey && survey->is_station_present(current) && !busy) timer->start();
}
void AutoMaskWorkflow::updateActions() {
  const bool hasData = survey && survey->n_stations(), hasCurrent = hasData && survey->is_station_present(current);
  preview->setEnabled(hasData && !busy); settings->setEnabled(!busy); details->setEnabled(hasData && !busy);
  for(auto *a : criteria) a->setEnabled(!busy);
  for(auto *a : components) a->setEnabled(!busy);
  pair->setEnabled(!busy);
  const auto names = hasData ? survey->get_stations_names() : std::vector<std::string>{};
  const auto it = std::find(names.begin(), names.end(), current);
  previous->setEnabled(!busy && hasCurrent && it != names.begin());
  next->setEnabled(!busy && hasCurrent && it != names.end() && std::next(it) != names.end());
  apply->setEnabled(!busy && preview->isChecked() && !pending.empty());
  applyNext->setEnabled(apply->isEnabled() && next->isEnabled());
  batch->setEnabled(!busy && hasData && std::any_of(selected.begin(), selected.end(), [this](const std::string &n) { return survey->is_station_present(n) && survey->is_active(n); }));
}
void AutoMaskWorkflow::analyseCurrent() {
  if(!preview->isChecked() || !survey || !survey->is_station_present(current)) return;
  try {
    const auto result = AutoMask::analyse(*survey, {current}, options);
    std::vector<bool> checked; QStringList reasons;
    for(const auto &p : result.points) {
      checked.push_back(p.candidate);
      if(p.candidate) {
        QStringList flags;
        for(unsigned k = 0; k < p.flagged.size(); ++k) if(p.flagged[k]) flags << AutoMaskSettings::criterionNames()[k];
        reasons << tr("%1 · %2 s: %3").arg(AutoMask::componentName(p.component)).arg(1./p.frequency, 0, 'g', 6).arg(flags.join(", "));
      }
    }
    pending = AutoMask::maskPlan(*survey, result.points, checked, pair->isChecked(), link->isChecked());
    paint(pending); status->setText(tr("Auto mask · %1: %n candidate(s)", nullptr, std::count(checked.begin(), checked.end(), true)).arg(QString::fromStdString(current)));
    status->setToolTip(reasons.isEmpty() ? tr("No candidates with the enabled checks. Unsupported tests provide no vote; open Detailed review for scores and reasons.") : reasons.join('\n'));
  } catch(const std::exception &e) { pending.clear(); paint({}); status->setText(tr("Auto mask: unavailable")); status->setToolTip(QString::fromUtf8(e.what())); }
  updateActions();
}
void AutoMaskWorkflow::applyCurrent(bool advance) {
  if(!apply->isEnabled() || !survey || pending.empty()) return;
  timer->stop(); const auto edits = pending; pending.clear();
  AutoMask::apply(*survey, edits); history.record(tr("Auto masking"));
  emit masksApplied();
  if(advance) emit navigateRequested(1);
}
void AutoMaskWorkflow::applySelected() {
  if(!batch->isEnabled() || !survey) return;
  timer->stop(); busy = true; updateActions();
  const auto generation = revision; const auto data = survey; const MTSurveyData snapshot(*data); const auto targets = selected;
  QProgressDialog progress(tr("Analyzing selected sites…"), tr("Cancel"), 0, targets.size(), &window);
  progress.setObjectName("autoMaskBatchProgress"); progress.setWindowModality(Qt::ApplicationModal); progress.setMinimumDuration(0);
  progress.setAutoClose(false); progress.setAutoReset(false);
  try {
    const auto result = AutoMask::analyse(snapshot, targets, options, [&](unsigned done, unsigned total) {
      progress.setMaximum(total); progress.setValue(done); QApplication::processEvents();
      return !progress.wasCanceled() && generation == revision;
    });
    progress.close();
    if(!result.cancelled && generation == revision) {
      std::vector<bool> checked; for(const auto &p : result.points) checked.push_back(p.candidate);
      const auto edits = AutoMask::maskPlan(snapshot, result.points, checked, pair->isChecked(), link->isChecked());
      if(!edits.empty()) { AutoMask::apply(*data, edits); history.record(tr("Auto masking selected sites")); emit masksApplied(); }
    }
  } catch(const std::exception &e) { window.statusBar()->showMessage(QString::fromUtf8(e.what()), 10000); }
  busy = false; refresh();
}
void AutoMaskWorkflow::review() { if(survey) emit reviewRequested(); }
bool AutoMaskWorkflow::pairDiagonals() const { return pair->isChecked(); }
