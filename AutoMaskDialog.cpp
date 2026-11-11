#include "AutoMaskDialog.h"
#include "AutoMaskSettings.h"
#include "GroupMaskPlots.h"
#include "include/MTSurveyData.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QHeaderView>
#include <QLabel>
#include <QMessageBox>
#include <QProgressDialog>
#include <QPushButton>
#include <QSaveFile>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>
#include <boost/archive/binary_oarchive.hpp>
#include <sstream>

namespace {
// Compare just station state: DEM assets and nested resampling snapshots do not
// affect these criteria. Kept in memory only, never written into an audit file.
std::string state(const MTSurveyData &survey) {
  std::ostringstream out; boost::archive::binary_oarchive archive(out);
  for(const auto &name : survey.get_stations_names()) archive << name << survey.get_station_data(name);
  return out.str();
}
}
bool reviewAutoMasks(QWidget *parent, MTSurveyData &survey, const std::set<std::string> &selected,
                     bool linkTensorMasks, const ImpedanceAnalysis::RhoPlusOptions &rhoPlus, AutoMask::Options *preferences, bool *pairPreference)
{
  QDialog dialog(parent); dialog.setObjectName("autoMaskDialog"); dialog.setWindowTitle(QObject::tr("Auto masking"));
  auto *outer = new QVBoxLayout(&dialog);
  auto *splitter = new QSplitter(&dialog); outer->addWidget(splitter, 1);
  auto *left = new QWidget(splitter); auto *layout = new QVBoxLayout(left); layout->setContentsMargins(0, 0, 0, 0);
  auto *scope = new QComboBox(left); scope->setObjectName("autoMaskScope");
  scope->addItems({QObject::tr("Selected sites (%1)").arg(selected.size()), QObject::tr("All enabled sites")});
  scope->setCurrentIndex(selected.empty() ? 1 : 0); layout->addWidget(scope);
  AutoMask::Options initial = preferences ? *preferences : AutoMask::Options{}; initial.rhoPlus = rhoPlus;
  auto *settings = new AutoMaskSettings(initial, left); layout->addWidget(settings);
  auto *actions = new QHBoxLayout;
  auto *analyse = new QPushButton(QObject::tr("Analyze"), left); analyse->setObjectName("autoMaskAnalyze"); actions->addWidget(analyse);
  auto *flaggedOnly = new QCheckBox(QObject::tr("Flagged only"), left); flaggedOnly->setObjectName("autoMaskFlaggedOnly"); flaggedOnly->setChecked(true); actions->addWidget(flaggedOnly);
  auto *exportAudit = new QPushButton(QObject::tr("CSV…"), left); exportAudit->setObjectName("autoMaskExport"); exportAudit->setEnabled(false); actions->addWidget(exportAudit); layout->addLayout(actions);
  auto *table = new QTableWidget(left); table->setObjectName("autoMaskCandidates"); table->setColumnCount(5 + AutoMask::criterionCount);
  table->setHorizontalHeaderLabels({QObject::tr("Use"), QObject::tr("Station"), QObject::tr("Component"), QObject::tr("Period (s)"), "ρ+ σ", QObject::tr("Period σ"), QObject::tr("Spatial σ"), QObject::tr("Error / limit"), QObject::tr("Flags")});
  table->setEditTriggers(QAbstractItemView::NoEditTriggers); table->verticalHeader()->hide(); table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents); table->horizontalHeader()->setStretchLastSection(true);
  table->setToolTip(QObject::tr("Check points to include them in the mask plan. Coloured scores triggered a criterion; hover for support and reasons. A dash means unavailable, not a passed test. Flagged only includes points with any flag, even below the minimum flag count.")); layout->addWidget(table, 1);
  auto *links = new QHBoxLayout;
  auto *paired = new QCheckBox(QObject::tr("Pair diagonals"), left); paired->setObjectName("autoMaskPair"); paired->setChecked(pairPreference ? *pairPreference : true);
  paired->setToolTip(QObject::tr("Also mask Zxx for each checked Zxy point, and Zyy for Zyx. These are associated edits, not independent detections.")); links->addWidget(paired);
  auto *linked = new QCheckBox(QObject::tr("Link Z / PT"), left); linked->setObjectName("autoMaskLink"); linked->setChecked(linkTensorMasks);
  linked->setToolTip(QObject::tr("Also mask all four PT entries at each checked period. Does not recursively mask other Z components. Tipper candidates affect only their own component; this link applies only to impedance candidates.")); links->addWidget(linked); layout->addLayout(links);
  auto *summary = new QLabel(QObject::tr("Analyze to preview candidates. No masks change until Apply."), left); summary->setWordWrap(true); summary->setObjectName("autoMaskSummary"); layout->addWidget(summary);
  GroupMaskPlots *plots = new GroupMaskPlots(splitter, survey, {}); splitter->setStretchFactor(1, 1); splitter->setSizes({560, 850});
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  auto *apply = buttons->button(QDialogButtonBox::Ok); apply->setText(QObject::tr("Apply checked")); apply->setObjectName("autoMaskApply"); apply->setEnabled(false); outer->addWidget(buttons);
  AutoMask::Result result; AutoMask::Options used; used.rhoPlus = rhoPlus;
  std::unique_ptr<MTSurveyData> snapshot;
  std::vector<bool> checked;
  GroupMaskEdits edits;
  std::string analysedState;
  const auto plan = [&] {
    edits = snapshot ? AutoMask::maskPlan(*snapshot, result.points, checked, paired->isChecked(), linked->isChecked()) : AutoMask::Edits{};
    const auto chosen = std::count(checked.begin(), checked.end(), true);
    unsigned total = 0; for(const auto &entry : edits) total += entry.second.size();
    std::array<unsigned, AutoMask::criterionCount> supported{};
    for(const auto &point : result.points) for(unsigned k = 0; k < supported.size(); ++k) supported[k] += std::isfinite(point.scores[k]);
    plots->setEdits(edits, false); apply->setEnabled(total > 0);
    summary->setText(QObject::tr("%1 checked / %2 assessed · %3 mask changes\nSupported: ρ+ %4 · period %5 · spatial %6 · uncertainty %7 · %8 disabled sites skipped")
      .arg(chosen).arg(result.points.size()).arg(total).arg(supported[0]).arg(supported[1]).arg(supported[2]).arg(supported[3]).arg(result.disabled));
  };
  const auto populate = [&] {
    const QSignalBlocker block(table); table->setRowCount(0);
    for(unsigned i = 0; i < result.points.size(); ++i) {
      const auto &point = result.points[i]; if(flaggedOnly->isChecked() && !point.votes()) continue;
      const int row = table->rowCount(); table->insertRow(row);
      auto *use = new QTableWidgetItem; use->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable);
      use->setCheckState(checked[i] ? Qt::Checked : Qt::Unchecked); use->setData(Qt::UserRole, i); table->setItem(row, 0, use);
      table->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(point.station)));
      table->setItem(row, 2, new QTableWidgetItem(AutoMask::componentName(point.component)));
      table->setItem(row, 3, new QTableWidgetItem(QString::number(1./point.frequency, 'g', 12)));
      for(unsigned k = 0; k < AutoMask::criterionCount; ++k) {
        auto *cell = new QTableWidgetItem(std::isfinite(point.scores[k]) ? QString::number(point.scores[k], 'g', 3) : QStringLiteral("—"));
        cell->setToolTip(QString::fromStdString(point.details[k])); if(point.flagged[k]) cell->setBackground(QColor("#fee2e2")); table->setItem(row, k+4, cell);
      }
      table->setItem(row, 4 + AutoMask::criterionCount, new QTableWidgetItem(QString::number(point.votes())));
    }
  };
  const auto invalidate = [&] {
    const QSignalBlocker block(table);
    result = {}; checked.clear(); edits.clear(); analysedState.clear(); table->setRowCount(0);
    plots->setEdits({}, false); apply->setEnabled(false); exportAudit->setEnabled(false);
    summary->setText(QObject::tr("Settings changed. Analyze to refresh candidates."));
  };
  const auto enableSettings = [&] { const auto o = settings->options(); analyse->setEnabled(std::any_of(o.criteria.begin(), o.criteria.end(), [](bool v) { return v; })); };
  QObject::connect(settings, &AutoMaskSettings::changed, &dialog, [&] { invalidate(); enableSettings(); });
  enableSettings();
  QObject::connect(scope, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, invalidate);
  for(auto *box : {paired, linked}) QObject::connect(box, &QCheckBox::toggled, &dialog, plan);
  QObject::connect(flaggedOnly, &QCheckBox::toggled, &dialog, populate);
  QObject::connect(table, &QTableWidget::itemChanged, &dialog, [&](QTableWidgetItem *item) {
    if(item->column() == 0) { checked.at(item->data(Qt::UserRole).toUInt()) = item->checkState() == Qt::Checked; plan(); }
  });
  QObject::connect(analyse, &QPushButton::clicked, &dialog, [&] {
    invalidate();
    std::set<std::string> targets = selected;
    if(scope->currentIndex() == 1) { targets.clear(); for(const auto &name : survey.get_stations_names()) if(survey.is_active(name)) targets.insert(name); }
    used = settings->options();
    if(targets.empty() || std::none_of(used.components.begin(), used.components.end(), [](bool v) { return v; })) { summary->setText(QObject::tr("Choose sites and components.")); return; }
    if(unsigned(std::count(used.criteria.begin(), used.criteria.end(), true)) < used.minimumFlags) {
      summary->setText(QObject::tr("Enable enough criteria or lower Minimum flags explicitly.")); return;
    }
    snapshot.reset(new MTSurveyData(survey)); analysedState = state(*snapshot);
    QProgressDialog progress(QObject::tr("Analyzing stations…"), QObject::tr("Cancel"), 0, targets.size(), &dialog);
    progress.setObjectName("autoMaskProgress"); progress.setWindowModality(Qt::ApplicationModal); progress.setMinimumDuration(0); progress.setAutoClose(false); progress.setAutoReset(false);
    try {
      result = AutoMask::analyse(*snapshot, targets, used, [&](unsigned done, unsigned total) {
        progress.setMaximum(total); progress.setValue(done); QApplication::processEvents(); return !progress.wasCanceled();
      });
    } catch(const std::exception &e) { invalidate(); summary->setText(QString::fromUtf8(e.what())); return; }
    progress.close();
    if(result.cancelled) { invalidate(); summary->setText(QObject::tr("Analysis cancelled. No masks changed.")); return; }
    for(const auto &point : result.points) checked.push_back(point.candidate);
    auto *replacement = new GroupMaskPlots(nullptr, *snapshot, targets); splitter->replaceWidget(1, replacement); delete plots; plots = replacement;
    populate(); plan(); exportAudit->setEnabled(!result.points.empty());
  });
  QObject::connect(exportAudit, &QPushButton::clicked, &dialog, [&] {
    const auto path = QFileDialog::getSaveFileName(&dialog, QObject::tr("Export auto-mask audit"), "auto-mask.csv", QObject::tr("CSV (*.csv)")); if(path.isEmpty()) return;
    QSaveFile file(path); if(!file.open(QIODevice::WriteOnly)) { QMessageBox::warning(&dialog, QObject::tr("Export"), file.errorString()); return; }
    QTextStream out(&file); out.setEncoding(QStringConverter::Utf8); out.setRealNumberPrecision(17);
    out << "# minimum_flags=" << used.minimumFlags << ",max_gap_ratio=" << used.maximumGapRatio << ",neighbours=" << used.neighbours << ",radius_km=" << used.radiusKm
        << ",period_tolerance=" << used.periodTolerance << ",min_phase_degrees=" << used.minimumPhaseDifference << ",agreement_degrees=" << used.phaseAgreement << ",rho_grid=" << used.rhoPlus.intervalsPerDecade
        << ",max_tipper_error=" << used.maximumTipperError << ",min_tipper_difference=" << used.minimumTipperDifference << ",tipper_agreement=" << used.tipperAgreement
        << ",max_relative_rho_error=" << used.maximumRelativeRhoError << ",max_phase_error_degrees=" << used.maximumPhaseError
        << ",pair_diagonals=" << paired->isChecked() << ",link_z_pt=" << linked->isChecked() << '\n';
    for(unsigned k = 0; k < AutoMask::criterionCount; ++k) out << "# criterion_" << k << "_enabled=" << used.criteria[k] << (k == 3 ? ",threshold_multiplier=" : ",threshold_sigma=") << used.threshold[k] << '\n';
    out << "station,component,period_s,rho_score,period_score,spatial_score,uncertainty_ratio,rho_flag,period_flag,spatial_flag,uncertainty_flag,suggested,checked,rho_detail,period_detail,spatial_detail,uncertainty_detail\n";
    auto quote = [](const std::string &s) { return '"' + QString::fromStdString(s).replace("\"", "\"\"") + '"'; };
    for(unsigned i = 0; i < result.points.size(); ++i) {
      const auto &p = result.points[i]; out << quote(p.station) << ',' << AutoMask::componentName(p.component) << ',' << 1./p.frequency;
      for(double score : p.scores) out << ',' << (std::isfinite(score) ? QString::number(score, 'g', 17) : QString());
      for(bool flag : p.flagged) out << ',' << flag;
      out << ',' << p.candidate << ',' << checked[i]; for(const auto &detail : p.details) out << ',' << quote(detail); out << '\n';
    }
    out.flush(); if(out.status() != QTextStream::Ok || !file.commit()) QMessageBox::warning(&dialog, QObject::tr("Export"), QObject::tr("Could not save the complete audit."));
  });
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, [&] {
    if(edits.empty()) return;
    if(state(survey) != analysedState) { invalidate(); summary->setText(QObject::tr("Survey changed. Analyze again before applying.")); return; }
    dialog.accept();
  });
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  dialog.resize(1450, 850);
  const auto accepted = dialog.exec() == QDialog::Accepted;
  if(preferences) *preferences = settings->options();
  if(pairPreference) *pairPreference = paired->isChecked();
  if(!accepted || edits.empty() || state(survey) != analysedState) return false;
  AutoMask::apply(survey, edits);
  return true;
}
