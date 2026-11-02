#include "EDIImportDialog.h"
#include "include/EDIPeriodMerge.h"
#include "include/PeriodResampling.h"
#include "include/ToolTips.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>
#include <cmath>

bool reviewEDIImport(QWidget *parent, MTSurveyData &candidate,
                     const std::set<std::string> &incoming, const std::vector<std::string> &duplicates)
{
  QDialog dialog(parent); dialog.setObjectName("ediImportDialog"); dialog.setWindowTitle(QObject::tr("Import EDIs"));
  auto *layout = new QVBoxLayout(&dialog);
  auto *summary = new QLabel(&dialog); summary->setObjectName("ediImportSummary"); summary->setWordWrap(true); layout->addWidget(summary);
  auto *controls = new QHBoxLayout;
  auto *mode = new QComboBox(&dialog); mode->setObjectName("ediPeriodMode");
  mode->addItems({QObject::tr("Keep original"), QObject::tr("Merge close periods"), QObject::tr("Interpolate to existing")});
  controls->addWidget(mode);
  ToolTips::set(mode, QObject::tr("Keep incoming frequencies, relabel nearby samples without interpolation, or resample only new stations onto the union of valid periods stored at existing stations, including disabled stations and masked periods. Existing stations are unchanged. Interpolation uses real/imaginary transfer functions versus log(period), never extrapolates and does not fill gaps above the chosen ratio."));
  auto *tolerance = new QDoubleSpinBox(&dialog); tolerance->setObjectName("ediPeriodTolerance"); tolerance->setDecimals(4); tolerance->setRange(0., 5.); tolerance->setSingleStep(.01); tolerance->setValue(.1);
  auto *toleranceLabel = ToolTips::label(tolerance, QObject::tr("Tolerance (%)"), QObject::tr("A group is merged only if (longest period - shortest period) / shortest period is strictly below this percentage. This limits the whole group, so chains of nearby values cannot bridge a larger gap."));
  controls->addWidget(toleranceLabel); controls->addWidget(tolerance);
  auto *gap = new QDoubleSpinBox(&dialog); gap->setObjectName("ediMaximumRatio"); gap->setRange(1.01, 1e6); gap->setValue(2.); gap->setKeyboardTracking(false);
  auto *gapLabel = ToolTips::label(gap, QObject::tr("Max gap ratio"), QObject::tr("Largest allowed ratio of upper to lower source period. Tests the whole bounding interval, even for targets near an endpoint. Outside each component's supported range, values remain unavailable."));
  controls->addWidget(gapLabel); controls->addWidget(gap);
  auto *barriers = new QCheckBox(QObject::tr("Respect masks"), &dialog); barriers->setObjectName("ediMaskBarriers"); barriers->setChecked(true);
  ToolTips::set(barriers, QObject::tr("Do not interpolate across masked samples. Exact masked samples remain masked even when this is off. Tipper and impedance masks are handled independently."));
  controls->addWidget(barriers); controls->addStretch(); layout->addLayout(controls);
  auto *periodStatus = new QLabel(&dialog); periodStatus->setObjectName("ediPeriodStatus"); periodStatus->setWordWrap(true); layout->addWidget(periodStatus);
  auto *table = new QTableWidget(&dialog); table->setObjectName("ediPeriodGroups"); table->setColumnCount(6);
  table->setEditTriggers(QAbstractItemView::NoEditTriggers); table->setSelectionBehavior(QAbstractItemView::SelectRows); table->verticalHeader()->hide();
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents); table->horizontalHeader()->setStretchLastSection(true); layout->addWidget(table);
  if(!duplicates.empty()) {
    QStringList names; for(const auto &name: duplicates) names.push_back(QString::fromStdString(name));
    auto *skipped = new QLabel(QObject::tr("%1 duplicate stations skipped").arg(duplicates.size()), &dialog);
    ToolTips::set(skipped, names.join('\n')); layout->addWidget(skipped);
  }
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  buttons->button(QDialogButtonBox::Ok)->setText(QObject::tr("Import")); buttons->button(QDialogButtonBox::Ok)->setEnabled(!incoming.empty()); layout->addWidget(buttons);
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept); QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  std::set<double> existingPeriods;
  for(const auto &name: candidate.get_stations_names()) if(!incoming.count(name))
    for(double f: candidate.get_station_data(name).frequencies())
      if(std::isfinite(f) && f > 0. && std::isfinite(1. / f) && 1. / f > 0.) existingPeriods.insert(1. / f);
  const std::vector<double> targets(existingPeriods.begin(), existingPeriods.end());
  // An interpolation preview is staged independently and applied only on acceptance.
  std::shared_ptr<MTSurveyData> interpolated;
  EDIPeriodMerge::Plan plan;
  const auto cells = [table](int row, const QStringList &values) {
    for(int c = 0; c < values.size(); ++c) table->setItem(row, c, new QTableWidgetItem(values[c]));
  };
  auto refresh = [&] {
    const bool interpolate = mode->currentIndex() == 2, merge = mode->currentIndex() == 1;
    tolerance->setVisible(!interpolate); toleranceLabel->setVisible(!interpolate);
    tolerance->setEnabled(merge);
    gap->setVisible(interpolate); gapLabel->setVisible(interpolate); barriers->setVisible(interpolate);
    interpolated.reset(); buttons->button(QDialogButtonBox::Ok)->setEnabled(!incoming.empty());
    table->clearContents();
    if(interpolate) {
      table->setColumnCount(5);
      table->setHorizontalHeaderLabels({QObject::tr("Period (s)"), QObject::tr("Existing"), QObject::tr("Interpolated"), QObject::tr("Derived"), QObject::tr("Unavailable")});
      ToolTips::set(table, QObject::tr("Counts include only new stations. Each value is one station/component: 4 complex impedances, 2 complex tippers and 4 real phase-tensor entries. Hover over Unavailable for causes. The common layout does not imply usable data at every target. Original incoming samples outside the target layout are omitted; the complete pre-interpolation survey and audit are retained in the project."));
      summary->setText(QObject::tr("%1 new / %2 total stations · %3 existing target periods").arg(incoming.size()).arg(candidate.n_stations()).arg(targets.size()));
      table->setRowCount(0); buttons->button(QDialogButtonBox::Ok)->setEnabled(false);
      if(targets.empty() || incoming.empty()) {
        periodStatus->setText(targets.empty() ? QObject::tr("No existing period layout. Choose Keep original or Merge close periods.") : QObject::tr("No new stations to import."));
        return;
      }
      try {
        interpolated = PeriodResampling::resample(candidate, targets, {gap->value(), barriers->isChecked()}, incoming);
        std::map<double, std::array<unsigned, 8>> counts;
        for(const auto &record: interpolated->resampling_info().records) ++counts[record.period][unsigned(record.status)];
        table->setRowCount(counts.size());
        unsigned retained = 0, estimated = 0, unavailable = 0; int row = 0;
        for(const auto &entry: counts) {
          const auto &n = entry.second;
          const unsigned missing = n[3] + n[4] + n[5] + n[6] + n[7];
          cells(row, {QString::number(entry.first, 'g', 12), QString::number(n[0]), QString::number(n[1]), QString::number(n[2]), QString::number(missing)});
          QStringList causes;
          for(unsigned s = 3; s < n.size(); ++s) if(n[s]) causes << QString::fromUtf8(PeriodResampling::status_name(PeriodResampling::Status(s))) + ": " + QString::number(n[s]);
          table->item(row, 4)->setToolTip(causes.join('\n'));
          retained += n[0]; estimated += n[1] + n[2]; unavailable += missing; ++row;
        }
        const bool usable = retained + estimated > 0;
        periodStatus->setText(QObject::tr("%1 retained · %2 estimated / derived · %3 unavailable values\n%4")
          .arg(retained).arg(estimated).arg(unavailable).arg(usable ? QObject::tr("New stations only; no extrapolation. Source and audit retained.") : QObject::tr("No usable target values. Adjust the gap limit or choose another import mode.")));
        buttons->button(QDialogButtonBox::Ok)->setEnabled(usable);
      } catch(const std::exception &error) { interpolated.reset(); periodStatus->setText(QString::fromUtf8(error.what())); }
      return;
    }
    table->setColumnCount(6);
    table->setHorizontalHeaderLabels({QObject::tr("Min (s)"), QObject::tr("Max (s)"), QObject::tr("Target (s)"), QObject::tr("Periods"), QObject::tr("Stations"), QObject::tr("Action")});
    ToolTips::set(table, QObject::tr("Period groups involving incoming stations. Merging assigns the lower median distinct period, or an existing survey period, without averaging or interpolation. Transfer functions, errors and masks are preserved; apparent resistivity is adjusted to the new frequency. Conflicting existing periods and same-station collisions stay unchanged. Keep original preserves all incoming frequencies."));
    plan = EDIPeriodMerge::analyze(candidate, incoming, tolerance->value() / 100.);
    summary->setText(QObject::tr("%1 new / %2 total stations · %3 stored samples\n%4 distinct periods → %5 after import · %6 mergeable / %7 blocked groups\nRange: %8 – %9 s · %10 samples %11")
      .arg(incoming.size()).arg(candidate.n_stations()).arg(plan.samples).arg(plan.before).arg(merge ? plan.after : plan.before)
      .arg(plan.mergeable).arg(plan.blocked).arg(plan.samples ? QString::number(plan.minimum, 'g', 9) : QObject::tr("N/A"))
      .arg(plan.samples ? QString::number(plan.maximum, 'g', 9) : QObject::tr("N/A")).arg(plan.changes.size())
      .arg(merge ? QObject::tr("will move") : QObject::tr("can move")) +
      (plan.invalid ? QObject::tr("\n%1 invalid periods excluded from statistics; kept unchanged.").arg(plan.invalid) : QString()));
    periodStatus->setText(incoming.empty() ? QObject::tr("No new stations to import.") : plan.groups.empty() ? QObject::tr("No valid periods to display.") :
      !plan.mergeable && !plan.blocked ? QObject::tr("No close periods to merge. Click Import to keep the original periods.") : QObject::tr("Period layout"));
    table->setRowCount(plan.groups.size());
    for(unsigned i = 0; i < plan.groups.size(); ++i) {
      const auto &g = plan.groups[i]; const bool blocked = g.collision || g.existingConflict;
      const QString action = g.collision ? QObject::tr("Keep: same-station samples") : g.existingConflict ? QObject::tr("Keep: existing periods") :
        g.periods == 1 ? (g.stations > 1 ? QObject::tr("Keep (already shared)") : QObject::tr("Keep")) :
        merge ? QObject::tr("Merge") : QObject::tr("Keep (merge available)");
      cells(i, {QString::number(g.minimum, 'g', 12), QString::number(g.maximum, 'g', 12), blocked ? QStringLiteral("—") : QString::number(g.target, 'g', 12), QString::number(g.periods), QString::number(g.stations), action});
    }
  };
  QObject::connect(tolerance, QOverload<double>::of(&QDoubleSpinBox::valueChanged), &dialog, refresh);
  QObject::connect(mode, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, refresh);
  QObject::connect(gap, QOverload<double>::of(&QDoubleSpinBox::valueChanged), &dialog, refresh);
  QObject::connect(barriers, &QCheckBox::toggled, &dialog, refresh); refresh(); dialog.resize(850, 460);
  if(dialog.exec() != QDialog::Accepted) return false;
  if(mode->currentIndex() == 1) EDIPeriodMerge::apply(candidate, plan);
  if(mode->currentIndex() == 2) {
    if(!interpolated || !buttons->button(QDialogButtonBox::Ok)->isEnabled()) return false;
    candidate.restore(*interpolated);
  }
  return true;
}
