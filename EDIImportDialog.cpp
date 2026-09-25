#include "EDIImportDialog.h"
#include "include/EDIPeriodMerge.h"
#include "include/HelpButton.h"
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

bool reviewEDIImport(QWidget *parent, MTSurveyData &candidate,
                     const std::set<std::string> &incoming, const std::vector<std::string> &duplicates)
{
  QDialog dialog(parent); dialog.setObjectName("ediImportDialog"); dialog.setWindowTitle(QObject::tr("Import EDIs"));
  auto *layout = new QVBoxLayout(&dialog);
  auto *summary = new QLabel(&dialog); summary->setObjectName("ediImportSummary"); summary->setWordWrap(true); layout->addWidget(summary);
  auto *controls = new QHBoxLayout;
  auto *merge = new QCheckBox(QObject::tr("Merge close periods"), &dialog); merge->setObjectName("ediMergePeriods"); controls->addWidget(merge);
  controls->addWidget(UiHelp::label(&dialog, QObject::tr("Tolerance (%)"), QObject::tr("A group is merged only if (longest period - shortest period) / shortest period is strictly below this percentage. This limits the whole group, so chains of nearby values cannot bridge a larger gap.")));
  auto *tolerance = new QDoubleSpinBox(&dialog); tolerance->setObjectName("ediPeriodTolerance"); tolerance->setDecimals(4); tolerance->setRange(0., 5.); tolerance->setSingleStep(.01); tolerance->setValue(.1); controls->addWidget(tolerance); controls->addStretch(); layout->addLayout(controls);
  layout->addWidget(UiHelp::label(&dialog, QObject::tr("Period groups"), QObject::tr("Statistics include existing survey stations and new stations. Each row shows a nearby-period group involving new data. Merging gives new stations a shared stored frequency, using the lower median distinct period or an existing survey period. It does not average measurements or interpolate. Z, tipper, phase tensor, errors and masks are preserved; apparent resistivity and its error are adjusted to the assigned frequency. Groups with multiple samples from one station or conflicting existing periods stay unchanged. Existing stations and EDI files are never modified. Leave Merge close periods off to keep all original frequencies."), "ediImportHelp"));
  auto *table = new QTableWidget(&dialog); table->setObjectName("ediPeriodGroups"); table->setColumnCount(6);
  table->setHorizontalHeaderLabels({QObject::tr("Min (s)"), QObject::tr("Max (s)"), QObject::tr("Target (s)"), QObject::tr("Periods"), QObject::tr("Stations"), QObject::tr("Action")});
  table->setEditTriggers(QAbstractItemView::NoEditTriggers); table->setSelectionBehavior(QAbstractItemView::SelectRows); table->verticalHeader()->hide();
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents); table->horizontalHeader()->setStretchLastSection(true); layout->addWidget(table);
  if(!duplicates.empty()) {
    QStringList names; for(const auto &name: duplicates) names.push_back(QString::fromStdString(name));
    layout->addWidget(UiHelp::label(&dialog, QObject::tr("%1 duplicate stations skipped").arg(duplicates.size()), names.join('\n')));
  }
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  buttons->button(QDialogButtonBox::Ok)->setText(QObject::tr("Import")); buttons->button(QDialogButtonBox::Ok)->setEnabled(!incoming.empty()); layout->addWidget(buttons);
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept); QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  EDIPeriodMerge::Plan plan;
  auto refresh = [&] {
    plan = EDIPeriodMerge::analyze(candidate, incoming, tolerance->value() / 100.);
    summary->setText(QObject::tr("%1 new / %2 total stations · %3 stored samples\n%4 distinct periods → %5 after import · %6 mergeable / %7 unchanged groups\nRange: %8 – %9 s · %10 samples %11")
      .arg(incoming.size()).arg(candidate.n_stations()).arg(plan.samples).arg(plan.before).arg(merge->isChecked() ? plan.after : plan.before)
      .arg(plan.mergeable).arg(plan.blocked).arg(plan.samples ? QString::number(plan.minimum, 'g', 9) : QObject::tr("N/A"))
      .arg(plan.samples ? QString::number(plan.maximum, 'g', 9) : QObject::tr("N/A")).arg(plan.changes.size())
      .arg(merge->isChecked() ? QObject::tr("will move") : QObject::tr("can move")) +
      (plan.invalid ? QObject::tr("\n%1 invalid periods excluded from statistics; kept unchanged.").arg(plan.invalid) : QString()));
    table->setRowCount(plan.groups.size());
    for(unsigned i = 0; i < plan.groups.size(); ++i) {
      const auto &g = plan.groups[i]; const bool blocked = g.collision || g.existingConflict;
      const QString action = g.collision ? QObject::tr("Keep: same-station samples") : g.existingConflict ? QObject::tr("Keep: existing periods") :
        merge->isChecked() ? QObject::tr("Merge") : QObject::tr("Keep (merge available)");
      const QStringList cells{QString::number(g.minimum, 'g', 12), QString::number(g.maximum, 'g', 12), blocked ? QStringLiteral("—") : QString::number(g.target, 'g', 12), QString::number(g.periods), QString::number(g.stations), action};
      for(int c = 0; c < cells.size(); ++c) table->setItem(i, c, new QTableWidgetItem(cells[c]));
    }
  };
  QObject::connect(tolerance, QOverload<double>::of(&QDoubleSpinBox::valueChanged), &dialog, refresh);
  QObject::connect(merge, &QCheckBox::toggled, &dialog, refresh); refresh(); dialog.resize(850, 430);
  if(dialog.exec() != QDialog::Accepted) return false;
  if(merge->isChecked()) EDIPeriodMerge::apply(candidate, plan);
  return true;
}
