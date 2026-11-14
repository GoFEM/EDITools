#include "ui/masking/GroupMaskDialog.h"
#include "ui/masking/GroupMaskPlots.h"
#include "core/MTSurveyData.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>
#include <cmath>

bool editGroupMasks(QWidget *parent, MTSurveyData &survey,
                    const std::set<std::string> &stations, bool linkTensorMasks)
{
  QDialog dialog(parent);
  dialog.setObjectName("groupMaskDialog"); dialog.setWindowTitle(QObject::tr("Group masking"));
  auto *outer = new QVBoxLayout(&dialog);
  auto *splitter = new QSplitter(&dialog); outer->addWidget(splitter, 1);
  auto *editor = new QWidget(splitter);
  auto *layout = new QVBoxLayout(editor); layout->setContentsMargins(0, 0, 0, 0);
  auto *scope = new QLabel(QObject::tr("%1 selected stations").arg(stations.size()), &dialog);
  scope->setToolTip(QObject::tr("Applies to the selected stations, including disabled stations. Select sites in the main station list before opening this dialog. Plot visibility does not change the component choices here."));
  layout->addWidget(scope);
  auto *controls = new QHBoxLayout;
  auto *operation = new QComboBox(&dialog); operation->setObjectName("groupMaskOperation");
  operation->addItems({QObject::tr("Mask"), QObject::tr("Unmask")}); controls->addWidget(operation);
  auto *mode = new QComboBox(&dialog); mode->setObjectName("groupMaskPeriods");
  mode->addItems({QObject::tr("All periods"), QObject::tr("Inside range"), QObject::tr("Outside range"), QObject::tr("Selected periods")});
  mode->setToolTip(QObject::tr("Ranges include their endpoints. Selected periods match exact stored periods; nearby values are not merged. Only samples already stored at each station are edited."));
  controls->addWidget(mode);
  auto *lower = new QLineEdit(&dialog), *upper = new QLineEdit(&dialog);
  lower->setObjectName("groupMaskLower"); upper->setObjectName("groupMaskUpper");
  lower->setPlaceholderText(QObject::tr("Min (s)")); upper->setPlaceholderText(QObject::tr("Max (s)"));
  lower->setToolTip(QObject::tr("Minimum period in seconds, inclusive. Scientific notation is accepted."));
  upper->setToolTip(QObject::tr("Maximum period in seconds, inclusive. Scientific notation is accepted."));
  controls->addWidget(lower); controls->addWidget(upper); layout->addLayout(controls);
  auto *periods = new QListWidget(&dialog); periods->setObjectName("groupMaskPeriodList");
  periods->setSelectionMode(QAbstractItemView::ExtendedSelection); periods->setMaximumHeight(160);
  periods->setToolTip(QObject::tr("Ctrl-click or Shift-click to select periods. Ctrl+A selects every listed period. Values are the exact union at the selected stations."));
  std::set<double> values;
  for(const auto &name : stations) if(survey.is_station_present(name))
    for(double f : survey.get_station_data(name).frequencies())
      if(std::isfinite(f) && f > 0. && std::isfinite(1. / f) && 1. / f > 0.) values.insert(1. / f);
  for(double p : values) {
    auto *item = new QListWidgetItem(QString::number(p, 'g', 17), periods);
    item->setData(Qt::UserRole, p);
  }
  if(!values.empty()) {
    lower->setText(QString::number(*values.begin(), 'g', 17));
    upper->setText(QString::number(*values.rbegin(), 'g', 17));
  }
  layout->addWidget(periods);
  const std::array<RealDataType, 10> types{{RealZxx, RealZxy, RealZyx, RealZyy, RealTzx, RealTzy, PTxx, PTxy, PTyx, PTyy}};
  const QStringList labels{"Zxx", "Zxy", "Zyx", "Zyy", "Tzx", "Tzy", "PTxx", "PTxy", "PTyx", "PTyy"};
  std::array<QCheckBox *, 10> components;
  auto *groups = new QHBoxLayout;
  for(unsigned first : {0u, 4u, 6u}) {
    auto *group = new QGroupBox(first == 0 ? QObject::tr("Z / ρ / phase") : first == 4 ? QObject::tr("Tipper") : QObject::tr("PT"), &dialog);
    auto *column = new QVBoxLayout(group);
    const unsigned end = first == 0 ? 4 : first == 4 ? 6 : 10;
    for(unsigned c = first; c < end; ++c) {
      auto *box = components[c] = new QCheckBox(labels[c], group);
      box->setObjectName("groupMask" + labels[c]); column->addWidget(box);
      box->setToolTip(QObject::tr("Edit this component at the chosen periods. Z shares its mask with apparent resistivity and phase; complex components share real/imaginary masks. Zxy also edits Zxx; Zyx also edits Zyy, independently of Link Z / PT."));
    }
    column->addStretch(); groups->addWidget(group);
  }
  layout->addLayout(groups);
  auto *selection = new QHBoxLayout;
  auto *all = new QPushButton(QObject::tr("All components"), &dialog);
  auto *none = new QPushButton(QObject::tr("Clear"), &dialog);
  all->setObjectName("groupMaskAll"); none->setObjectName("groupMaskClear");
  selection->addWidget(all); selection->addWidget(none); selection->addStretch();
  auto *linked = new QCheckBox(QObject::tr("Link Z / PT"), &dialog); linked->setObjectName("groupMaskLink");
  linked->setChecked(linkTensorMasks);
  linked->setToolTip(QObject::tr("Starts with the main Z/PT setting; changes here affect this operation only. A selected Z component also edits all PT entries; a selected PT entry also edits all Z components. Links do not propagate recursively. Tippers remain independent."));
  selection->addWidget(linked); layout->addLayout(selection);
  auto *table = new QTableWidget(&dialog); table->setObjectName("groupMaskPreview");
  table->setColumnCount(5); table->setHorizontalHeaderLabels({QObject::tr("Station"), QObject::tr("Periods"), "Z", "Tipper", "PT"});
  table->setEditTriggers(QAbstractItemView::NoEditTriggers); table->verticalHeader()->hide();
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  table->setToolTip(QObject::tr("Periods count stored samples matching the condition. Z, Tipper and PT count mask flags that will change, including linked components, once each. Flags also exist for missing values; unmasking never creates data. Already matching flags are excluded."));
  layout->addWidget(table);
  auto *summary = new QLabel(&dialog); summary->setObjectName("groupMaskSummary"); summary->setWordWrap(true); layout->addWidget(summary);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
  auto *apply = buttons->button(QDialogButtonBox::Ok); apply->setText(QObject::tr("Apply")); apply->setObjectName("groupMaskApply");
  outer->addWidget(buttons);
  auto *plots = new GroupMaskPlots(splitter, survey, stations);
  splitter->setStretchFactor(1, 1); splitter->setSizes({420, 900});
  GroupMaskEdits edits;
  const auto refresh = [&] {
    edits.clear(); table->setRowCount(0); apply->setEnabled(false);
    const int condition = mode->currentIndex();
    for(auto *edit : {lower, upper}) edit->setVisible(condition == 1 || condition == 2);
    periods->setVisible(condition == 3);
    bool validLow, validHigh;
    const double low = lower->text().toDouble(&validLow), high = upper->text().toDouble(&validHigh);
    if((condition == 1 || condition == 2) && (!validLow || !validHigh || !std::isfinite(low) || !std::isfinite(high) || low <= 0. || high < low)) {
      plots->setEdits({}, operation->currentIndex() == 1);
      summary->setText(QObject::tr("Enter a positive period range: minimum ≤ maximum.")); return;
    }
    std::set<double> chosen;
    for(auto *item : periods->selectedItems()) chosen.insert(item->data(Qt::UserRole).toDouble());
    std::set<RealDataType> selected;
    for(unsigned c = 0; c < types.size(); ++c) if(components[c]->isChecked()) {
      const auto targets = MTStationData::mask_edit_types(types[c], linked->isChecked());
      selected.insert(targets.begin(), targets.end());
    }
    size_t total = 0;
    const bool enabled = operation->currentIndex() == 1;
    for(const auto &name : stations) if(survey.is_station_present(name)) {
      const auto &station = survey.get_station_data(name);
      const std::array<std::vector<std::vector<bool>>, 3> masks{{station.impedance_mask(true), station.tipper_mask(true), station.phase_tensor_mask(true)}};
      std::array<unsigned, 3> counts{}; unsigned matching = 0;
      std::vector<std::pair<RealDataType, double>> changes;
      for(unsigned i = 0; i < station.frequencies().size(); ++i) {
        const double f = station.frequencies()[i], p = 1. / f;
        if(!std::isfinite(f) || f <= 0. || !std::isfinite(p) || p <= 0.) continue;
        const bool inside = p >= low && p <= high;
        if((condition == 1 && !inside) || (condition == 2 && inside) || (condition == 3 && !chosen.count(p))) continue;
        ++matching;
        for(unsigned c = 0; c < types.size(); ++c) if(selected.count(types[c])) {
          const unsigned family = c < 4 ? 0 : c < 6 ? 1 : 2, column = c < 4 ? c : c < 6 ? c - 4 : c - 6;
          if(masks[family][column][i] != enabled) { changes.emplace_back(types[c], f); ++counts[family]; }
        }
      }
      total += changes.size();
      if(!changes.empty()) edits.emplace(name, std::move(changes));
      const int row = table->rowCount(); table->insertRow(row);
      const QStringList cells{QString::fromStdString(name), QString::number(matching), QString::number(counts[0]), QString::number(counts[1]), QString::number(counts[2])};
      for(int c = 0; c < cells.size(); ++c) table->setItem(row, c, new QTableWidgetItem(cells[c]));
    }
    summary->setText(selected.empty() ? QObject::tr("Choose components to edit.") :
      QObject::tr("%1 stations affected · %2 mask changes").arg(edits.size()).arg(total));
    apply->setEnabled(total > 0);
    plots->setEdits(edits, enabled);
  };
  for(auto *box : components) QObject::connect(box, &QCheckBox::toggled, &dialog, refresh);
  QObject::connect(linked, &QCheckBox::toggled, &dialog, refresh);
  for(auto *combo : {operation, mode}) QObject::connect(combo, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, refresh);
  for(auto *edit : {lower, upper}) QObject::connect(edit, &QLineEdit::textChanged, &dialog, refresh);
  QObject::connect(periods, &QListWidget::itemSelectionChanged, &dialog, refresh);
  for(auto *button : {all, none}) QObject::connect(button, &QPushButton::clicked, &dialog, [&, button] {
    for(auto *box : components) { const QSignalBlocker block(box); box->setChecked(button == all); } refresh();
  });
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  refresh(); dialog.resize(1320, 780);
  if(dialog.exec() != QDialog::Accepted) return false;
  refresh(); // Apply precisely the preview; never expand linked entries recursively.
  for(const auto &station : edits) for(const auto &edit : station.second)
    survey.get_station_data(station.first).set_data_mask(edit.first, edit.second, operation->currentIndex() == 1);
  return !edits.empty();
}
