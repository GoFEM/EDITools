#include "ExportNativeMTDialog.h"
#include "include/MTSurveyData.h"

#include <QCheckBox>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QSaveFile>
#include <QToolButton>
#include <QVBoxLayout>

#include <cmath>
#include <memory>

ExportNativeMTDialog::ExportNativeMTDialog(MTSurveyData &data, const QString &path, QWidget *parent)
  : QDialog(parent), survey(data), directory(path)
{
  setWindowTitle(tr("Export to MT Inversion"));
  resize(630, 590);
  auto *layout = new QVBoxLayout(this);
  auto *selection = new QWidget(this);
  layout->addWidget(selection);
  auto *selection_layout = new QHBoxLayout(selection);
  auto *period_layout = new QVBoxLayout;
  period_layout->addWidget(new QLabel(tr("Periods to export (seconds)")));
  periods = new QListWidget(selection);
  periods->setObjectName("nativePeriods");
  periods->setSelectionMode(QAbstractItemView::MultiSelection);
  // The legacy period list merges with an absolute tolerance. Keep exact
  // per-station values here so native export can diagnose frequency collisions.
  std::set<double> period_values;
  for(const auto &name: survey.get_stations_names()) {
    const auto &station = survey.get_station_data(name);
    if(!station.active()) continue;
    for(double f: station.frequencies())
      if(std::isfinite(f) && f > 0 && std::isfinite(1. / f)) period_values.insert(1. / f);
  }
  for(double period: period_values) {
    auto *item = new QListWidgetItem(QString::number(period, 'g', 17), periods);
    item->setData(Qt::UserRole, period);
    item->setSelected(true);
  }
  period_layout->addWidget(periods);
  auto *all_periods = new QPushButton(tr("Select all periods"));
  connect(all_periods, &QPushButton::clicked, periods, &QListWidget::selectAll);
  period_layout->addWidget(all_periods);
  selection_layout->addLayout(period_layout);
  auto *type_layout = new QVBoxLayout;
  type_layout->addWidget(new QLabel(tr("Scalar observations")));
  types = new QListWidget(selection);
  types->setObjectName("nativeTypes");
  for(const auto &m: NativeMT::mappings()) {
    auto *item = new QListWidgetItem(QString("%1  %2").arg(m.observable, m.component), types);
    item->setData(Qt::UserRole, int(m.type));
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(int(m.type) >= int(RealZxx) && int(m.type) <= int(ImagZyy) ? Qt::Checked : Qt::Unchecked);
  }
  type_layout->addWidget(types);
  auto *clear_types = new QPushButton(tr("Clear components"));
  connect(clear_types, &QPushButton::clicked, this, [this] {
    for(int i = 0; i < types->count(); ++i) types->item(i)->setCheckState(Qt::Unchecked);
  });
  type_layout->addWidget(clear_types);
  selection_layout->addLayout(type_layout);
  auto *coordinate_help = new QLabel(tr("Uses the survey coordinates. Convert or center them via Coordinates → Convert latitude/longitude to UTM. If needed, UTM is chosen automatically."), this);
  coordinate_help->setWordWrap(true);
  layout->addWidget(coordinate_help);
  auto *advanced_button = new QToolButton(this);
  advanced_button->setText(tr("Advanced"));
  advanced_button->setCheckable(true);
  advanced_button->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  advanced_button->setArrowType(Qt::RightArrow);
  layout->addWidget(advanced_button);
  auto *advanced = new QWidget(this);
  auto *advanced_layout = new QVBoxLayout(advanced);
  conjugate = new QCheckBox(tr("Conjugate responses (source uses exp(−iωt))"), advanced);
  conjugate->setObjectName("nativeConjugate");
  reverse_vertical = new QCheckBox(tr("Reverse tipper sign (source magnetic z is upward)"), advanced);
  distortion = new QCheckBox(tr("Check distortion-estimation compatibility"), advanced);
  advanced_layout->addWidget(conjugate);
  advanced_layout->addWidget(reverse_vertical);
  advanced_layout->addWidget(distortion);
  layout->addWidget(advanced);
  advanced->hide();
  connect(advanced_button, &QToolButton::toggled, this, [advanced, advanced_button](bool on) {
    advanced->setVisible(on);
    advanced_button->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
  });

  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel, this);
  buttons->button(QDialogButtonBox::Save)->setText(tr("Export…"));
  connect(buttons, &QDialogButtonBox::accepted, this, &ExportNativeMTDialog::exportFiles);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  layout->addWidget(buttons);
}

void ExportNativeMTDialog::exportFiles()
{
  try {
    NativeMT::Options options;
    for(auto *item: periods->selectedItems()) options.periods.push_back(item->data(Qt::UserRole).toDouble());
    for(int i = 0; i < types->count(); ++i)
      if(types->item(i)->checkState() == Qt::Checked)
        options.types.push_back(static_cast<RealDataType>(types->item(i)->data(Qt::UserRole).toInt()));
    options.horizontal_axes_confirmed = true;
    options.conventions_confirmed = true;
    options.negative_time_convention = conjugate->isChecked();
    options.vertical_magnetic_up = reverse_vertical->isChecked();
    options.distortion = distortion->isChecked();
    survey.ensure_utm_coordinates();
    const auto names = survey.get_stations_names();
    const auto positions = survey.get_stations_locations();
    for(unsigned i = 0; i < names.size(); ++i)
      options.receivers[names[i]] = {{positions[i][0], positions[i][1], -positions[i][2]}};
    options.coordinate_description = survey.coordinates().description() +
      "\nNative z = -EDI elevation (metres); vertical origin elevation = 0, source elevation datum retained.";
    std::vector<const MTStationData *> stations;
    for(const auto &name: survey.get_stations_names()) stations.push_back(&survey.get_station_data(name));
    const auto result = NativeMT::prepare(stations, options);

    QString summary = tr("%1 observations, %2 receivers, %3 frequencies.\nOmitted: %4 masked, %5 missing, %6 undefined phases at zero impedance.")
      .arg(qulonglong(result.observations.size())).arg(qulonglong(result.receivers.size())).arg(qulonglong(result.frequencies.size()))
      .arg(qulonglong(result.masked)).arg(qulonglong(result.missing)).arg(qulonglong(result.undefined_phase));
    for(const auto &warning: result.warnings) summary += "\n\n" + QString::fromStdString(warning);
    QString path = QFileDialog::getSaveFileName(this, tr("Export native MT observations"), directory, tr("MT observations (*.data);;All files (*)"));
    if(path.isEmpty()) return;
    if(QFileInfo(path).suffix().isEmpty()) path += ".data";
    QFileInfo info(path);
    QString base = info.absolutePath() + "/" + info.completeBaseName();
    const std::array<QString, 3> paths{{path, base + ".recvs", base + ".freqs"}};
    if(paths[0] == paths[1] || paths[0] == paths[2])
      throw std::runtime_error("Choose an observation filename different from the .recvs and .freqs companions.");
    QStringList existing;
    for(unsigned i = 1; i < paths.size(); ++i) if(QFileInfo::exists(paths[i])) existing << paths[i];
    if(!existing.empty() && QMessageBox::question(this, tr("Replace companion files?"), existing.join("\n")) != QMessageBox::Yes) return;
    const std::array<std::string, 3> texts{{result.data_text, result.receiver_text, result.frequency_text}};
    std::array<std::unique_ptr<QSaveFile>, 3> files;
    // Stage all three files before replacing any destination.
    for(unsigned i = 0; i < paths.size(); ++i) {
      files[i].reset(new QSaveFile(paths[i]));
      if(!files[i]->open(QIODevice::WriteOnly) ||
         files[i]->write(texts[i].data(), qint64(texts[i].size())) != qint64(texts[i].size()))
        throw std::runtime_error("Cannot write " + paths[i].toStdString() + ": " + files[i]->errorString().toStdString());
    }
    for(unsigned i = 0; i < paths.size(); ++i)
      if(!files[i]->commit()) throw std::runtime_error("Cannot finalize " + paths[i].toStdString() +
        ". Some companion files may already have been replaced; retry the complete export. " + files[i]->errorString().toStdString());
    exported_file = path;
    QMessageBox::information(this, tr("MT export complete"), summary + "\n\n" + paths[0] + "\n" + paths[1] + "\n" + paths[2]);
    accept();
  } catch(const std::exception &e) {
    QMessageBox::warning(this, tr("Cannot export MT data"), QString::fromUtf8(e.what()));
  }
}
