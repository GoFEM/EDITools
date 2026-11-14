#include "ui/dialogs/ExportNativeMTDialog.h"
#include "ui/dialogs/ExportPreview.h"
#include <QTimer>
#include "core/MTSurveyData.h"

#include <QCheckBox>
#include <QDialogButtonBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QToolButton>
#include <QVBoxLayout>

#include <cmath>

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

  preview = new ExportPreview(this); layout->addWidget(preview);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel, this);
  save=buttons->button(QDialogButtonBox::Save); save->setText(tr("Export…"));
  auto *timer=new QTimer(this); timer->setSingleShot(true);
  connect(timer,&QTimer::timeout,this,&ExportNativeMTDialog::updatePreview);
  auto changed=[timer] { timer->start(0); };
  connect(preview,&ExportPreview::optionsChanged,this,changed);
  connect(periods,&QListWidget::itemSelectionChanged,this,changed);
  connect(types,&QListWidget::itemChanged,this,changed);
  for(auto *box:{conjugate,reverse_vertical,distortion})connect(box,&QCheckBox::toggled,this,changed);
  connect(buttons, &QDialogButtonBox::accepted, this, &ExportNativeMTDialog::exportFiles);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  layout->addWidget(buttons);
  updatePreview();
  resize(780,800);
}

NativeMT::Options ExportNativeMTDialog::options() const {
  NativeMT::Options result;
  for(auto *item:periods->selectedItems())result.periods.push_back(item->data(Qt::UserRole).toDouble());
  for(int i=0;i<types->count();++i)if(types->item(i)->checkState()==Qt::Checked)
    result.types.push_back(static_cast<RealDataType>(types->item(i)->data(Qt::UserRole).toInt()));
  result.negative_time_convention=conjugate->isChecked();
  result.vertical_magnetic_up=reverse_vertical->isChecked();result.distortion=distortion->isChecked();
  return result;
}
void ExportNativeMTDialog::updatePreview() {
  save->setEnabled(preview->update(survey,options(),SurveyExport::Format::Native));
}
void ExportNativeMTDialog::exportFiles() {
  try {
    updatePreview();if(!save->isEnabled())return;
    exported_file=saveSurveyExport(this,preview->result(),SurveyExport::Format::Native,directory);
    if(exported_file.isEmpty())return;
    QMessageBox::information(this,tr("MT export complete"),tr("Observations and receiver/frequency companions saved."));
    accept();
  } catch(const std::exception &e) { QMessageBox::warning(this,tr("Cannot export MT data"),QString::fromUtf8(e.what())); }
}
