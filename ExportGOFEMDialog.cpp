#include "ExportGOFEMDialog.h"
#include "ui_ExportGOFEMDialog.h"

ExportGOFEMDialog::ExportGOFEMDialog(std::vector<double> periods, QWidget *parent) :
  QDialog(parent),
  ui(new Ui::ExportGOFEMDialog),
  m_periods(periods)
{
  ui->setupUi(this);

  for(double period: periods)
    ui->periodListWidget->addItem(QString("%1").arg(period));

//  for(int i = 0; i < ui->periodListWidget->count(); ++i)
//  {
//      QListWidgetItem* item = ui->periodListWidget->item(i);
//      item->setSelected(true);
//  }
}

ExportGOFEMDialog::~ExportGOFEMDialog()
{
  delete ui;
}

std::vector<RealDataType> ExportGOFEMDialog::getSelectedDataTypes() const
{
  std::vector<RealDataType> selected;

  if(ui->zxxCheckBox->isChecked())
    selected.push_back(RealZxx), selected.push_back(ImagZxx);
  if(ui->zxyCheckBox->isChecked())
    selected.push_back(RealZxy), selected.push_back(ImagZxy);
  if(ui->zyxCheckBox->isChecked())
    selected.push_back(RealZyx), selected.push_back(ImagZyx);
  if(ui->zyyCheckBox->isChecked())
    selected.push_back(RealZyy), selected.push_back(ImagZyy);

  if(ui->pxxCheckBox->isChecked())
    selected.push_back(PTxx);
  if(ui->pxyCheckBox->isChecked())
    selected.push_back(PTxy);
  if(ui->pyxCheckBox->isChecked())
    selected.push_back(PTyx);
  if(ui->pyyCheckBox->isChecked())
    selected.push_back(PTyy);

  if(ui->rxxCheckBox->isChecked())
    selected.push_back(RhoZxx);
  if(ui->rxyCheckBox->isChecked())
    selected.push_back(RhoZxy);
  if(ui->ryxCheckBox->isChecked())
    selected.push_back(RhoZyx);
  if(ui->ryyCheckBox->isChecked())
    selected.push_back(RhoZyy);

  if(ui->phixxCheckBox->isChecked())
    selected.push_back(PhsZxx);
  if(ui->phixyCheckBox->isChecked())
    selected.push_back(PhsZxy);
  if(ui->phiyxCheckBox->isChecked())
    selected.push_back(PhsZyx);
  if(ui->phiyyCheckBox->isChecked())
    selected.push_back(PhsZyy);

  if(ui->tzxCheckBox->isChecked())
    selected.push_back(RealTzx), selected.push_back(ImagTzx);
  if(ui->tzyCheckBox->isChecked())
    selected.push_back(RealTzy), selected.push_back(ImagTzy);

  return selected;
}

std::vector<double> ExportGOFEMDialog::getSelectedPeriods() const
{
  auto selected = ui->periodListWidget->selectionModel()->selectedIndexes();

  std::vector<double> selected_periods;
  for(auto index: selected)
    selected_periods.push_back(m_periods[index.row()]);

  return selected_periods;
}

void ExportGOFEMDialog::on_impGroupBox_toggled(bool on)
{
  if(on)
  {
    ui->zxxCheckBox->setChecked(true);
    ui->zxyCheckBox->setChecked(true);
    ui->zyxCheckBox->setChecked(true);
    ui->zyyCheckBox->setChecked(true);
  }
  else
  {
    ui->zxxCheckBox->setChecked(false);
    ui->zxyCheckBox->setChecked(false);
    ui->zyxCheckBox->setChecked(false);
    ui->zyyCheckBox->setChecked(false);
  }
}

void ExportGOFEMDialog::on_ptGroupBox_toggled(bool on)
{
  if(on)
  {
    ui->pxxCheckBox->setChecked(true);
    ui->pxyCheckBox->setChecked(true);
    ui->pyxCheckBox->setChecked(true);
    ui->pyyCheckBox->setChecked(true);
  }
  else
  {
    ui->pxxCheckBox->setChecked(false);
    ui->pxyCheckBox->setChecked(false);
    ui->pyxCheckBox->setChecked(false);
    ui->pyyCheckBox->setChecked(false);
  }
}

void ExportGOFEMDialog::on_tipperGroupBox_toggled(bool on)
{
  if(on)
  {
    ui->tzxCheckBox->setChecked(true);
    ui->tzyCheckBox->setChecked(true);
  }
  else
  {
    ui->tzxCheckBox->setChecked(false);
    ui->tzyCheckBox->setChecked(false);
  }
}

void ExportGOFEMDialog::on_periodListWidget_itemSelectionChanged()
{
  auto selected = ui->periodListWidget->selectionModel()->selectedIndexes();
  ui->statusLabel->setText(QString("Selected periods: %1").arg(selected.size()));
}

void ExportGOFEMDialog::on_rhoGroupBox_toggled(bool on)
{
  if(on)
  {
    ui->rxxCheckBox->setChecked(true);
    ui->rxyCheckBox->setChecked(true);
    ui->ryxCheckBox->setChecked(true);
    ui->ryyCheckBox->setChecked(true);
  }
  else
  {
    ui->rxxCheckBox->setChecked(false);
    ui->rxyCheckBox->setChecked(false);
    ui->ryxCheckBox->setChecked(false);
    ui->ryyCheckBox->setChecked(false);
  }
}

void ExportGOFEMDialog::on_phaseGroupBox_toggled(bool on)
{
  if(on)
  {
    ui->phixxCheckBox->setChecked(true);
    ui->phixyCheckBox->setChecked(true);
    ui->phiyxCheckBox->setChecked(true);
    ui->phiyyCheckBox->setChecked(true);
  }
  else
  {
    ui->phixxCheckBox->setChecked(false);
    ui->phixyCheckBox->setChecked(false);
    ui->phiyxCheckBox->setChecked(false);
    ui->phiyyCheckBox->setChecked(false);
  }
}
