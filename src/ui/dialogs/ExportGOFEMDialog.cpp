/*
 * The EDI Tools application.
 *
 * Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "ui/dialogs/ExportGOFEMDialog.h"
#include "ui_ExportGOFEMDialog.h"
#include "ui/dialogs/ExportPreview.h"
#include "core/MTSurveyData.h"
#include <QPushButton>
#include <QTimer>

ExportGOFEMDialog::ExportGOFEMDialog(const MTSurveyData &data, QWidget *parent) :
  QDialog(parent),
  ui(new Ui::ExportGOFEMDialog),
  survey(data),
  m_periods(data.get_unique_periods())
{
  ui->setupUi(this);

  for(double period: m_periods)
    ui->periodListWidget->addItem(QString("%1").arg(period));

  preview=new ExportPreview(this);
  ui->gridLayout_4->addWidget(preview,0,2,4,1);
  auto *timer=new QTimer(this);timer->setSingleShot(true);
  connect(timer,&QTimer::timeout,this,&ExportGOFEMDialog::updatePreview);
  auto changed=[timer] { timer->start(0); };
  connect(preview,&ExportPreview::optionsChanged,this,changed);
  connect(ui->periodListWidget,&QListWidget::itemSelectionChanged,this,changed);
  for(auto *box:findChildren<QCheckBox*>())
    if(!preview->isAncestorOf(box)) connect(box,&QCheckBox::toggled,this,changed);
  ui->periodListWidget->selectAll();updatePreview();resize(1020,610);
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

namespace {
void checkComponents(QGroupBox *group, bool on) {
  for(auto *box:group->findChildren<QCheckBox*>())box->setChecked(on);
}
}
void ExportGOFEMDialog::on_impGroupBox_toggled(bool on) { checkComponents(ui->impGroupBox,on); }
void ExportGOFEMDialog::on_ptGroupBox_toggled(bool on) { checkComponents(ui->ptGroupBox,on); }
void ExportGOFEMDialog::on_tipperGroupBox_toggled(bool on) { checkComponents(ui->tipperGroupBox,on); }
void ExportGOFEMDialog::on_rhoGroupBox_toggled(bool on) { checkComponents(ui->rhoGroupBox,on); }
void ExportGOFEMDialog::on_phaseGroupBox_toggled(bool on) { checkComponents(ui->phaseGroupBox,on); }

void ExportGOFEMDialog::on_periodListWidget_itemSelectionChanged()
{
  auto selected = ui->periodListWidget->selectionModel()->selectedIndexes();
  ui->statusLabel->setText(QString("Selected periods: %1").arg(selected.size()));
}

void ExportGOFEMDialog::updatePreview() {
  NativeMT::Options options;options.types=getSelectedDataTypes();options.periods=getSelectedPeriods();
  ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(preview->update(survey,options,SurveyExport::Format::GoFEM));
}

void ExportGOFEMDialog::accept() {
  updatePreview();
  if(ui->buttonBox->button(QDialogButtonBox::Ok)->isEnabled()) QDialog::accept();
}
