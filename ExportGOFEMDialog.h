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

#ifndef EXPORTGOFEMDIALOG_H
#define EXPORTGOFEMDIALOG_H

#include <QDialog>

#include "include/datum.h"

#include <set>

namespace Ui {
class ExportGOFEMDialog;
}

class ExportGOFEMDialog : public QDialog
{
  Q_OBJECT

public:
  explicit ExportGOFEMDialog(std::vector<double> periods, QWidget *parent = nullptr);
  ~ExportGOFEMDialog();

  std::vector<RealDataType> getSelectedDataTypes() const;
  std::vector<double> getSelectedPeriods() const;

private slots:
  void on_impGroupBox_toggled(bool on);
  void on_ptGroupBox_toggled(bool on);
  void on_tipperGroupBox_toggled(bool on);
  void on_rhoGroupBox_toggled(bool on);
  void on_phaseGroupBox_toggled(bool on);

  void on_periodListWidget_itemSelectionChanged();

private:
  Ui::ExportGOFEMDialog *ui;

  std::vector<double> m_periods;
};

#endif // EXPORTGOFEMDIALOG_H
