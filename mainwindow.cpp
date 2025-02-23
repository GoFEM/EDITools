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

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <fstream>
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

#include <QMessageBox>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>

#include "include/ApparentResistivityPlot.h"
#include "include/PhasePlot.h"
#include "include/PhaseTensorPlot.h"
#include "include/TipperPlot.h"
#include "ExportGOFEMDialog.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  ui->actionShow_error_bars->setChecked(true);

  QLabel* stationInfoLabel = new QLabel(ui->splitter_4);
  stationInfoLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
  stationInfoLabel->setMaximumHeight(20);
  stationInfoLabel->setAlignment(Qt::AlignHCenter);
  ui->splitter_4->insertWidget(0, stationInfoLabel);

  plot11 = new MyCustomPlot(ui->splitter_2);
  plot11->setObjectName(QStringLiteral("plot11"));
  ui->splitter_2->addWidget(plot11);
  plot12 = new MyCustomPlot(ui->splitter_2);
  plot12->setObjectName(QStringLiteral("plot12"));
  ui->splitter_2->addWidget(plot12);
  plot21 = new MyCustomPlot(ui->splitter);
  plot21->setObjectName(QStringLiteral("plot21"));
  ui->splitter->addWidget(plot21);
  plot22 = new MyCustomPlot(ui->splitter);
  plot22->setObjectName(QStringLiteral("plot22"));
  ui->splitter->addWidget(plot22);

  plotHandlers.push_back(std::shared_ptr<MTDataPlot>(new ApparentResistivityPlot(plot11)));
  plotHandlers.push_back(std::shared_ptr<MTDataPlot>(new PhasePlot(plot12)));
  plotHandlers.push_back(std::shared_ptr<MTDataPlot>(new TipperPlot(plot21)));
  plotHandlers.push_back(std::shared_ptr<MTDataPlot>(new PhaseTensorPlot(plot22)));

  plotHandlers[0]->set_associated_plot(*plotHandlers[1]);
  plotHandlers[1]->set_associated_plot(*plotHandlers[0]);

  mapHandler.reset(new MapPlot(ui->mapPlot, this));

  setupListContextMenu();

  lastDirectory = "/home/ag/USArray_edi";

  this->showMaximized();
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::on_actionLoad_EDI_triggered()
{
  if(mtSurvey == nullptr)
  {
    QString surveyName = QInputDialog::getText(this, tr("New survey"),
          tr("Survey name:"), QLineEdit::Normal, "Survey");

    mtSurvey.reset(new MTSurveyData(surveyName.toStdString()));

    this->setWindowTitle("MT Survey [" + surveyName + "]");
  }

  QStringList files = QFileDialog::getOpenFileNames(
                            this,
                            "Select one or more files to open",
                            lastDirectory,
                            "EDI files (*.edi)");

  if(files.size() > 0)
  {
    std::vector<std::string> file_list;
    for(auto &file: files)
      file_list.push_back(file.toStdString());

    try
    {
      auto duplicates = mtSurvey->load_from_edi(file_list);

      if(duplicates.size() != 0)
      {
        QString dupStr;
        for(auto &s: duplicates)
          dupStr += QString(s.c_str()) + "\n";

        QMessageBox msgBox;
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.setText(QString("Following stations had duplicates and were ignored: ") + dupStr);
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
      }
    }
    catch(std::exception &e)
    {
      QMessageBox msgBox;
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.setText(QString("Exception on file reading. The error message: ") + e.what());
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.exec();
    }
  }

  projectFile = "";

  createStationsList();
  updateMap();
}

void MainWindow::createStationsList()
{
  ui->stationList->blockSignals(true);
  ui->stationList->clear();
  ui->stationList->blockSignals(false);

  const auto names = mtSurvey->get_stations_names();
  for(unsigned i = 0; i < names.size(); ++i)
  {
    QListWidgetItem* item = new QListWidgetItem(QString(names[i].c_str()), ui->stationList);
    item->setCheckState(mtSurvey->is_active(names[i]) ? Qt::Checked : Qt::Unchecked);
  }

  ui->label->setText("Stations (" + QString("%1").arg(mtSurvey->n_stations()) + "):");
}

void MainWindow::updateMap()
{
  mapHandler->set_data(mtSurvey->get_stations_locations());
}

void MainWindow::updatePlots()
{
  QListWidgetItem* item = ui->stationList->currentItem();
  if(item == nullptr)
    return;

  MTStationData& data = mtSurvey->get_station_data(item->text().toStdString());
  for(auto &ph: plotHandlers)
    ph->set_observed_data(data);

  // Display calculated responses if any were loaded
  QListWidgetItem* resp_item = ui->responsesList->currentItem();
  if(resp_item != nullptr)
  {
    auto it = mtResponses.find(resp_item->toolTip().toStdString());
    if(it != mtResponses.end())
    {
      MTSurveyData &survey = it->second;
      if(survey.is_station_present(item->text().toStdString()))
      {
        MTStationData& resp_data = survey.get_station_data(item->text().toStdString());

        for(auto &ph: plotHandlers)
          ph->set_predicted_data(resp_data);

        ui->label->setText(QString("RMS = %1").arg(data.rms(resp_data)));
      }
    }
  }

  ui->statusBar->showMessage("# frequencies: " + QString("%1").arg(data.frequencies().size()) +
                             "\t\t\t# path: " + data.file_path().c_str());
}

void MainWindow::stationSelected(const QCPDataSelection &selection)
{
  if(!selection.isEmpty())
  {
    auto range = selection.dataRange();

    if(range.length() > 0)
    {
      double lon, lat;
      mapHandler->get_point_value(range.begin(), lon, lat);

      QString name(mtSurvey->closest_station_name(lat, lon).c_str());

      auto items = ui->stationList->findItems(name, Qt::MatchExactly);
      if(items.size() > 0)
        ui->stationList->setCurrentItem(items[0]);
    }
  }
}

void MainWindow::maskDataType(bool on)
{
  const std::string station_name = pointedItem->text().toStdString();
  QAction* action = qobject_cast<QAction*>(sender());
  const RealDataType type = static_cast<RealDataType>(action->data().toInt());

  mtSurvey->set_active_flag(station_name, type, on);

  updatePlots();
}

void MainWindow::deleteStation()
{
  ui->stationList->blockSignals(true);
  QListWidgetItem* item = ui->stationList->takeItem(ui->stationList->currentRow());
  ui->stationList->blockSignals(false);

  mtSurvey->remove_station(item->text().toStdString());

  ui->label->setText("Stations (" + QString("%1").arg(mtSurvey->n_stations()) + "):");

  updatePlots();
}

void MainWindow::renameStation()
{
  QListWidgetItem* item = ui->stationList->currentItem();

  bool ok;
  QString newName = QInputDialog::getText(this, tr("Name"),
                                          tr("New name:"), QLineEdit::Normal,
                                          item->text(), &ok);
  if (ok && !newName.isEmpty())
  {
    mtSurvey->rename_station(item->text().toStdString(), newName.toStdString());
    item->setText(newName);
  }
}

void MainWindow::on_actionSave_project_triggered()
{
  if(projectFile.length() < 1)
  {
    projectFile = QFileDialog::getSaveFileName(
          this,
          "Save project data",
          lastDirectory,
          "Project file (*.mtd)");

    if (!projectFile.endsWith(".mtd"))
      projectFile += ".mtd";
  }

  std::ofstream ofs(projectFile.toStdString());

  if(ofs.is_open())
  {
    boost::archive::binary_oarchive oa(ofs);
    oa << mtSurvey;
    oa << mtResponses;
  }

  ofs.close();
}

void MainWindow::on_actionLoad_project_triggered()
{
  projectFile = QFileDialog::getOpenFileName(
                            this,
                            "Open project file",
                            lastDirectory,
                            "Project file (*.mtd)");

  if(projectFile.length() == 0)
    return;

  std::ifstream ifs(projectFile.toStdString());

  if(ifs.is_open())
  {
    boost::archive::binary_iarchive ia(ifs);
    ia >> mtSurvey;
    ia >> mtResponses;
  }

  ifs.close();

  createStationsList();
  updateMap();
}

void MainWindow::on_stationList_currentRowChanged(int /*currentRow*/)
{
  updatePlots();

  std::vector<std::string> names;

  QListWidgetItem* item = ui->stationList->currentItem();
  names.push_back(item->text().toStdString());
}

void MainWindow::on_stationList_itemChanged(QListWidgetItem *item)
{
  bool checked = item->checkState() == Qt::Checked;
  mtSurvey->set_active_flag(item->text().toStdString(), checked);

  if(checked)
    item->setForeground(Qt::black);
  else
    item->setForeground(Qt::gray);

  auto selection = ui->stationList->selectedItems();
  for(auto sitem: selection)
  {
    sitem->setCheckState(item->checkState());
    mtSurvey->set_active_flag(sitem->text().toStdString(), checked);
  }

  updatePlots();
}

void MainWindow::on_actionDecimate_triggered()
{
  if(mtSurvey == nullptr)
    return;

  mtSurvey->decimate();
  updatePlots();
}

void MainWindow::on_stationList_itemSelectionChanged()
{
  std::vector<std::string> names;

  auto selection = ui->stationList->selectedItems();
  for(auto sitem: selection)
    names.push_back(sitem->text().toStdString());

  mapHandler->set_selected_points(mtSurvey->get_stations_locations(names));
}

void MainWindow::on_actionExport_in_GoFEM_triggered()
{
  if(mtSurvey == nullptr)
    return;

  QString exportFile = QFileDialog::getSaveFileName(
        this,
        "Export data",
        lastDirectory,
        "Project file (*.*)");

  if(exportFile.length() == 0)
    return;

  std::vector<double> periods = mtSurvey->get_unique_periods();

  ExportGOFEMDialog* dlg = new ExportGOFEMDialog(periods, this);

  if(dlg->exec())
  {
    auto selected_data_types = dlg->getSelectedDataTypes();
    auto selected_periods = dlg->getSelectedPeriods();
    mtSurvey->write_gofem(exportFile.toStdString(), selected_data_types, selected_periods);
  }
}

void MainWindow::on_actionSet_error_floor_triggered()
{
  bool ok;
  double error_floor = QInputDialog::getDouble(this, "Specify error floor",
                                               "Error floor (in %):",
                                               5, 0, 100, 2, &ok);
  if(ok && mtSurvey != nullptr)
    mtSurvey->set_error_floor(error_floor / 100.);

  updatePlots();
}

void MainWindow::on_actionShow_error_bars_toggled(bool on)
{
  for(auto &plot: plotHandlers)
    plot->set_error_bars_visible(on);

  updatePlots();
}

void MainWindow::on_stationList_customContextMenuRequested(const QPoint &pos)
{
  pointedItem = ui->stationList->itemAt(pos);
  if(pointedItem)
  {
    std::set<RealDataType> types = mtSurvey->get_active_types(pointedItem->text().toStdString());

    for(auto action: listContextActions)
    {
      RealDataType type = static_cast<RealDataType>(action->data().toInt());
      if(types.find(type) == types.end())
        action->setChecked(false);
      else
        action->setChecked(true);
    }

    listContextMenu->exec(ui->stationList->mapToGlobal(pos));
  }
}

void MainWindow::setupListContextMenu()
{
  listContextMenu = new QMenu;
  QMenu* mask_menu = new QMenu(tr("Mask"));
  listContextMenu->addMenu(mask_menu);
  listContextMenu->addAction(tr("Delete"), this, &MainWindow::deleteStation);
  listContextMenu->addAction(tr("Rename"), this, &MainWindow::renameStation);

  std::vector<RealDataType> types = {RealZxx, RealZxy, RealZyx, RealZyy,
                                     RealTzx, RealTzy,
                                     PTxx, PTxy, PTyx, PTyy};
  for(auto type: types)
  {
    QString type_name;
    if(Datum::get_parent_type(type) == InvalidComplexType)
      type_name = Datum::convert_type_to_string(type).c_str();
    else
      type_name = Datum::convert_type_to_string(Datum::get_parent_type(type)).c_str();

    QAction* action = mask_menu->addAction(type_name);
    action->setCheckable(true);
    action->setChecked(true);
    action->setData(int(type));
    connect(action, &QAction::triggered, this, &MainWindow::maskDataType);

    listContextActions.push_back(action);
  }
}

void MainWindow::on_actionActivate_masking_mode_toggled(bool on)
{
  for(auto &plot: plotHandlers)
    plot->set_masking_mode(on);
}

void MainWindow::on_actionSave_as_PDF_triggered()
{
  QListWidgetItem* item = ui->stationList->currentItem();

  if(!item)
    return;

  QString sname = item->text();

  QString exportFile = QFileDialog::getSaveFileName(
        this,
        "Save image",
        "./" + sname + ".png",
        "Image file (*.png)");

  if(exportFile.length() == 0)
    return;

  ui->splitter_4->grab().save(exportFile, "PNG", 100);
}

void MainWindow::on_actionLoad_GoFEM_responses_triggered()
{
  QString responseFile = QFileDialog::getOpenFileName(
                            this,
                            "Open GoFEM data file",
                            lastDirectory,
                            "Project file (*.*)");

  if(responseFile.length() == 0)
    return;

  QString fileName = QFileInfo(responseFile).fileName();

  try
  {
    MTSurveyData responses;
    responses.load_from_gofem(responseFile.toStdString());
    mtResponses.insert(std::make_pair(responseFile.toStdString(), responses));

    QListWidgetItem* item = new QListWidgetItem(fileName);
    item->setToolTip(responseFile);
    ui->responsesList->addItem(item);
  }
  catch(std::exception &e)
  {
    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Critical);
    msgBox.setText(QString("Exception on file reading. The error message: ") + e.what());
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
  }
}

void MainWindow::on_responsesList_currentRowChanged(int /*currentRow*/)
{
  updatePlots();
}

void MainWindow::on_actionAbout_EDI_Tools_triggered()
{
    QMessageBox::about(
        nullptr,
        "About",
        "<h3>Program Name: EDI Tools v1.0</h3>"
        "<p><b>Creator:</b> Alexander Grayver</p>"
        "<p><b>License:</b> GNU GPL 3.0</p>"
        );
}

