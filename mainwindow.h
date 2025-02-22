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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QListWidgetItem>
#include <QMenu>
#include <QAction>

#include <memory>

#include "MTSurveyData.h"
#include "MTDataPlot.h"
#include "MapPlot.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

private slots:
  void createStationsList();
  void updateMap();
  void updatePlots();
  void stationSelected(const QCPDataSelection &selection);

  void maskDataType(bool on);
  void deleteStation();
  void renameStation();

  void on_actionLoad_EDI_triggered();
  void on_actionSave_project_triggered();
  void on_actionLoad_project_triggered();
  void on_actionDecimate_triggered();
  void on_actionExport_in_GoFEM_triggered();
  void on_actionSet_error_floor_triggered();
  void on_actionShow_error_bars_toggled(bool on);

  void on_actionActivate_masking_mode_toggled(bool on);
  void on_actionSave_as_PDF_triggered();

  void on_stationList_currentRowChanged(int currentRow);
  void on_stationList_itemChanged(QListWidgetItem *item);
  void on_stationList_itemSelectionChanged();
  void on_stationList_customContextMenuRequested(const QPoint &pos);

  void on_actionLoad_GoFEM_responses_triggered();

  void on_responsesList_currentRowChanged(int currentRow);

  void on_actionAbout_EDI_Tools_triggered();

  private:
  void setupListContextMenu();

private:
  Ui::MainWindow *ui;

  std::shared_ptr<MTSurveyData> mtSurvey;
  std::map<std::string, MTSurveyData> mtResponses;
  std::vector<std::shared_ptr<MTDataPlot>> plotHandlers;
  std::shared_ptr<MapPlot> mapHandler;

  QMenu* listContextMenu;
  QVector<QAction*> listContextActions;

  QListWidgetItem* pointedItem;

  QCustomPlot *plot11;
  QCustomPlot *plot12;
  QCustomPlot *plot21;
  QCustomPlot *plot22;

  QString projectFile, lastDirectory;

  QLabel *stationInfoLabel;
};

#endif // MAINWINDOW_H
