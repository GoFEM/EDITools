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
#include <vector>

#include "include/MTSurveyData.h"
#include "include/MapPlot.h"
#include "include/MTDataPlot.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = nullptr);
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
  void on_actionPhase_wrap_toggled(bool on);
  void on_actionShow_station_names_toggled(bool on);
  void on_actionTipper_arrows_toggled(bool on);
  void on_actionSave_as_PDF_triggered();
  void on_actionPlot_axis_ranges_triggered();

  void on_stationList_currentRowChanged(int currentRow);
  void on_stationList_itemChanged(QListWidgetItem *item);
  void on_stationList_itemSelectionChanged();
  void on_stationList_customContextMenuRequested(const QPoint &pos);

  void on_actionLoad_GoFEM_responses_triggered();

  void on_responsesList_currentRowChanged(int currentRow);

  void on_actionAbout_EDI_Tools_triggered();

  private:
  struct PlotAxisOptions
  {
    bool autoscale = true;
    double lower = 0;
    double upper = 0;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & autoscale;
      ar & lower;
      ar & upper;
    }
  };

  struct PlotOptions
  {
    bool phaseWrap = true;
    bool showStationNames = true;
    bool tipperArrows = false;
    std::vector<PlotAxisOptions> axes;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & phaseWrap;
      ar & showStationNames;
      ar & tipperArrows;
      ar & axes;
    }
  };

  struct PlotOptionsV2
  {
    bool phaseWrap = true;
    std::vector<PlotAxisOptions> axes;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & phaseWrap;
      ar & axes;
    }
  };

  struct PlotOptionsV3
  {
    bool phaseWrap = true;
    bool showStationNames = true;
    std::vector<PlotAxisOptions> axes;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & phaseWrap;
      ar & showStationNames;
      ar & axes;
    }
  };

  void setupListContextMenu();
  void rememberDirectory(const QString &path);
  PlotOptions currentPlotOptions() const;
  void applyPlotOptions(const PlotOptions &options);

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
