#include "include/MapBackground.h"
#include "TopographyDialog.h"
#include "include/FileLabels.h"
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
#include "ImpedanceAnalysisWindows.h"
#include "include/RhoPlusOverlay.h"
#include "ui_mainwindow.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>

#include <QMessageBox>
#include <QSettings>
#include <QFileInfo>
#include <QDir>
#include <QCheckBox>
#include <QDialog>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QSignalBlocker>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>

#include "include/ApparentResistivityPlot.h"
#include "include/PhasePlot.h"
#include "include/PhaseTensorPlot.h"
#include "include/TipperPlot.h"
#include "ExportGOFEMDialog.h"
#include "ExportNativeMTDialog.h"
#include "ExportPreview.h"
#include "StationComparisonWindow.h"
#include "PolarityWindow.h"
#include "GroupMaskDialog.h"
#include "AutoMaskWorkflow.h"
#include "SurveyCoordinatesDialog.h"
#include "include/SurveyHistory.h"
#include "include/CommandRibbon.h"
#include "include/TFViewControls.h"
#include "include/StationListDelegate.h"
#include "include/CommandIcons.h"
#include <QUndoView>
#include <QVBoxLayout>
#include "FitStatisticsWindow.h"
#include "PeriodMapWindow.h"
#include "ProfileWindow.h"
#include "PeriodLayoutWindow.h"
#include "EDIImportDialog.h"
#include "include/SurveyReport.h"

namespace
{
// Project files created before plot options were added start directly with the
// serialized survey object. New files start with this marker and a format
// version so the loader can distinguish both layouts and keep old .mtd files
// readable.
const std::uint32_t projectFileMagic = 0x45444954; // "EDIT"
const std::uint32_t projectFileVersion = 9;
}

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  ui->stationList->setItemDelegate(new StationListDelegate(ui->stationList));
  ui->stationList->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  history.reset(new SurveyHistory([this] { restoreSurveyView(); }));
  auto *edit = new QMenu(tr("Edit"), this);
  ui->menuBar->insertMenu(ui->menuData->menuAction(), edit);
  auto *undo = history->stack().createUndoAction(this, tr("Undo"));
  undo->setObjectName("actionUndo"); undo->setShortcut(QKeySequence::Undo); edit->addAction(undo);
  auto *redo = history->stack().createRedoAction(this, tr("Redo"));
  redo->setObjectName("actionRedo"); redo->setShortcuts(QKeySequence::Redo); edit->addAction(redo);
  auto *historyDialog = new QDialog(this, Qt::Window);
  historyDialog->setObjectName("surveyHistory"); historyDialog->setWindowTitle(tr("Survey history"));
  auto *historyLayout = new QVBoxLayout(historyDialog);
  auto *view = new QUndoView(&history->stack(), historyDialog);
  view->setEmptyLabel(tr("Start of retained history")); historyLayout->addWidget(view);
  view->setToolTip(tr("Select a step to undo or redo edits. The latest 20 edits are retained for this session. Loading a project starts a new history. Resampling opens a separate survey with its own history."));
  auto *original = new QPushButton(tr("Open original survey"), historyDialog); historyLayout->addWidget(original);
  connect(original, &QPushButton::clicked, this, [this] { if(auto data=history->original()) openSurveyCopy(std::move(data)); });
  historyDialog->resize(420, 360);
  auto *historyAction = edit->addAction(tr("History…")); historyAction->setObjectName("actionSurvey_history");
  connect(historyAction, &QAction::triggered, historyDialog, [historyDialog] { historyDialog->show(); historyDialog->raise(); });
  ui->responsesList->setTextElideMode(Qt::ElideLeft);
  ui->responsesList->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  auto *responsesPanel = ui->responsesList->parentWidget();
  ui->splitter_3->setCollapsible(ui->splitter_3->indexOf(responsesPanel), false);
  const auto setResponsesExpanded = [this, responsesPanel](bool expanded) {
    ui->responsesToggle->setArrowType(expanded ? Qt::DownArrow : Qt::RightArrow);
    ui->responsesList->setVisible(expanded);
    responsesPanel->setMaximumHeight(expanded ? QWIDGETSIZE_MAX : responsesPanel->sizeHint().height());
    if(expanded) {
      auto sizes = ui->splitter_3->sizes();
      const int index = ui->splitter_3->indexOf(responsesPanel);
      sizes[index] = std::max(sizes[index], 150);
      ui->splitter_3->setSizes(sizes);
    }
  };
  connect(ui->responsesToggle, &QToolButton::toggled, this, setResponsesExpanded);
  setResponsesExpanded(false);
  ui->actionShow_error_bars->setChecked(true);
  ui->actionPhase_wrap->setChecked(true);

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
  rhoPlusAction = ui->menuPlot->addAction(tr("Show rho+ guides"));
  rhoPlusAction->setObjectName("actionShow_rho_plus"); rhoPlusAction->setCheckable(true);
  const auto rhoPlusHelp = tr("Thin XY/YX rho+ fits in the component colors: dashed for stable fits, dotted for provisional fits. Fits use enabled, unmasked observations and refresh after masking. Component checkboxes also hide the guides. These are 1-D consistency guides, not automatic validity judgments. Change grid density in the Rho+ analysis window.");
  rhoPlusAction->setToolTip(rhoPlusHelp); ui->menuPlot->setToolTipsVisible(true);
  rhoPlusOverlay.reset(new RhoPlusOverlay(plot11,plot12,[this,rhoPlusHelp](const QString &status) {
    rhoPlusAction->setToolTip(rhoPlusHelp + (status.isEmpty()?QString():"\n\n"+status));
  }));
  connect(rhoPlusAction,&QAction::toggled,this,[this]{updatePlots(false);});
  for(const auto &plot: plotHandlers) {
    connect(plot.get(), &MTDataPlot::observationsChanged, this, [this] { history->record(tr("Edit masks")); updatePlots(false); refreshAnalysisWindows(); });
    connect(plot.get(), &MTDataPlot::componentVisibilityChanged, this, [this] { updatePlots(); });
  }

  mapHandler.reset(new MapPlot(ui->mapPlot, this));
  mapBackground = new MapBackground(ui->mapPlot);
  auto *topography = ui->menuCoordinates->addAction(tr("Topography…"));
  topography->setObjectName("actionTopography");
  connect(topography, &QAction::triggered, this, [this] {
    if(!mtSurvey) { QMessageBox::information(this,tr("Topography"),tr("Load a survey first.")); return; }
    try {
      TopographyDialog dialog(this,*mtSurvey,[this] {
        mapBackground->setTopography(mtSurvey->topography());ui->mapPlot->replot();refreshAnalysisWindows();
      });
      dialog.exec();
      history->record(tr("Topography / elevations"));
    } catch(const std::exception &e) { QMessageBox::warning(this,tr("Topography"),QString::fromUtf8(e.what())); }
  });
  auto *mapControls = new QWidget(ui->splitter_3);
  auto *mapControlsLayout = new QVBoxLayout(mapControls);
  mapControlsLayout->setContentsMargins(4, 0, 4, 0);
  auto *mapActions = new QHBoxLayout;
  mapStationNames = new QCheckBox(tr("Station names"), mapControls);
  mapStationNames->setObjectName("stationMapNames");
  mapActions->addWidget(mapStationNames);
  auto *layersButton = MapBackground::button(mapControls, mapBackground->layers(), [this](int layers) { mapBackground->setLayers(layers); ui->mapPlot->replot(); });
  mapActions->addWidget(layersButton);
  auto *layersAction = ui->menuCoordinates->addAction(tr("Map layers…"));
  layersAction->setObjectName("actionMap_layers");
  layersAction->setMenu(layersButton->menu());
  layersAction->setToolTip(layersButton->toolTip());
  auto *terrainButton = new QPushButton(tr("Topography…"),mapControls);
  terrainButton->setObjectName("openTopography");mapActions->addWidget(terrainButton);
  connect(terrainButton,&QPushButton::clicked,topography,&QAction::trigger);
  connect(topography, &QAction::changed, terrainButton, [terrainButton, topography] { terrainButton->setEnabled(topography->isEnabled()); });
  mapActions->addStretch();
  auto *mapsButton = new QPushButton(tr("Maps…"), mapControls);
  mapsButton->setObjectName("openPeriodMaps");
  mapsButton->setToolTip(tr("Induction vectors, phase tensors, apparent resistivity and phase at a selected period"));
  mapActions->addWidget(mapsButton);
  mapControlsLayout->addLayout(mapActions);
  connect(mapStationNames, &QCheckBox::toggled, ui->actionShow_station_names, &QAction::setChecked);
  connect(mapsButton, &QPushButton::clicked, ui->actionPeriod_maps, &QAction::trigger);
  connect(ui->actionPeriod_maps, &QAction::changed, mapsButton, [this, mapsButton] { mapsButton->setEnabled(ui->actionPeriod_maps->isEnabled()); });
  auto *profileAction = ui->menuCoordinates->addAction(tr("2D profile analysis…"));
  profileAction->setObjectName("actionProfile_analysis");
  profileAction->setToolTip(tr("Draw a profile corridor and view resistivity, phase, tensor sections and strike roses."));
  auto *profileButton = new QPushButton(tr("2D profile…"), mapControls);
  profileButton->setObjectName("openProfile");
  profileButton->setToolTip(profileAction->toolTip());
  mapActions->addWidget(profileButton);
  connect(profileButton, &QPushButton::clicked, profileAction, &QAction::trigger);
  connect(profileAction, &QAction::changed, profileButton, [profileButton, profileAction] { profileButton->setEnabled(profileAction->isEnabled()); });
  connect(profileAction, &QAction::triggered, this, [this] {
    const bool first = !profile;
    if(first) profile = new ProfileWindow(this, [this](std::shared_ptr<MTSurveyData> data) { openSurveyCopy(std::move(data)); });
    profile->setData(mtSurvey, mtResponses);
    if(first) profile->setMapLayers(mapBackground->layers());
    profile->show(); profile->raise(); profile->activateWindow();
  });
  auto *periodLayoutAction = ui->menuData->addAction(tr("Period layout / Resampling…"));
  periodLayoutAction->setObjectName("actionPeriod_layout");
  connect(periodLayoutAction, &QAction::triggered, this, &MainWindow::openPeriodLayout);
  auto *compare = ui->menuData->addAction(tr("Compare neighbouring stations…"));
  compare->setObjectName("actionCompare_stations");
  connect(compare, &QAction::triggered, this, [this] {
    if(!stationComparison) stationComparison = new StationComparisonWindow(this);
    QStringList selected;
    if(auto *item=ui->stationList->currentItem()) selected << item->text();
    for(auto *item:ui->stationList->selectedItems()) if(!selected.contains(item->text()))selected << item->text();
    stationComparison->setData(mtSurvey,selected);
    stationComparison->show();stationComparison->raise();stationComparison->activateWindow();
  });
  auto *validityAction = ui->menuData->addAction(tr("Impedance validity / Rho+…"));
  auto *polarityAction = ui->menuData->addAction(tr("Polarity review…"));
  polarityAction->setObjectName("actionPolarity_review");
  polarityAction->setToolTip(tr("Review likely polarity reversals against neighbouring stations, preview channel corrections and apply with Undo support."));
  connect(polarityAction, &QAction::triggered, this, [this] {
    if(!polarity) polarity = new PolarityWindow(this, *history, [this] { updatePlots(false); refreshAnalysisWindows(); });
    polarity->setData(mtSurvey, ui->stationList->currentItem() ? ui->stationList->currentItem()->text() : QString());
    polarity->show(); polarity->raise(); polarity->activateWindow();
  });
  validityAction->setObjectName("actionImpedance_validity");
  connect(validityAction, &QAction::triggered, this, [this] {
    if(!impedanceValidity) impedanceValidity = new ImpedanceValidityWindow(this, [this](const ImpedanceAnalysis::RhoPlusOptions &options) {
      rhoPlusOptions=options; updatePlots(false);
    }, [this] { updatePlots(false); refreshAnalysisWindows(); });
    impedanceValidity->setLinkTensorMasks(linkTensorMasksAction->isChecked());
    impedanceValidity->setFitOptions(rhoPlusOptions);
    impedanceValidity->setHistory(history.get());
    impedanceValidity->setData(mtSurvey, ui->stationList->currentItem() ? ui->stationList->currentItem()->text() : QString());
    impedanceValidity->show(); impedanceValidity->raise(); impedanceValidity->activateWindow();
  });
  auto *averagesAction = ui->menuData->addAction(tr("Survey impedance averages…"));
  averagesAction->setObjectName("actionSurvey_averages");
  connect(averagesAction, &QAction::triggered, this, [this] {
    if(!surveyAverage) surveyAverage = new SurveyAverageWindow(this);
    surveyAverage->setData(mtSurvey);
    surveyAverage->show(); surveyAverage->raise(); surveyAverage->activateWindow();
  });
  linkTensorMasksAction = ui->menuData->addAction(tr("Link Z / PT masks"));
  linkTensorMasksAction->setObjectName("actionLink_tensor_masks");
  linkTensorMasksAction->setCheckable(true); linkTensorMasksAction->setChecked(true);
  linkTensorMasksAction->setToolTip(tr("Link mask and unmask edits at the same period: a selected Z component affects all PT entries; a selected PT entry affects all Z components. Phi = Re(Z)^-1 Im(Z), so there is no one-to-one mapping in general 3D. Independently, Zxy edits also affect Zxx, and Zyx edits also affect Zyy. Tippers remain independent. Toggling this option does not alter existing masks."));
  connect(linkTensorMasksAction, &QAction::toggled, this, [this](bool on) {
    for(const auto &plot: plotHandlers) plot->set_link_tensor_masks(on);
    if(periodMaps) periodMaps->setLinkTensorMasks(on);
    if(impedanceValidity) impedanceValidity->setLinkTensorMasks(on);
  });
  auto *groupMasks = edit->addAction(tr("Group masking…"));
  groupMasks->setObjectName("actionGroup_masking");
  groupMasks->setToolTip(tr("Mask or unmask chosen periods and components across the selected stations, with one Undo step."));
  connect(groupMasks, &QAction::triggered, this, [this] {
    if(!mtSurvey) return;
    std::set<std::string> names;
    for(const auto *item : ui->stationList->selectedItems()) names.insert(item->text().toStdString());
    if(names.empty() || !editGroupMasks(this, *mtSurvey, names, linkTensorMasksAction->isChecked())) return;
    history->record(tr("Group masking"));
    updatePlots(false); refreshAnalysisWindows();
  });
  autoMaskWorkflow = new AutoMaskWorkflow(*this, *history, plotHandlers, linkTensorMasksAction);
  connect(autoMaskWorkflow, &AutoMaskWorkflow::masksApplied, this, [this] { updatePlots(false); refreshAnalysisWindows(); });
  connect(autoMaskWorkflow, &AutoMaskWorkflow::navigateRequested, this, [this](int direction) {
    const int row = ui->stationList->currentRow() + direction;
    if(row >= 0 && row < ui->stationList->count()) ui->stationList->setCurrentRow(row, QItemSelectionModel::NoUpdate);
  });
  auto *autoMenu = new QMenu(tr("Auto mask"), this); ui->menuBar->insertMenu(ui->menuData->menuAction(), autoMenu);
  for(const char *name : {"actionAuto_mask_preview", "actionAuto_mask_settings", "actionAuto_masking", "actionAuto_mask_previous", "actionAuto_mask_next",
       "actionAuto_mask_apply", "actionAuto_mask_apply_next", "actionAuto_mask_selected", "actionAuto_mask_check0", "actionAuto_mask_check1",
       "actionAuto_mask_check2", "actionAuto_mask_check3", "actionAuto_mask_XY", "actionAuto_mask_YX", "actionAuto_mask_Tzx", "actionAuto_mask_Tzy", "actionAuto_mask_pair"})
    autoMenu->addAction(findChild<QAction *>(name));
  auto *reportAction = ui->menuFile->addAction(tr("Export survey report…"));
  reportAction->setObjectName("actionExport_survey_report");
  connect(reportAction, &QAction::triggered, this, &MainWindow::exportSurveyReport);
  mapCoordinateSwitch = new QComboBox(mapControls);
  mapCoordinateSwitch->setObjectName("mapCoordinateSwitch");
  mapCoordinateSwitch->addItems({tr("Lat/Lon"), tr("UTM")});
  mapCoordinateSwitch->setToolTip(tr("Map coordinates only; export coordinates and the UTM origin are preserved."));
  auto *mapSwitchLayout = new QHBoxLayout;
  mapSwitchLayout->addWidget(new QLabel(tr("Map coordinates:"), mapControls));
  mapSwitchLayout->addWidget(mapCoordinateSwitch);
  mapCenterUtm = new QCheckBox(tr("Center UTM"), mapControls);
  mapCenterUtm->setObjectName("mapCenterUtm");
  mapCenterUtm->setEnabled(false);
  mapCenterUtm->setToolTip(tr("Subtract the midpoint of the projected survey bounds. The origin is saved with the survey and shared by maps, reports and exports. The chosen UTM zone is preserved; centering does not change projection distortion."));
  mapSwitchLayout->addWidget(mapCenterUtm);
  mapControlsLayout->addLayout(mapSwitchLayout);
  coordinateInfoLabel = new QLabel(tr("Latitude / longitude (degrees)"), mapControls);
  coordinateInfoLabel->setObjectName("mapCoordinateInfo");
  coordinateInfoLabel->setWordWrap(true);
  coordinateInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
  coordinateInfoLabel->setMaximumHeight(85);
  mapControlsLayout->addWidget(coordinateInfoLabel);
  mapControls->setMaximumHeight(155);
  ui->splitter_3->addWidget(mapControls);
  ui->splitter_3->setCollapsible(ui->splitter_3->indexOf(mapControls), false);
  connect(mapCoordinateSwitch, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] {
    try {
      updateMap();
    } catch(const std::exception &e) {
      const QSignalBlocker blocker(mapCoordinateSwitch);
      mapCoordinateSwitch->setCurrentIndex(mapCoordinates.utm ? 1 : 0);
      QMessageBox::warning(this, tr("Map coordinates"), QString::fromUtf8(e.what()));
    }
  });
  connect(mapCenterUtm, &QCheckBox::toggled, this, [this](bool centered) {
    if(!mtSurvey) return;
    try {
      const auto geographic = mtSurvey->geographic_locations();
      const auto settings = mtSurvey->coordinates().utm ? mtSurvey->coordinates() :
        SurveyCoordinates::suggested(geographic);
      mtSurvey->set_coordinates(SurveyCoordinates::calculate(geographic, settings.zone,
                                                            settings.north, centered));
      history->record(tr("Center UTM coordinates"));
      updateMap();
      refreshAnalysisWindows();
    } catch(const std::exception &e) {
      const QSignalBlocker blocker(mapCenterUtm);
      mapCenterUtm->setChecked(mtSurvey->coordinates().centered);
      QMessageBox::warning(this, tr("Map coordinates"), QString::fromUtf8(e.what()));
    }
  });
  ui->actionShow_station_names->setChecked(false);

  setupListContextMenu();

  lastDirectory = QSettings().value("lastDirectory", QDir::homePath()).toString();

  setupTransferFunctionView();
  for(auto *plot: {plot11, plot12, plot21, plot22})
    connect(plot->xAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, [this](const QCPRange &range) {
      if(syncingPeriods) return;
      syncingPeriods = true;
      for(auto *other: {plot11, plot12, plot21, plot22}) {
        other->xAxis->setRange(range); other->replot(QCustomPlot::rpQueuedReplot);
      }
      syncingPeriods = false;
      if(viewControls) viewControls->refreshRange();
    });
  setupRibbon();
  for(auto *plot: {plot11, plot12}) for(auto *rect: plot->axisRects())
    connect(rect->axis(QCPAxis::atLeft), qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this,
            [this] { viewControls->refreshRange(); });
  this->showMaximized();
}

void MainWindow::setupTransferFunctionView()
{
  tfView = new QWidget;
  tfView->setObjectName("transferFunctionView");
  tfView->setSizePolicy(ui->splitter_4->sizePolicy());
  ui->splitter_5->replaceWidget(ui->splitter_5->indexOf(ui->splitter_4), tfView);
  auto *layout = new QVBoxLayout(tfView);
  layout->setContentsMargins(0, 0, 0, 0); layout->setSpacing(0);
  stationHeading = new QLabel(tr("No station selected"), tfView);
  stationHeading->setObjectName("tfStationName");
  stationHeading->setTextFormat(Qt::PlainText);
  stationHeading->setAlignment(Qt::AlignCenter);
  stationHeading->setMargin(5);
  stationHeading->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Fixed);
  auto font = stationHeading->font(); font.setBold(true);
  if(font.pointSizeF() > 0) font.setPointSizeF(font.pointSizeF() + 2);
  stationHeading->setFont(font);
  layout->addWidget(stationHeading); layout->addWidget(ui->splitter_4, 1);
  ui->splitter_2->setOrientation(Qt::Vertical);
  ui->splitter->setOrientation(Qt::Vertical);
  for(int i = 0; i < 2; ++i) plotHandlers[i]->set_split_components(true);
  for(int column = 0; column < 2; ++column) {
    auto *margins = new QCPMarginGroup(plot11);
    for(auto *plot: {plot11, plot12})
      plot->axisRect(column)->setMarginGroup(QCP::msLeft | QCP::msRight, margins);
  }
  static_cast<PhaseTensorPlot *>(plotHandlers[3].get())->set_ellipse_mode(true);
  static_cast<TipperPlot *>(plotHandlers[2].get())->set_arrow_mode(true);
  ui->splitter_4->setSizes({650, 350});
  ui->splitter_2->setSizes({350, 300});
  ui->splitter->setSizes({170, 180});
}

void MainWindow::setupRibbon()
{
  const auto command = [this](const char *name, const QString &label, const char *glyph = nullptr, bool needsSurvey = false) {
    auto *action = findChild<QAction *>(QLatin1String(name));
    Q_ASSERT(action);
    action->setIconText(label);
    if(glyph) action->setIcon(commandIcon(glyph));
    if(needsSurvey) surveyActions.push_back(action);
    return action;
  };
  auto *open = command("actionLoad_project", tr("Open"), "folder-open");
  auto *save = command("actionSave_project", tr("Save"), "save", true);
  auto *undo = command("actionUndo", tr("Undo"), "arrow-undo");
  auto *redo = command("actionRedo", tr("Redo"), "arrow-redo");
  viewControls = new TFViewControls(this, [this](int i) { return viewRange(i); },
    [this](int i, TFViewRange range) { setViewRange(i, range); }, [this] { applyViewStyle(); });
  auto *fit = ui->menuPlot->addAction(tr("Fit all TF axes"));
  fit->setObjectName("actionFit_TF_axes"); fit->setIcon(commandIcon("arrow-expand"));
  connect(fit, &QAction::triggered, this, [this] {
    for(int i = 0; i < 5; ++i) { auto range = viewRange(i); range.automatic = true; setViewRange(i, range); }
  });
  auto *reset = ui->menuPlot->addAction(tr("Reset TF appearance"));
  reset->setObjectName("actionReset_TF_appearance"); reset->setIcon(commandIcon("arrow-undo"));
  connect(reset, &QAction::triggered, this, [this] { viewControls->setStyle({}); applyViewStyle(); });
  ui->actionPlot_axis_ranges->setText(tr("Configure TF view"));
  auto viewGroups = viewControls->groups();
  viewGroups.push_back({tr("Display"), {command("actionShow_error_bars", tr("Error bars"), "ruler"),
    command("actionPhase_wrap", tr("Wrap phase"), "arrow-repeat-all"), fit, reset}});
  const std::vector<CommandRibbon::Page> pages{
    {tr("&Survey"), {
      {tr("Project"), {open, save}},
      {tr("Import"), {command("actionLoad_EDI", tr("EDI files"), "document-arrow-down"), command("actionLoad_responses", tr("Responses"), "data-line")}},
      {tr("Export data"), {command("actionExport_native_MT", tr("Native MT"), "document-arrow-up", true), command("actionExport_in_GoFEM", tr("GoFEM"), "document-arrow-up", true)}}}},
    {tr("&Edit"), {
      {tr("History"), {undo, redo, command("actionSurvey_history", tr("History"), "history")}},
      {tr("Masks"), {command("actionActivate_masking_mode", tr("Masking mode"), "edit"), command("actionLink_tensor_masks", tr("Link Z / PT"), "link"), command("actionGroup_masking", tr("Group masking"), "select-all-on", true)}},
      {tr("Periods"), {command("actionPeriod_layout", tr("Layout / resample"), "table", true), command("actionDecimate", tr("Decimate"), "filter", true)}},
      {tr("Errors"), {command("actionSet_error_floor", tr("Error floors"), "options", true)}}}},
    {tr("Auto &mask"), {
      {tr("Preview"), {command("actionAuto_mask_preview", tr("Preview"), "eye"), command("actionAuto_mask_settings", tr("Settings"), "settings"), command("actionAuto_masking", tr("Detailed review"), "table")}},
      {tr("Station"), {command("actionAuto_mask_previous", tr("Previous"), "chevron-left"), command("actionAuto_mask_next", tr("Next"), "chevron-right"), command("actionAuto_mask_apply", tr("Apply"), "checkmark"), command("actionAuto_mask_apply_next", tr("Apply + next"), "checkmark-circle")}},
      {tr("Batch"), {command("actionAuto_mask_selected", tr("Apply to selected"), "select-all-on"), command("actionAuto_mask_pair", tr("Pair diagonals"), "link"), command("actionLink_tensor_masks", tr("Link Z / PT"), "link")}},
      {tr("Checks"), {command("actionAuto_mask_check0", tr("Rho+")), command("actionAuto_mask_check1", tr("Period outliers")), command("actionAuto_mask_check2", tr("Neighbours")), command("actionAuto_mask_check3", tr("Uncertainty"))}},
      {tr("Components"), {command("actionAuto_mask_XY", tr("XY")), command("actionAuto_mask_YX", tr("YX")), command("actionAuto_mask_Tzx", tr("Tzx")), command("actionAuto_mask_Tzy", tr("Tzy"))}}}},
    {tr("&Analysis"), {
      {tr("Validity"), {command("actionImpedance_validity", tr("Rho+ analysis"), "data-line", true), command("actionShow_rho_plus", tr("Rho+ guides"), "data-line"), command("actionPolarity_review", tr("Polarity review"), "arrow-swap", true)}},
      {tr("Compare"), {command("actionCompare_stations", tr("Neighbours"), "people-community", true), command("actionSurvey_averages", tr("Survey averages"), "math-symbols", true)}},
      {tr("Data fit"), {command("actionFit_statistics", tr("Statistics / RMS"), "data-histogram", true)}}}},
    {tr("&Maps"), {
      {tr("Coordinates"), {command("actionConvert_UTM", tr("UTM / origin"), "my-location", true), command("actionGeographic_coordinates", tr("Geographic"), "globe", true)}},
      {tr("Survey maps"), {command("actionPeriod_maps", tr("PT / induction"), "map", true), command("actionProfile_analysis", tr("2D profile"), "map", true), command("actionShow_station_names", tr("Station names"), "text-field")}},
      {tr("Background"), {command("actionMap_layers", tr("Map layers"), "layer"), command("actionTopography", tr("Topography"), "mountain-location-top", true)}}}},
    {tr("&Report"), {
      {tr("Export figures"), {command("actionExport_survey_report", tr("Survey PDF"), "document-pdf", true), command("actionSave_as_PDF", tr("Plot image"), "image")}},
      {tr("Plot display"), {command("actionPlot_axis_ranges", tr("Configure TF view"), "options")}}}},
    {tr("&View"), viewGroups}
  };
  ui->actionAbout_EDI_Tools->setIcon(commandIcon("info"));
  ribbon = new CommandRibbon(*this, {open, save, undo, redo}, pages, ui->actionAbout_EDI_Tools);
  updateActionAvailability();
}

void MainWindow::updateActionAvailability()
{
  const bool hasData = mtSurvey && mtSurvey->n_stations();
  for(auto *action : surveyActions) action->setEnabled(hasData);
  if(auto *action = findChild<QAction *>("actionGroup_masking"))
    action->setEnabled(hasData && !ui->stationList->selectedItems().isEmpty());
  ui->actionSave_as_PDF->setEnabled(hasData && ui->stationList->currentItem());
  if(autoMaskWorkflow) {
    std::set<std::string> selected;
    for(const auto *item : ui->stationList->selectedItems()) selected.insert(item->text().toStdString());
    autoMaskWorkflow->setSelection(selected);
  }
}

void MainWindow::exportSurveyReport()
{
  if(!mtSurvey) { QMessageBox::information(this, tr("Survey report"), tr("Load a survey first.")); return; }
  SurveyReport::Options options;
  options.errorBars = ui->actionShow_error_bars->isChecked();
  const auto plots = currentPlotOptions();
  options.phaseWrap = plots.phaseWrap;
  options.tipperArrows = plots.tipperArrows;
  options.mapLayers = mapBackground->layers();
  if(periodMaps) { options.mapSettings = periodMaps->displaySettings(); options.useMapSettings = true; }
  options.components = plots.componentVisible;
  for(unsigned i = 0; i < options.axes.size() && i < plots.axes.size(); ++i)
    options.axes[i] = {plots.axes[i].autoscale, plots.axes[i].lower, plots.axes[i].upper};
  const auto *selected = ui->responsesList->currentItem();
  SurveyReport::show_dialog(this, *mtSurvey, mtResponses, selected ? selected->toolTip() : QString(), options, lastDirectory);
  QSettings().setValue("lastDirectory", lastDirectory);
}

void MainWindow::rememberDirectory(const QString &path)
{
  if(path.isEmpty())
    return;

  const QFileInfo info(path);
  const QString dir = info.isDir() ? info.absoluteFilePath() : info.absolutePath();
  if(dir.isEmpty())
    return;

  lastDirectory = dir;
  QSettings().setValue("lastDirectory", lastDirectory);
}

MainWindow::PlotOptions MainWindow::currentPlotOptions() const
{
  PlotOptions options;
  options.style = viewControls->style();
  options.period = periodRange;
  options.errorBars = ui->actionShow_error_bars->isChecked();
  for(unsigned i = 0; i < 2; ++i) {
    const auto range = plotHandlers[i]->fixed_y_axis_range(true);
    options.diagonals[i] = {plotHandlers[i]->y_axis_autoscale(true), range.lower, range.upper};
  }
  options.rhoPlus = rhoPlusAction->isChecked();
  options.rhoPlusGridDensity = rhoPlusOptions.intervalsPerDecade;
  options.linkTensorMasks = linkTensorMasksAction->isChecked();
  for(unsigned i = 0; i < plotHandlers.size(); ++i)
    options.componentVisible[i] = plotHandlers[i]->component_visibility();
  options.axes.reserve(plotHandlers.size());

  for(const auto &plot: plotHandlers)
  {
    const QCPRange range = plot->fixed_y_axis_range();
    PlotAxisOptions axisOptions;
    axisOptions.autoscale = plot->y_axis_autoscale();
    axisOptions.lower = range.lower;
    axisOptions.upper = range.upper;
    options.axes.push_back(axisOptions);
  }

  if(plotHandlers.size() > 1)
  {
    if(auto *phasePlot = dynamic_cast<PhasePlot*>(plotHandlers[1].get()))
      options.phaseWrap = phasePlot->phase_wrap();
  }
  options.showStationNames = ui->actionShow_station_names->isChecked();
  if(plotHandlers.size() > 2)
  {
    if(auto *tipperPlot = dynamic_cast<TipperPlot*>(plotHandlers[2].get()))
      options.tipperArrows = tipperPlot->arrow_mode();
  }

  return options;
}

void MainWindow::applyPlotOptions(const PlotOptions &options)
{
  viewControls->setStyle(options.style);
  applyViewStyle();
  periodRange = options.period;
  if(!std::isfinite(periodRange.lower) || !std::isfinite(periodRange.upper) || periodRange.lower <= 0. ||
     !QCPRange::validRange(periodRange.lower, periodRange.upper) || periodRange.lower >= periodRange.upper)
    periodRange = {};
  for(unsigned i = 0; i < 2; ++i) {
    plotHandlers[i]->set_y_axis_range(options.diagonals[i].lower, options.diagonals[i].upper, true);
    plotHandlers[i]->set_y_axis_autoscale(options.diagonals[i].automatic, true);
  }
  const QSignalBlocker errorsBlock(ui->actionShow_error_bars);
  ui->actionShow_error_bars->setChecked(options.errorBars);
  for(auto &plot: plotHandlers) plot->set_error_bars_visible(options.errorBars);
  const QSignalBlocker rhoPlusBlock(rhoPlusAction);
  rhoPlusAction->setChecked(options.rhoPlus);
  rhoPlusOptions.intervalsPerDecade = options.rhoPlusGridDensity >= 4 && options.rhoPlusGridDensity <= 48 ? options.rhoPlusGridDensity : 12;
  if(impedanceValidity) impedanceValidity->setFitOptions(rhoPlusOptions);
  linkTensorMasksAction->setChecked(options.linkTensorMasks);
  const QSignalBlocker phaseBlock(ui->actionPhase_wrap);
  const QSignalBlocker namesBlock(ui->actionShow_station_names);
  const size_t nAxes = std::min(options.axes.size(), plotHandlers.size());
  for(size_t i = 0; i < nAxes; ++i)
  {
    plotHandlers[i]->set_y_axis_range(options.axes[i].lower, options.axes[i].upper);
    plotHandlers[i]->set_y_axis_autoscale(options.axes[i].autoscale);
  }

  ui->actionPhase_wrap->setChecked(options.phaseWrap);
  ui->actionShow_station_names->setChecked(options.showStationNames);
  {
    const QSignalBlocker block(mapStationNames);
    mapStationNames->setChecked(options.showStationNames);
  }
  if(plotHandlers.size() > 1)
  {
    if(auto *phasePlot = dynamic_cast<PhasePlot*>(plotHandlers[1].get()))
      phasePlot->set_phase_wrap(options.phaseWrap);
  }
  // Keep reading the legacy tipperArrows archive field; the main view always uses arrows.
  mapHandler->set_station_names_visible(options.showStationNames);
  for(unsigned i = 0; i < plotHandlers.size(); ++i)
    plotHandlers[i]->set_component_visibility(options.componentVisible[i]);
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::on_actionLoad_EDI_triggered()
{
  const auto files = QFileDialog::getOpenFileNames(this, tr("Select EDI files"), lastDirectory, tr("EDI files (*.edi)"));
  if(files.isEmpty()) return;
  QString surveyName;
  if(!mtSurvey) {
    bool accepted = false;
    surveyName = QInputDialog::getText(this, tr("New survey"), tr("Survey name:"), QLineEdit::Normal, "Survey", &accepted);
    if(!accepted) return;
  }
  rememberDirectory(files.first());
  try {
    // Parse and review a copy: cancellation or a failed file leaves the open project intact.
    auto candidate = mtSurvey ? std::make_shared<MTSurveyData>(*mtSurvey) : std::make_shared<MTSurveyData>(surveyName.toStdString());
    const auto oldNames = candidate->get_stations_names();
    const std::set<std::string> existing(oldNames.begin(), oldNames.end());
    std::vector<std::string> paths;
    for(const auto &file: files) paths.push_back(file.toStdString());
    const auto duplicates = candidate->load_from_edi(paths);
    std::set<std::string> incoming;
    for(const auto &name: candidate->get_stations_names()) if(!existing.count(name)) incoming.insert(name);
    if(!reviewEDIImport(this, *candidate, incoming, duplicates)) return;
    if(mtSurvey) { mtSurvey->restore(*candidate); history->record(tr("Import EDI stations")); }
    else { mtSurvey = std::move(candidate); history->reset(mtSurvey); }
    projectFile.clear();
    setWindowTitle(tr("MT Survey [%1]").arg(QString::fromStdString(mtSurvey->get_survey_name())));
    createStationsList();
    updateMap();
    refreshAnalysisWindows();
    if(ui->stationList->count()) ui->stationList->setCurrentRow(0);
  } catch(const std::exception &error) {
    QMessageBox::critical(this, tr("EDI import"), QString::fromUtf8(error.what()));
  }
}

void MainWindow::createStationsList()
{
  const QSignalBlocker blocker(ui->stationList);
  ui->stationList->clear();

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
  if(!mtSurvey) return;
  SurveyCoordinates view;
  if(mapCoordinateSwitch->currentIndex() == 1 && mtSurvey->n_stations())
    view = mtSurvey->coordinates().utm ? mtSurvey->coordinates() :
      SurveyCoordinates::suggested(mtSurvey->geographic_locations());
  const auto locations = view.transform(mtSurvey->geographic_locations());
  mapCoordinates = view;
  {
    const QSignalBlocker blocker(mapCenterUtm);
    mapCenterUtm->setEnabled(view.utm);
    mapCenterUtm->setChecked(mtSurvey->coordinates().centered);
  }
  mapBackground->setCoordinates(view);
  mapBackground->setTopography(mtSurvey->topography());
  const auto &coordinates = view.utm ? view : mtSurvey->coordinates();
  mapHandler->set_coordinate_labels(view.utm, view.centered);
  coordinateInfoLabel->setText(coordinates.utm ?
    tr("WGS84 / UTM %1%2\nOrigin E: %3 m\nOrigin N: %4 m")
      .arg(coordinates.zone).arg(coordinates.north ? "N" : "S")
      .arg(coordinates.origin_easting, 0, 'f', 3).arg(coordinates.origin_northing, 0, 'f', 3) :
    tr("Latitude / longitude (degrees)"));
  coordinateInfoLabel->setToolTip(QString::fromStdString(coordinates.description()));
  mapHandler->set_data(locations,
                       mtSurvey->get_stations_names());
  const QSignalBlocker blocker(ui->stationList);
  const auto names = mtSurvey->get_stations_names();
  std::map<std::string, std::array<double, 3>> byName;
  for(unsigned i = 0; i < names.size(); ++i) byName.emplace(names[i], locations[i]);
  for(int i = 0; i < ui->stationList->count(); ++i) {
    auto *item = ui->stationList->item(i);
    const auto &p = byName.at(item->text().toStdString());
    const auto text = !std::isfinite(p[0]) || !std::isfinite(p[1]) ? tr("Coordinates unavailable") :
      view.utm ? tr("E %1  N %2 m").arg(p[1], 0, 'f', 1).arg(p[0], 0, 'f', 1) :
                 tr("Lat %1°  Lon %2°").arg(p[0], 0, 'f', 4).arg(p[1], 0, 'f', 4);
    item->setData(StationListDelegate::Coordinates, text);
    item->setData(StationListDelegate::CoordinateReference, QString::fromStdString(view.description()));
  }
  updateStationDetails();
  on_stationList_itemSelectionChanged();
}

void MainWindow::updatePlots(bool rescaleAxes)
{
  // Designer-connected actions can fire while the plots are being constructed.
  if(!rhoPlusOverlay) return;
  updateStationDetails();
  updateActionAvailability();
  rhoPlusOverlay->hide();
  QListWidgetItem* item = ui->stationList->currentItem();
  if(stationHeading) {
    stationHeading->setText(item ? item->text() : tr("No station selected"));
    stationHeading->setToolTip(stationHeading->text());
  }
  if(autoMaskWorkflow) autoMaskWorkflow->setData(mtSurvey, item ? item->text() : QString());
  if(item == nullptr) {
    for(auto &plot: plotHandlers) plot->clear_observed_data();
    rhoPlusOverlay->update(nullptr,rhoPlusAction->isChecked(),ui->actionPhase_wrap->isChecked(),
      plotHandlers[0]->component_visibility(),plotHandlers[1]->component_visibility(),rhoPlusOptions);
    return;
  }

  MTStationData& data = mtSurvey->get_station_data(item->text().toStdString());
  for(auto &ph: plotHandlers) {
    ph->clear_predicted_data();
    ph->set_observed_data(data, rescaleAxes);
  }
  ui->label->setText(tr("Stations (%1):").arg(mtSurvey->n_stations()));

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

        const double rms = data.rms(resp_data);
        ui->label->setText(std::isfinite(rms) ? tr("RMS = %1").arg(rms) :
                                             tr("RMS: no matching data"));
      }
      else
        ui->label->setText(tr("No response for this station"));
    }
  }

  ui->statusBar->showMessage("# frequencies: " + QString("%1").arg(data.frequencies().size()) +
                             "\t\t\t# path: " + data.file_path().c_str());
  rhoPlusOverlay->update(&data,rhoPlusAction->isChecked(),ui->actionPhase_wrap->isChecked(),
    plotHandlers[0]->component_visibility(),plotHandlers[1]->component_visibility(),rhoPlusOptions);
  if(rescaleAxes && periodRange.automatic && !data.frequencies().empty()) {
    const auto bounds = std::minmax_element(data.frequencies().begin(), data.frequencies().end());
    if(*bounds.first > 0. && std::isfinite(*bounds.second))
      plot11->xAxis->setRange(.5 / *bounds.second, 2. / *bounds.first);
  }
  if(!periodRange.automatic) plot11->xAxis->setRange(periodRange.lower, periodRange.upper);
  if(viewControls) viewControls->refreshRange();
}

void MainWindow::stationSelected(const QCPDataSelection &selection)
{
  if(!selection.isEmpty())
  {
    auto range = selection.dataRange();

    if(range.length() > 0)
    {
      const QString name = QString::fromStdString(mapHandler->station_name(range.begin()));

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

  for(auto target: MTStationData::mask_edit_types(type, linkTensorMasksAction->isChecked()))
    mtSurvey->get_station_data(station_name).mask_type(target, on);
  history->record(tr("Edit component masks"));

  updatePlots();
  refreshAnalysisWindows();
}

void MainWindow::deleteStations()
{
  const auto selected = ui->stationList->selectedItems();
  if(selected.isEmpty() || !mtSurvey) return;
  for(const auto *item : selected)
    mtSurvey->remove_station(item->text().toStdString());
  history->record(tr("Delete %n station(s)", nullptr, selected.size()));
  restoreSurveyView();
}

void MainWindow::restoreSurveyView()
{
  const auto *item = ui->stationList->currentItem();
  const QString selected = item ? item->text() : QString();
  std::set<QString> selection;
  for(const auto *entry : ui->stationList->selectedItems()) selection.insert(entry->text());
  {
    const QSignalBlocker blocker(ui->stationList);
    createStationsList();
    const auto matches = ui->stationList->findItems(selected, Qt::MatchExactly);
    if(!matches.isEmpty()) ui->stationList->setCurrentItem(matches.front(), QItemSelectionModel::NoUpdate);
    else if(ui->stationList->count()) ui->stationList->setCurrentRow(0, QItemSelectionModel::NoUpdate);
    for(int i = 0; i < ui->stationList->count(); ++i)
      ui->stationList->item(i)->setSelected(selection.count(ui->stationList->item(i)->text()) || (selected.isEmpty() && selection.empty() && i == 0));
  }
  updateMap(); updatePlots(false); refreshAnalysisWindows();
}

void MainWindow::renameStation()
{
  QListWidgetItem* item = ui->stationList->currentItem();

  if(!item || !mtSurvey) return;
  bool ok;
  QString newName = QInputDialog::getText(this, tr("Name"),
                                          tr("New name:"), QLineEdit::Normal,
                                          item->text(), &ok);
  if (ok && !newName.isEmpty())
  {
    try {
      mtSurvey->rename_station(item->text().toStdString(), newName.toStdString());
      { const QSignalBlocker block(ui->stationList); item->setText(newName); }
      history->record(tr("Rename station"));
      updateMap();
      updatePlots(false);
      refreshAnalysisWindows();
    } catch(const std::exception &e) {
      QMessageBox::warning(this,tr("Rename station"),QString::fromUtf8(e.what()));
    }
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

    if (projectFile.isEmpty())
      return;

    if (!projectFile.endsWith(".mtd"))
      projectFile += ".mtd";

    rememberDirectory(projectFile);
  }

  std::ofstream ofs(projectFile.toStdString(), std::ios::binary);

  if(ofs.is_open())
  {
    boost::archive::binary_oarchive oa(ofs);
    oa << projectFileMagic;
    oa << projectFileVersion;
    oa << mtSurvey;
    oa << mtResponses;
    oa << currentPlotOptions();
  }

  ofs.close();
}

void MainWindow::on_actionLoad_project_triggered()
{
  const QString path = QFileDialog::getOpenFileName(
                            this,
                            "Open project file",
                            lastDirectory,
                            "Project file (*.mtd)");

  if(path.isEmpty())
    return;

  rememberDirectory(path);

  bool loaded = false;
  std::shared_ptr<MTSurveyData> loadedSurvey;
  std::map<std::string,MTSurveyData> loadedResponses;

  try
  {
    std::ifstream ifs(path.toStdString(), std::ios::binary);
    if(ifs.is_open())
    {
      boost::archive::binary_iarchive ia(ifs);
      std::uint32_t magic = 0;
      std::uint32_t version = 0;
      ia >> magic;
      ia >> version;

      if(magic != projectFileMagic)
        throw std::runtime_error("Legacy project format");

      ia >> loadedSurvey;
      ia >> loadedResponses;

      PlotOptions plotOptions;
      if(version == 2)
      {
        PlotOptionsV2 legacyPlotOptions;
        ia >> legacyPlotOptions;
        plotOptions.phaseWrap = legacyPlotOptions.phaseWrap;
        plotOptions.axes = legacyPlotOptions.axes;
      }
      else if(version == 3)
      {
        PlotOptionsV3 legacyPlotOptions;
        ia >> legacyPlotOptions;
        plotOptions.phaseWrap = legacyPlotOptions.phaseWrap;
        plotOptions.showStationNames = legacyPlotOptions.showStationNames;
        plotOptions.axes = legacyPlotOptions.axes;
      }
      else if(version == 4)
      {
        PlotOptionsV4 legacyPlotOptions;
        ia >> legacyPlotOptions;
        plotOptions.phaseWrap = legacyPlotOptions.phaseWrap;
        plotOptions.showStationNames = legacyPlotOptions.showStationNames;
        plotOptions.tipperArrows = legacyPlotOptions.tipperArrows;
        plotOptions.axes = legacyPlotOptions.axes;
      }
      else if(version == 5)
      {
        PlotOptionsV5 legacyPlotOptions;
        ia >> legacyPlotOptions;
        plotOptions.phaseWrap = legacyPlotOptions.phaseWrap;
        plotOptions.showStationNames = legacyPlotOptions.showStationNames;
        plotOptions.tipperArrows = legacyPlotOptions.tipperArrows;
        plotOptions.axes = legacyPlotOptions.axes;
        plotOptions.componentVisible = {{legacyPlotOptions.impedanceVisible, legacyPlotOptions.impedanceVisible,
                                        legacyPlotOptions.tipperVisible, legacyPlotOptions.impedanceVisible}};
      }
      else if(version == 6)
      {
        PlotOptionsV6 legacyPlotOptions;
        ia >> legacyPlotOptions;
        static_cast<PlotOptionsV6 &>(plotOptions) = legacyPlotOptions;
      }
      else if(version == 7)
      {
        PlotOptionsV7 legacyPlotOptions;
        ia >> legacyPlotOptions;
        static_cast<PlotOptionsV7 &>(plotOptions) = legacyPlotOptions;
      }
      else if(version == 8)
      {
        PlotOptionsV8 legacyPlotOptions;
        ia >> legacyPlotOptions;
        static_cast<PlotOptionsV8 &>(plotOptions) = legacyPlotOptions;
      }
      else if(version >= 9)
      {
        ia >> plotOptions;
      }
      if(version < 9) for(unsigned i = 0; i < 2 && i < plotOptions.axes.size(); ++i) {
        const auto &axis = plotOptions.axes[i];
        plotOptions.diagonals[i] = {axis.autoscale, axis.lower, axis.upper};
      }
      applyPlotOptions(plotOptions);
      loaded = true;
    }
  }
  catch(const std::exception&)
  {
    // No valid marker/version was found, so retry as the original project
    // layout: mtSurvey followed by mtResponses, with no plot options.
    try
    {
      std::ifstream ifs(path.toStdString(), std::ios::binary);
      if(ifs.is_open())
      {
        boost::archive::binary_iarchive ia(ifs);
        ia >> loadedSurvey;
        ia >> loadedResponses;

        PlotOptions defaultPlotOptions;
        for(const auto &plot: plotHandlers)
        {
          const QCPRange range = plot->y_axis_range();
          PlotAxisOptions axisOptions;
          axisOptions.autoscale = true;
          axisOptions.lower = range.lower;
          axisOptions.upper = range.upper;
          defaultPlotOptions.axes.push_back(axisOptions);
        }
        applyPlotOptions(defaultPlotOptions);
        loaded = true;
      }
    }
    catch(const std::exception&)
    {
      loaded = false;
    }
  }

  if(!loaded || !loadedSurvey)
  {
    QMessageBox::critical(this, tr("Load project"),
                          tr("Failed to load project file."));
    return;
  }

  mtSurvey=std::move(loadedSurvey);
  mtResponses=std::move(loadedResponses);
  projectFile=path;
  history->reset(mtSurvey);
  createResponsesList();
  restoreSurveyView();
}

void MainWindow::on_stationList_currentRowChanged(int /*currentRow*/)
{
  updatePlots();
}

void MainWindow::on_stationList_itemChanged(QListWidgetItem *item)
{
  const QSignalBlocker blocker(ui->stationList);
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

  history->record(tr("Enable / disable stations"));
  updatePlots();
  refreshAnalysisWindows();
}

void MainWindow::on_actionDecimate_triggered()
{
  if(mtSurvey == nullptr)
    return;

  mtSurvey->decimate();
  history->record(tr("Decimate periods"));
  updatePlots();
  refreshAnalysisWindows();
}

void MainWindow::on_stationList_itemSelectionChanged()
{
  updateActionAvailability();
  if(!mtSurvey) return;
  std::vector<std::array<double, 3>> geographic;

  auto selection = ui->stationList->selectedItems();
  for(auto sitem: selection)
    geographic.push_back(mtSurvey->get_station_data(sitem->text().toStdString()).position());

  mapHandler->set_selected_points(mapCoordinates.transform(geographic));
}

void MainWindow::on_actionExport_in_GoFEM_triggered()
{
  if(mtSurvey == nullptr)
    return;

  ExportGOFEMDialog dialog(*mtSurvey,this);
  if(dialog.exec()!=QDialog::Accepted) return;
  try {
    const auto path=saveSurveyExport(this,dialog.exportPreview().result(),SurveyExport::Format::GoFEM,lastDirectory);
    if(!path.isEmpty())rememberDirectory(path);
  } catch(const std::exception &e) { QMessageBox::warning(this,tr("Export data"),QString::fromUtf8(e.what())); }
}

void MainWindow::on_actionExport_native_MT_triggered()
{
  if(mtSurvey == nullptr) return;
  ExportNativeMTDialog dialog(*mtSurvey, lastDirectory, this);
  if(dialog.exec() == QDialog::Accepted) rememberDirectory(dialog.exportedFile());
  updateMap();
}

void MainWindow::on_actionConvert_UTM_triggered()
{
  if(mtSurvey == nullptr) return;
  try {
    SurveyCoordinatesDialog dialog(*mtSurvey, this);
    if(dialog.exec() == QDialog::Accepted) {
      history->record(tr("UTM projection"));
      const QSignalBlocker blocker(mapCoordinateSwitch);
      mapCoordinateSwitch->setCurrentIndex(1);
      updateMap();
      refreshAnalysisWindows();
    }
  } catch(const std::exception &e) {
    QMessageBox::warning(this, tr("Coordinate conversion"), QString::fromUtf8(e.what()));
  }
}

void MainWindow::on_actionGeographic_coordinates_triggered()
{
  if(mtSurvey == nullptr) return;
  mapCoordinateSwitch->setCurrentIndex(0);
}

void MainWindow::on_actionSet_error_floor_triggered()
{
  bool ok;
  double error_floor = QInputDialog::getDouble(this, "Specify error floor",
                                               "Error floor (in %):",
                                               5, 0, 100, 2, &ok);
  if(ok && mtSurvey != nullptr)
    mtSurvey->set_error_floor(error_floor / 100.);
  history->record(tr("Error floor"));

  updatePlots();
  refreshAnalysisWindows();
}

void MainWindow::on_actionShow_error_bars_toggled(bool on)
{
  if(plotHandlers.empty()) return;
  for(auto &plot: plotHandlers)
    plot->set_error_bars_visible(on);

  for(auto *plot: {plot11, plot12, plot21, plot22})
    plot->replot(QCustomPlot::rpQueuedReplot);
}

void MainWindow::on_actionPhase_wrap_toggled(bool on)
{
  if(plotHandlers.size() > 1)
  {
    if(auto *phasePlot = dynamic_cast<PhasePlot*>(plotHandlers[1].get()))
      phasePlot->set_phase_wrap(on);
  }

  updatePlots();
}

void MainWindow::on_actionShow_station_names_toggled(bool on)
{
  const QSignalBlocker block(mapStationNames);
  mapStationNames->setChecked(on);
  mapHandler->set_station_names_visible(on);
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
  listContextMenu = new QMenu(this);
  QMenu* mask_menu = listContextMenu->addMenu(tr("Mask"));
  listContextMenu->addAction(findChild<QAction *>("actionGroup_masking"));
  listContextMenu->addAction(findChild<QAction *>("actionAuto_masking"));
  auto *remove = listContextMenu->addAction(tr("Delete"), this, &MainWindow::deleteStations);
  remove->setObjectName("actionDelete_stations");
  remove->setIcon(commandIcon("delete"));
  remove->setToolTip(tr("Delete all selected stations. Undo restores them together."));
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
        lastDirectory + "/" + sname + ".png",
        "Image file (*.png)");

  if(exportFile.length() == 0)
    return;

  rememberDirectory(exportFile);

  tfView->grab().save(exportFile, "PNG", 100);
}

void MainWindow::on_actionPlot_axis_ranges_triggered()
{
  if(ribbon) ribbon->showPage(6);
  if(viewControls) viewControls->refreshRange();
}

void MainWindow::on_actionLoad_responses_triggered()
{
  loadResponseFiles(QFileDialog::getOpenFileNames(
    this, tr("Open computed responses"), lastDirectory,
    tr("Response files (*.txt *.data *.dat *.gofem);;All files (*)")));
}

void MainWindow::loadResponseFiles(const QStringList &files)
{
  if(files.isEmpty()) return;
  rememberDirectory(files.front());
  QStringList errors;
  QListWidgetItem *selected = nullptr;
  auto sorted = files;
  sorted.sort();
  {
    const QSignalBlocker blocker(ui->responsesList);
    for(const auto &file: sorted) {
      const QString path = QFileInfo(file).absoluteFilePath();
      try {
        MTSurveyData responses;
        responses.load_responses(path.toStdString());
        mtResponses[path.toStdString()] = std::move(responses);

        QListWidgetItem *item = nullptr;
        for(int i = 0; i < ui->responsesList->count(); ++i)
          if(ui->responsesList->item(i)->toolTip() == path)
            item = ui->responsesList->item(i);
        if(!item) {
          item = new QListWidgetItem(FileLabels::fileName(path), ui->responsesList);
          item->setToolTip(path);
        }
        selected = item;
      } catch(const std::exception &e) {
        errors.push_back(path + "\n" + QString::fromUtf8(e.what()));
      }
    }
    if(selected) ui->responsesList->setCurrentItem(selected);
  }
  if(selected) ui->responsesToggle->setChecked(true);
  updatePlots();
  if(!errors.isEmpty())
    QMessageBox::warning(this, tr("Load responses"), errors.join("\n\n"));
  refreshAnalysisWindows();
}

void MainWindow::on_actionFit_statistics_triggered()
{
  if(!fitStatistics)
    fitStatistics = new FitStatisticsWindow(this, [this] { refreshAnalysisWindows(); });
  refreshAnalysisWindows();
  fitStatistics->show();
  fitStatistics->raise();
  fitStatistics->activateWindow();
}

void MainWindow::refreshAnalysisWindows()
{
  updateActionAvailability();
  if(fitStatistics) fitStatistics->setData(mtSurvey, mtResponses);
  if(periodMaps) periodMaps->setData(mtSurvey, mtResponses);
  if(profile) profile->setData(mtSurvey, mtResponses);
  if(periodLayout) periodLayout->setData(mtSurvey);
  if(impedanceValidity) impedanceValidity->setData(mtSurvey);
  if(surveyAverage) surveyAverage->setData(mtSurvey);
  if(stationComparison) stationComparison->setData(mtSurvey);
  if(polarity) polarity->setData(mtSurvey);
}

void MainWindow::openPeriodLayout()
{
  if(!periodLayout) periodLayout = new PeriodLayoutWindow(this, [this](std::shared_ptr<MTSurveyData> data) { openSurveyCopy(std::move(data)); });
  periodLayout->setData(mtSurvey);
  periodLayout->show(); periodLayout->raise(); periodLayout->activateWindow();
}

void MainWindow::openSurveyCopy(std::shared_ptr<MTSurveyData> data)
{
  auto *window = new MainWindow;
  window->setAttribute(Qt::WA_DeleteOnClose);
  window->setObjectName("resampledSurveyWindow");
  window->mtSurvey = std::move(data);
  window->history->reset(window->mtSurvey);
  window->mtResponses = mtResponses;
  window->lastDirectory = lastDirectory;
  window->projectFile.clear();
  window->applyPlotOptions(currentPlotOptions());
  window->createStationsList(); window->createResponsesList(); window->updateMap();
  if(window->ui->stationList->count()) window->ui->stationList->setCurrentRow(0);
  const auto name = QString::fromStdString(window->mtSurvey->get_survey_name()).trimmed();
  window->setWindowTitle(tr("EDITools — %1").arg(name.isEmpty() ? tr("Survey copy") : name));
  window->show();
}

void MainWindow::on_actionPeriod_maps_triggered()
{
  const bool firstOpen = !periodMaps;
  if(firstOpen) periodMaps = new PeriodMapWindow(this, [this] { history->record(tr("Map masks")); updatePlots(); refreshAnalysisWindows(); });
  periodMaps->setLinkTensorMasks(linkTensorMasksAction->isChecked());
  periodMaps->setData(mtSurvey, mtResponses);
  if(firstOpen) {
    auto settings = periodMaps->displaySettings();
    settings.mapLayers = mapBackground->layers();
    settings.names = mapStationNames->isChecked();
    periodMaps->setDisplaySettings(settings);
  }
  periodMaps->show();
  periodMaps->raise();
  periodMaps->activateWindow();
}

void MainWindow::createResponsesList()
{
  const QSignalBlocker blocker(ui->responsesList);
  ui->responsesList->clear();
  for(const auto &response: mtResponses) {
    const QString path = QString::fromStdString(response.first);
    auto *item = new QListWidgetItem(FileLabels::fileName(path), ui->responsesList);
    item->setToolTip(path);
  }
  if(ui->responsesList->count())
    ui->responsesList->setCurrentRow(ui->responsesList->count() - 1);
  ui->responsesToggle->setChecked(ui->responsesList->count() > 0);
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
