#include "ui/masking/AutoMaskDialog.h"
#include "ui/masking/AutoMaskWorkflow.h"
#include "ui/analysis/FitStatisticsWindow.h"
#include "ui/analysis/ImpedanceAnalysisWindows.h"
#include "ui/analysis/PeriodLayoutWindow.h"
#include "ui/analysis/PeriodMapWindow.h"
#include "ui/analysis/PolarityWindow.h"
#include "ui/analysis/ProfileWindow.h"
#include "ui/analysis/StationComparisonWindow.h"
#include "ui/maps/MapBackground.h"
#include "ui/workspace/SurveyHistory.h"
#include "ui/workspace/WorkspaceTabs.h"
#include "ui/workspace/WorkspaceView.h"
#include "app/MainWindow.h"
#include "ui_MainWindow.h"
#include <QAbstractButton>
#include <QFileDialog>
#include <QJsonArray>
#include <QMessageBox>
#include <QRegularExpression>
#include <QStackedWidget>
#include <algorithm>
#include <iterator>

namespace {
struct WorkspaceKind {
  const char *id;
  const char *title;
};
const WorkspaceKind kinds[]{
    {"profile", QT_TR_NOOP("2D profile")},       {"map", QT_TR_NOOP("Period map")},
    {"fit", QT_TR_NOOP("Fit statistics")},       {"comparison", QT_TR_NOOP("Station comparison")},
    {"validity", QT_TR_NOOP("Rho+ validity")},   {"average", QT_TR_NOOP("Survey average")},
    {"polarity", QT_TR_NOOP("Polarity review")}, {"periods", QT_TR_NOOP("Period layout")},
    {"automask", QT_TR_NOOP("Auto-mask review")}};
} // namespace

void MainWindow::setupWorkspace() {
  sharedMapHost = new QStackedWidget;
  sharedMapHost->setObjectName("sharedMapHost");
  sharedMapHost->setSizePolicy(ui->mapPlot->sizePolicy());
  const int mapIndex = ui->splitter_3->indexOf(ui->mapPlot);
  ui->splitter_3->replaceWidget(mapIndex, sharedMapHost);
  sharedMapHost->addWidget(ui->mapPlot);
  ui->splitter_3->setStretchFactor(mapIndex, 3);
  workspaces = new WorkspaceTabs;
  ui->splitter_5->replaceWidget(ui->splitter_5->indexOf(tfView), workspaces);
  ui->splitter_5->setStretchFactor(0, 0);
  ui->splitter_5->setStretchFactor(1, 1);
  ui->splitter_5->setCollapsible(0, false);
  workspaces->addView(tfView, "tf", tr("TF view"), "tf");
  connect(workspaces, &QTabWidget::currentChanged, this, [this] { activateWorkspace(); });
  connect(workspaces, &WorkspaceTabs::closing, this, [this](QWidget *view) {
    if (auto *analysis = dynamic_cast<WorkspaceView *>(view))
      analysis->setSharedMapActive(false, mapCoordinates);
    if (!restoringWorkspace) {
      closedViews.append(workspaceSnapshot(view));
      if (closedViews.size() > 20)
        closedViews.removeFirst();
    }
    if (auto *v = dynamic_cast<WorkspaceView *>(view))
      if (v->centralMaps().contains(sharedMapHost->currentWidget()))
        sharedMapHost->setCurrentWidget(ui->mapPlot);
    if (ribbon)
      ribbon->releaseContext(workspaces->id(view));
  });
  connect(workspaces, &WorkspaceTabs::titleChanged, this,
          [](QWidget *view, const QString &title) { view->setProperty("viewTitle", title); });
  connect(workspaces, &WorkspaceTabs::duplicateRequested, this, [this](QWidget *view) {
    auto saved = workspaceSnapshot(view);
    saved.remove("id");
    saved["title"] = saved["title"].toString() + tr(" copy");
    openWorkspace(workspaces->type(view), true, saved);
  });
  connect(workspaces, &WorkspaceTabs::reopenRequested, this, [this] {
    if (closedViews.isEmpty())
      return;
    const auto saved = closedViews.takeLast();
    openWorkspace(saved["type"].toString(), true, saved);
  });
  auto *menu = menuBar()->addMenu(tr("Views"));
  auto *newView = menu->addAction(tr("New view"));
  newView->setObjectName("actionNew_view");
  auto *choices = new QMenu(menu);
  newView->setMenu(choices);
  for (const auto &kind : kinds) {
    const QString type = QString::fromLatin1(kind.id);
    connect(choices->addAction(tr(kind.title)), &QAction::triggered, this,
            [this, type] { openWorkspace(type, true); });
  }
  auto *duplicate = menu->addAction(tr("Duplicate view"));
  duplicate->setObjectName("actionDuplicate_view");
  connect(duplicate, &QAction::triggered, this,
          [this] { emit workspaces->duplicateRequested(workspaces->currentWidget()); });
  auto *rename = menu->addAction(tr("Rename tab"));
  rename->setObjectName("actionRename_view");
  connect(rename, &QAction::triggered, workspaces, &WorkspaceTabs::renameCurrent);
  auto *close = menu->addAction(tr("Close tab"));
  close->setObjectName("actionClose_view");
  connect(close, &QAction::triggered, workspaces, &WorkspaceTabs::closeCurrent);
  auto *reopen = menu->addAction(tr("Reopen closed tab"));
  reopen->setObjectName("actionReopen_view");
  connect(reopen, &QAction::triggered, this, [this] { emit workspaces->reopenRequested(); });
  auto *navigator = menu->addAction(tr("Station navigator"));
  navigator->setObjectName("actionStation_navigator");
  navigator->setCheckable(true);
  navigator->setChecked(true);
  connect(navigator, &QAction::toggled, ui->splitter_3, &QWidget::setVisible);
  for (const auto &format : {QString("png"), QString("pdf"), QString("csv")}) {
    auto *action = menu->addAction(tr("Active view %1…").arg(format.toUpper()));
    action->setObjectName("actionView_" + format);
    connect(action, &QAction::triggered, this, [this, format] { exportWorkspace(format); });
  }
}

WorkspaceView *MainWindow::createWorkspaceView(const QString &type) {
  const auto changed = [this] {
    updatePlots(false);
    refreshAnalysisWindows();
  };
  const auto copy = [this](std::shared_ptr<MTSurveyData> data) { openSurveyCopy(std::move(data)); };
  if (type == "profile") {
    auto *view = new ProfileWindow(this, copy, ui->mapPlot);
    view->setMapLayers(mapBackground->layers());
    return view;
  }
  if (type == "map") {
    auto *view = new PeriodMapWindow(
        this,
        [this, changed] {
          history->record(tr("Map masks"));
          changed();
        },
        true);
    auto settings = view->displaySettings();
    settings.mapLayers = mapBackground->layers();
    settings.names = mapStationNames->isChecked();
    view->setDisplaySettings(settings);
    view->setLinkTensorMasks(linkTensorMasksAction->isChecked());
    lastPeriodMap = view;
    return view;
  }
  if (type == "fit")
    return new FitStatisticsWindow(
        this, [this] { refreshAnalysisWindows(); }, true);
  if (type == "comparison") {
    auto *view = new StationComparisonWindow(this, ui->mapPlot);
    view->setData(mtSurvey, {currentStationName()});
    return view;
  }
  if (type == "validity") {
    auto *view = new ImpedanceValidityWindow(
        this,
        [this](const ImpedanceAnalysis::RhoPlusOptions &options) {
          rhoPlusOptions = options;
          updatePlots(false);
        },
        changed);
    view->setLinkTensorMasks(linkTensorMasksAction->isChecked());
    view->setFitOptions(rhoPlusOptions);
    view->setHistory(history.get());
    view->setData(mtSurvey, currentStationName());
    lastImpedanceValidity = view;
    return view;
  }
  if (type == "average")
    return new SurveyAverageWindow(this);
  if (type == "polarity") {
    auto *view = new PolarityWindow(this, *history, changed);
    view->setData(mtSurvey, currentStationName());
    return view;
  }
  if (type == "periods")
    return new PeriodLayoutWindow(this, copy);
  if (type == "automask") {
    std::set<std::string> selection;
    for (const auto *item : ui->stationList->selectedItems())
      selection.insert(item->text(0).toStdString());
    if (selection.empty() && ui->stationList->currentItem())
      selection.insert(currentStationName().toStdString());
    return new AutoMaskReview(this, mtSurvey, selection, linkTensorMasksAction->isChecked(),
                              autoMaskWorkflow->currentOptions(), autoMaskWorkflow->pairDiagonals(),
                              [this, changed] {
                                history->record(tr("Auto masking"));
                                changed();
                              });
  }
  return nullptr;
}

QWidget *MainWindow::openWorkspace(const QString &type, bool create, const QJsonObject &saved) {
  if (!workspaces || !mtSurvey ||
      !std::any_of(std::begin(kinds), std::end(kinds),
                   [&](const WorkspaceKind &kind) { return type == kind.id; }))
    return nullptr;
  if (!create)
    for (int i = workspaces->count() - 1; i > 0; --i)
      if (workspaces->type(workspaces->widget(i)) == type) {
        workspaces->setCurrentIndex(i);
        return workspaces->widget(i);
      }
  auto *view = createWorkspaceView(type);
  if (!view)
    return nullptr;
  for (auto *map : view->centralMaps())
    sharedMapHost->addWidget(map);
  connect(view, &WorkspaceView::centralMapChanged, this, [this, view] {
    if (workspaces->currentWidget() == view)
      syncCentralMap();
  });
  view->setEmbedded();
  view->markSurveyDirty();
  refreshWorkspace(view);
  QString title = saved["title"].toString();
  if (title.isEmpty()) {
    title = view->windowTitle();
    int number = 1;
    for (int i = 1; i < workspaces->count(); ++i)
      if (workspaces->type(workspaces->widget(i)) == type)
        ++number;
    if (number > 1)
      title += QString(" %1").arg(number);
  }
  view->setProperty("viewTitle", title);
  if (saved.contains("state")) {
    view->restoreViewState(saved["state"].toObject());
    view->markSurveyDirty();
  }
  connect(view, &WorkspaceView::toolsChanged, this, [this, view] {
    if (view == activeAnalysis())
      updateWorkspaceActions();
  });
  connect(view, &WorkspaceView::stationActivated, this, [this](const QString &name) {
    const auto items = ui->stationList->findItems(name, Qt::MatchExactly);
    if (!items.isEmpty())
      ui->stationList->setCurrentItem(items.front(), 0, QItemSelectionModel::NoUpdate);
  });
  workspaces->addView(view, type, title, saved["id"].toString());
  if (ribbon && !restoringWorkspace)
    ribbon->setContext(workspaces->id(view), view->toolPages(), true);
  return view;
}

QString MainWindow::currentStationName() const {
  const auto *item = ui->stationList->currentItem();
  return item ? item->text(0) : QString();
}

WorkspaceView *MainWindow::activeAnalysis() const {
  return workspaces ? dynamic_cast<WorkspaceView *>(workspaces->currentWidget()) : nullptr;
}

void MainWindow::refreshWorkspace(QWidget *widget) {
  if (auto *view = dynamic_cast<WorkspaceView *>(widget)) {
    view->refreshSurvey(mtSurvey, mtResponses);
    if (view == activeAnalysis())
      syncCentralMap();
  }
}

void MainWindow::activateWorkspace() {
  if (!ribbon || !workspaces->currentWidget())
    return;
  auto *view = workspaces->currentWidget();
  refreshWorkspace(view);
  // Release the previous overlay before the next view captures map interactions.
  for (int i = 0; i < workspaces->count(); ++i)
    if (auto *analysis = dynamic_cast<WorkspaceView *>(workspaces->widget(i)))
      if (analysis != view)
        analysis->setSharedMapActive(false, mapCoordinates);
  if (auto *analysis = activeAnalysis()) {
    analysis->setSharedMapActive(true, mapCoordinates);
    analysis->selectStation(currentStationName());
  }
  syncCentralMap();
  auto *analysis = activeAnalysis();
  ribbon->setContext(workspaces->id(view), analysis ? analysis->toolPages() : tfTools);
  updateWorkspaceActions();
}

void MainWindow::updateWorkspaceActions() {
  auto *analysis = activeAnalysis();
  for (const auto &name : {"actionClose_view", "actionDuplicate_view"})
    findChild<QAction *>(name)->setEnabled(analysis);
  findChild<QAction *>("actionReopen_view")->setEnabled(!closedViews.isEmpty());
  for (const auto &format : {QString("pdf"), QString("csv")}) {
    auto *button = analysis ? analysis->exportButton(format) : nullptr;
    findChild<QAction *>("actionView_" + format)->setEnabled(button && button->isEnabled());
  }
}

void MainWindow::syncCentralMap() {
  if (!sharedMapHost || !workspaces)
    return;
  auto *view = activeAnalysis();
  auto *map = view ? view->centralMap() : nullptr;
  sharedMapHost->setCurrentWidget(map ? map : ui->mapPlot);
  if (auto *controls = findChild<QWidget *>("surveyMapControls"))
    controls->setVisible(!map);
  if (mtSurvey)
    updateStationCoordinates(map ? view->centralMapCoordinates() : mapCoordinates);
}

QJsonObject MainWindow::workspaceSnapshot(QWidget *view) const {
  QJsonObject state;
  if (auto *analysis = dynamic_cast<WorkspaceView *>(view))
    state = analysis->saveViewState();
  else
    state = {{"station", currentStationName()},
             {"response",
              ui->responsesList->currentItem() ? ui->responsesList->currentItem()->toolTip() : QString()},
             {"autoMask", autoMaskWorkflow->settingsView()->saveViewState()},
             {"coordinates", mapCoordinateSwitch->currentIndex()},
             {"responses", ui->responsesToggle->isChecked()},
             {"navigator", findChild<QAction *>("actionStation_navigator")->isChecked()},
             {"splitter", QString::fromLatin1(ui->splitter_5->saveState().toBase64())},
             {"plots", WorkspaceView::savePlotRanges(tfView)}};
  return {{"id", workspaces->id(view)},
          {"type", workspaces->type(view)},
          {"title", workspaces->tabText(workspaces->indexOf(view))},
          {"state", state}};
}
QJsonObject MainWindow::saveWorkspace() const {
  QJsonArray views;
  for (int i = 0; i < workspaces->count(); ++i)
    views.append(workspaceSnapshot(workspaces->widget(i)));
  return {{"schema", 1}, {"views", views}, {"active", workspaces->id(workspaces->currentWidget())}};
}
void MainWindow::clearWorkspace() {
  restoringWorkspace = true;
  workspaces->clearViews();
  closedViews.clear();
  restoringWorkspace = false;
}
void MainWindow::restoreWorkspace(const QJsonObject &state) {
  restoringWorkspace = true;
  for (const auto &entry : state["views"].toArray()) {
    const auto saved = entry.toObject();
    if (saved["type"].toString() == "tf") {
      const auto settings = saved["state"].toObject();
      autoMaskWorkflow->settingsView()->restoreViewState(settings["autoMask"].toObject());
      workspaces->setTabText(0, saved["title"].toString(tr("TF view")));
      findChild<QAction *>("actionStation_navigator")->setChecked(settings["navigator"].toBool(true));
      mapCoordinateSwitch->setCurrentIndex(settings["coordinates"].toInt());
      const auto stations = ui->stationList->findItems(settings["station"].toString(), Qt::MatchExactly);
      if (!stations.isEmpty())
        ui->stationList->setCurrentItem(stations.front());
      for (int i = 0; i < ui->responsesList->count(); ++i)
        if (ui->responsesList->item(i)->toolTip() == settings["response"].toString())
          ui->responsesList->setCurrentRow(i);
      ui->responsesToggle->setChecked(settings["responses"].toBool());
      ui->splitter_5->restoreState(QByteArray::fromBase64(settings["splitter"].toString().toLatin1()));
      WorkspaceView::restorePlotRanges(tfView, settings["plots"].toObject());
    } else
      openWorkspace(saved["type"].toString(), true, saved);
  }
  workspaces->setCurrentIndex(0);
  for (int i = 0; i < workspaces->count(); ++i)
    if (workspaces->id(workspaces->widget(i)) == state["active"].toString())
      workspaces->setCurrentIndex(i);
  restoringWorkspace = false;
  activateWorkspace();
}

void MainWindow::exportWorkspace(const QString &format) {
  auto *view = workspaces->currentWidget();
  if (format != "png") {
    if (auto *analysis = dynamic_cast<WorkspaceView *>(view))
      if (auto *button = analysis->exportButton(format))
        button->click();
    return;
  }
  auto title = workspaces->tabText(workspaces->currentIndex());
  title.replace(QRegularExpression("[^\\w .-]"), "_");
  const auto path = QFileDialog::getSaveFileName(
      this, tr("Export active view"), lastDirectory + "/" + title + ".png", tr("PNG image (*.png)"));
  if (path.isEmpty())
    return;
  rememberDirectory(path);
  auto image = (workspaces->type(view) == "tf" ? tfView : view)->grab();
  if (workspaces->type(view) == "map")
    image = sharedMapHost->currentWidget()->grab();
  else if (workspaces->type(view) == "profile" || workspaces->type(view) == "comparison" ||
           workspaces->type(view) == "fit") {
    const auto map = sharedMapHost->currentWidget()->grab();
    QPixmap combined(map.width() + image.width(), std::max(map.height(), image.height()));
    combined.setDevicePixelRatio(image.devicePixelRatio());
    combined.fill(palette().color(QPalette::Window));
    {
      QPainter painter(&combined);
      painter.drawPixmap(0, 0, map);
      painter.drawPixmap(QPointF(map.width() / map.devicePixelRatio(), 0), image);
    }
    image = combined;
  }
  if (!image.save(path, "PNG"))
    QMessageBox::warning(this, tr("Export"), tr("Could not save the image."));
}
