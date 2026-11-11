#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include/WorkspaceTabs.h"
#include "include/WorkspaceView.h"
#include "include/MapBackground.h"
#include "include/SurveyHistory.h"
#include "ProfileWindow.h"
#include "PeriodMapWindow.h"
#include "FitStatisticsWindow.h"
#include "StationComparisonWindow.h"
#include "PolarityWindow.h"
#include "PeriodLayoutWindow.h"
#include "ImpedanceAnalysisWindows.h"
#include "AutoMaskDialog.h"
#include "AutoMaskWorkflow.h"
#include <QAbstractButton>
#include <QFileDialog>
#include <QJsonArray>
#include <QMessageBox>
#include <QRegularExpression>

namespace {
const QStringList kinds{"profile", "map", "fit", "comparison", "validity", "average", "polarity", "periods", "automask"};
}

void MainWindow::setupWorkspace()
{
  workspaces = new WorkspaceTabs;
  ui->splitter_5->replaceWidget(ui->splitter_5->indexOf(tfView), workspaces);
  ui->splitter_5->setStretchFactor(0, 0);
  ui->splitter_5->setStretchFactor(1, 1);
  ui->splitter_5->setCollapsible(0, false);
  workspaces->addView(tfView, "tf", tr("TF view"), "tf");
  connect(workspaces, &QTabWidget::currentChanged, this, [this] { activateWorkspace(); });
  connect(workspaces, &WorkspaceTabs::closing, this, [this](QWidget *view) {
    if(!restoringWorkspace) { closedViews.append(workspaceSnapshot(view)); if(closedViews.size() > 20) closedViews.removeFirst(); }
    if(ribbon) ribbon->releaseContext(workspaces->id(view));
  });
  connect(workspaces, &WorkspaceTabs::titleChanged, this, [](QWidget *view, const QString &title) { view->setProperty("viewTitle", title); });
  connect(workspaces, &WorkspaceTabs::duplicateRequested, this, [this](QWidget *view) {
    auto saved = workspaceSnapshot(view); saved.remove("id"); saved["title"] = saved["title"].toString() + tr(" copy");
    openWorkspace(workspaces->type(view), true, saved);
  });
  connect(workspaces, &WorkspaceTabs::reopenRequested, this, [this] {
    if(closedViews.isEmpty()) return;
    const auto saved = closedViews.takeLast(); openWorkspace(saved["type"].toString(), true, saved);
  });
  auto *menu = menuBar()->addMenu(tr("Views"));
  auto *newView = menu->addAction(tr("New view")); newView->setObjectName("actionNew_view");
  auto *choices = new QMenu(menu); newView->setMenu(choices);
  const QStringList titles{tr("2D profile"), tr("Period map"), tr("Fit statistics"), tr("Station comparison"), tr("Rho+ validity"), tr("Survey average"), tr("Polarity review"), tr("Period layout"), tr("Auto-mask review")};
  for(int i = 0; i < kinds.size(); ++i) connect(choices->addAction(titles[i]), &QAction::triggered, this, [this, i] { openWorkspace(kinds[i], true); });
  auto *duplicate = menu->addAction(tr("Duplicate view")); duplicate->setObjectName("actionDuplicate_view");
  connect(duplicate, &QAction::triggered, this, [this] { emit workspaces->duplicateRequested(workspaces->currentWidget()); });
  auto *rename = menu->addAction(tr("Rename tab")); rename->setObjectName("actionRename_view");
  connect(rename, &QAction::triggered, workspaces, &WorkspaceTabs::renameCurrent);
  auto *close = menu->addAction(tr("Close tab")); close->setObjectName("actionClose_view");
  connect(close, &QAction::triggered, workspaces, &WorkspaceTabs::closeCurrent);
  auto *reopen = menu->addAction(tr("Reopen closed tab")); reopen->setObjectName("actionReopen_view");
  connect(reopen, &QAction::triggered, this, [this] { emit workspaces->reopenRequested(); });
  auto *navigator = menu->addAction(tr("Station navigator")); navigator->setObjectName("actionStation_navigator"); navigator->setCheckable(true); navigator->setChecked(true);
  connect(navigator, &QAction::toggled, ui->splitter_3, &QWidget::setVisible);
  for(const auto &format: {QString("png"), QString("pdf"), QString("csv")}) {
    auto *action = menu->addAction(tr("Active view %1…").arg(format.toUpper())); action->setObjectName("actionView_" + format);
    connect(action, &QAction::triggered, this, [this, format] { exportWorkspace(format); });
  }
}

QWidget *MainWindow::openWorkspace(const QString &type, bool create, const QJsonObject &saved)
{
  if(!workspaces || !mtSurvey || !kinds.contains(type)) return nullptr;
  if(!create) for(int i = workspaces->count() - 1; i > 0; --i) if(workspaces->type(workspaces->widget(i)) == type) {
    workspaces->setCurrentIndex(i); return workspaces->widget(i);
  }
  WorkspaceView *view = nullptr;
  const auto changed = [this] { updatePlots(false); refreshAnalysisWindows(); };
  const auto copy = [this](std::shared_ptr<MTSurveyData> data) { openSurveyCopy(std::move(data)); };
  if(type == "profile") { profile = new ProfileWindow(this, copy); profile->setMapLayers(mapBackground->layers()); view = profile; }
  else if(type == "map") {
    periodMaps = new PeriodMapWindow(this, [this, changed] { history->record(tr("Map masks")); changed(); });
    auto settings = periodMaps->displaySettings(); settings.mapLayers = mapBackground->layers(); settings.names = mapStationNames->isChecked();
    periodMaps->setDisplaySettings(settings); periodMaps->setLinkTensorMasks(linkTensorMasksAction->isChecked()); view = periodMaps;
  } else if(type == "fit") { fitStatistics = new FitStatisticsWindow(this, [this] { refreshAnalysisWindows(); }); view = fitStatistics; }
  else if(type == "comparison") { stationComparison = new StationComparisonWindow(this); view = stationComparison; }
  else if(type == "validity") {
    impedanceValidity = new ImpedanceValidityWindow(this, [this](const ImpedanceAnalysis::RhoPlusOptions &options) { rhoPlusOptions = options; updatePlots(false); }, changed);
    impedanceValidity->setLinkTensorMasks(linkTensorMasksAction->isChecked()); impedanceValidity->setFitOptions(rhoPlusOptions); impedanceValidity->setHistory(history.get()); view = impedanceValidity;
  } else if(type == "average") { surveyAverage = new SurveyAverageWindow(this); view = surveyAverage; }
  else if(type == "polarity") { polarity = new PolarityWindow(this, *history, changed); view = polarity; }
  else if(type == "periods") { periodLayout = new PeriodLayoutWindow(this, copy); view = periodLayout; }
  else if(type == "automask") {
    std::set<std::string> selection;
    for(const auto *item: ui->stationList->selectedItems()) selection.insert(item->text().toStdString());
    if(selection.empty() && ui->stationList->currentItem()) selection.insert(ui->stationList->currentItem()->text().toStdString());
    view = new AutoMaskReview(this, mtSurvey, selection, linkTensorMasksAction->isChecked(), autoMaskWorkflow->currentOptions(), autoMaskWorkflow->pairDiagonals(), [this, changed] { history->record(tr("Auto masking")); changed(); });
  }
  if(!view) return nullptr;
  view->setEmbedded(); view->setProperty("workspaceDirty", true);
  refreshWorkspace(view);
  const auto station = ui->stationList->currentItem() ? ui->stationList->currentItem()->text() : QString();
  if(type == "comparison") stationComparison->setData(mtSurvey, {station});
  if(type == "polarity") polarity->setData(mtSurvey, station);
  if(type == "validity") impedanceValidity->setData(mtSurvey, station);
  QString title = saved["title"].toString();
  if(title.isEmpty()) {
    title = view->windowTitle(); int number = 1;
    for(int i = 1; i < workspaces->count(); ++i) if(workspaces->type(workspaces->widget(i)) == type) ++number;
    if(number > 1) title += QString(" %1").arg(number);
  }
  view->setProperty("viewTitle", title);
  if(saved.contains("state")) { view->restoreViewState(saved["state"].toObject()); view->setProperty("workspaceDirty", true); }
  connect(view, &WorkspaceView::toolsChanged, this, [this, view] {
    if(workspaces->currentWidget() != view) return;
    for(const auto &format: {QString("pdf"), QString("csv")}) if(auto *action = findChild<QAction *>("actionView_" + format)) {
      auto *button = view->exportButton(format); action->setEnabled(button && button->isEnabled());
    }
  });
  connect(view, &WorkspaceView::stationActivated, this, [this](const QString &name) {
    const auto items = ui->stationList->findItems(name, Qt::MatchExactly);
    if(!items.isEmpty()) ui->stationList->setCurrentItem(items.front(), QItemSelectionModel::NoUpdate);
  });
  workspaces->addView(view, type, title, saved["id"].toString());
  if(ribbon && !restoringWorkspace) ribbon->setContext(workspaces->id(view), view->toolPages(), true);
  return view;
}

void MainWindow::refreshWorkspace(QWidget *view)
{
  if(!view || !view->property("workspaceDirty").toBool()) return;
  const auto ranges = WorkspaceView::savePlotRanges(view);
  const bool initialized = view->property("workspaceInitialized").toBool();
  const auto coordinates = mtSurvey ? QString::fromStdString(mtSurvey->coordinates().description()) : QString();
  const bool sameCoordinates = view->property("workspaceCoordinates").toString() == coordinates;
  view->setProperty("workspaceCoordinates", coordinates);
  view->setProperty("workspaceDirty", false);
  if(auto *v = dynamic_cast<ProfileWindow *>(view)) v->setData(mtSurvey, mtResponses);
  else if(auto *v = dynamic_cast<PeriodMapWindow *>(view)) v->setData(mtSurvey, mtResponses);
  else if(auto *v = dynamic_cast<FitStatisticsWindow *>(view)) v->setData(mtSurvey, mtResponses);
  else if(auto *v = dynamic_cast<StationComparisonWindow *>(view)) v->setData(mtSurvey);
  else if(auto *v = dynamic_cast<PolarityWindow *>(view)) v->setData(mtSurvey);
  else if(auto *v = dynamic_cast<ImpedanceValidityWindow *>(view)) v->setData(mtSurvey);
  else if(auto *v = dynamic_cast<SurveyAverageWindow *>(view)) v->setData(mtSurvey);
  else if(auto *v = dynamic_cast<PeriodLayoutWindow *>(view)) v->setData(mtSurvey);
  else if(auto *v = dynamic_cast<AutoMaskReview *>(view)) v->setData(mtSurvey);
  view->setProperty("workspaceInitialized", true);
  if(initialized && sameCoordinates) {
    WorkspaceView::restorePlotRanges(view, ranges);
    if(dynamic_cast<ImpedanceValidityWindow *>(view) || dynamic_cast<SurveyAverageWindow *>(view)) view->setProperty("pendingPlotRanges", ranges);
  }
}

void MainWindow::activateWorkspace()
{
  if(!ribbon || !workspaces->currentWidget()) return;
  auto *view = workspaces->currentWidget(); refreshWorkspace(view);
  auto *analysis = dynamic_cast<WorkspaceView *>(view);
  ribbon->setContext(workspaces->id(view), analysis ? analysis->toolPages() : tfTools);
  for(const auto &name: {"actionClose_view", "actionDuplicate_view"}) findChild<QAction *>(name)->setEnabled(analysis);
  findChild<QAction *>("actionReopen_view")->setEnabled(!closedViews.isEmpty());
  for(const auto &format: {QString("pdf"), QString("csv")}) {
    auto *action = findChild<QAction *>("actionView_" + format);
    auto *button = analysis ? analysis->exportButton(format) : nullptr;
    action->setEnabled(button && button->isEnabled());
  }
}

QJsonObject MainWindow::workspaceSnapshot(QWidget *view) const
{
  QJsonObject state;
  if(auto *analysis = dynamic_cast<WorkspaceView *>(view)) state = analysis->saveViewState();
  else state = {{"station", ui->stationList->currentItem() ? ui->stationList->currentItem()->text() : QString()},
    {"response", ui->responsesList->currentItem() ? ui->responsesList->currentItem()->toolTip() : QString()},
    {"autoMask", autoMaskWorkflow->settingsView()->saveViewState()},
    {"coordinates", mapCoordinateSwitch->currentIndex()}, {"responses", ui->responsesToggle->isChecked()},
    {"navigator", findChild<QAction *>("actionStation_navigator")->isChecked()},
    {"splitter", QString::fromLatin1(ui->splitter_5->saveState().toBase64())}, {"plots", WorkspaceView::savePlotRanges(tfView)}};
  return {{"id", workspaces->id(view)}, {"type", workspaces->type(view)}, {"title", workspaces->tabText(workspaces->indexOf(view))}, {"state", state}};
}
QJsonObject MainWindow::saveWorkspace() const
{
  QJsonArray views;
  for(int i = 0; i < workspaces->count(); ++i) views.append(workspaceSnapshot(workspaces->widget(i)));
  return {{"schema", 1}, {"views", views}, {"active", workspaces->id(workspaces->currentWidget())}};
}
void MainWindow::clearWorkspace()
{
  restoringWorkspace = true; workspaces->clearViews(); closedViews.clear(); restoringWorkspace = false;
}
void MainWindow::restoreWorkspace(const QJsonObject &state)
{
  restoringWorkspace = true;
  for(const auto &entry: state["views"].toArray()) {
    const auto saved = entry.toObject();
    if(saved["type"].toString() == "tf") {
      const auto settings = saved["state"].toObject();
      autoMaskWorkflow->settingsView()->restoreViewState(settings["autoMask"].toObject());
      workspaces->setTabText(0, saved["title"].toString(tr("TF view")));
      findChild<QAction *>("actionStation_navigator")->setChecked(settings["navigator"].toBool(true));
      mapCoordinateSwitch->setCurrentIndex(settings["coordinates"].toInt());
      const auto stations = ui->stationList->findItems(settings["station"].toString(), Qt::MatchExactly);
      if(!stations.isEmpty()) ui->stationList->setCurrentItem(stations.front());
      for(int i = 0; i < ui->responsesList->count(); ++i) if(ui->responsesList->item(i)->toolTip() == settings["response"].toString()) ui->responsesList->setCurrentRow(i);
      ui->responsesToggle->setChecked(settings["responses"].toBool());
      ui->splitter_5->restoreState(QByteArray::fromBase64(settings["splitter"].toString().toLatin1()));
      WorkspaceView::restorePlotRanges(tfView, settings["plots"].toObject());
    } else openWorkspace(saved["type"].toString(), true, saved);
  }
  workspaces->setCurrentIndex(0);
  for(int i = 0; i < workspaces->count(); ++i) if(workspaces->id(workspaces->widget(i)) == state["active"].toString()) workspaces->setCurrentIndex(i);
  restoringWorkspace = false; activateWorkspace();
}

void MainWindow::exportWorkspace(const QString &format)
{
  auto *view = workspaces->currentWidget();
  if(format != "png") {
    if(auto *analysis = dynamic_cast<WorkspaceView *>(view)) if(auto *button = analysis->exportButton(format)) button->click();
    return;
  }
  auto title = workspaces->tabText(workspaces->currentIndex()); title.replace(QRegularExpression("[^\\w .-]"), "_");
  const auto path = QFileDialog::getSaveFileName(this, tr("Export active view"), lastDirectory + "/" + title + ".png", tr("PNG image (*.png)"));
  if(path.isEmpty()) return;
  rememberDirectory(path);
  if(!(workspaces->type(view) == "tf" ? tfView : view)->grab().save(path, "PNG")) QMessageBox::warning(this, tr("Export"), tr("Could not save the image."));
}
