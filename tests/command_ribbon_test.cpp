#include "ui/workspace/CommandRibbon.h"
#include "ui/workspace/CommandIcons.h"
#include "app/MainWindow.h"
#include "ui/workspace/WorkspaceView.h"
#include "ui/maps/SharedMapInteraction.h"
#include <QListWidget>
#include <QTreeWidget>
#include <QApplication>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMenuBar>
#include <QPalette>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QStackedWidget>
#include <QTabBar>
#include <QTemporaryDir>
#include <QToolButton>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition, const char *message)
{
  if(!condition) throw std::runtime_error(message);
}
template<class T> T *control(QObject &parent, const char *name)
{
  auto *result = parent.findChild<T *>(name);
  check(result, name);
  return result;
}
QToolButton *button(QObject &parent, QAction *action)
{
  for(auto *candidate : parent.findChildren<QToolButton *>())
    if(candidate->defaultAction() == action) return candidate;
  throw std::runtime_error("Missing action button");
}
void key(QWidget &window, int code, Qt::KeyboardModifiers modifiers)
{
  QKeyEvent press(QEvent::KeyPress, code, modifiers), release(QEvent::KeyRelease, code, modifiers);
  QApplication::sendEvent(&window, &press);
  QApplication::sendEvent(&window, &release);
  QApplication::processEvents();
}
void iconRendering()
{
  const auto original = QApplication::palette();
  const auto icon = commandIcon("folder-open");
  check(!icon.isNull(), "Bundled SVG resource available");
  auto palette = original;
  const auto verify = [&](QIcon::Mode mode, QColor expected, qreal scale) {
    const auto pixmap = icon.pixmap(QSize(20, 20), scale, mode);
    check(pixmap.size() == QSize(20 * scale, 20 * scale), "High DPI icon resolution");
    const auto image = pixmap.toImage();
    bool ink = false, transparent = false;
    for(int y = 0; y < image.height(); ++y)
      for(int x = 0; x < image.width(); ++x) {
        const auto pixel = image.pixelColor(x, y);
        transparent |= pixel.alpha() == 0;
        if(pixel.alpha() > 250) {
          ink = true;
          check(pixel.rgb() == expected.rgb(), "Icon uses current palette colour");
        }
      }
    check(ink && transparent, "Icon contains visible glyph and transparent background");
  };
  palette.setColor(QPalette::Active, QPalette::ButtonText, QColor(25, 30, 40));
  QApplication::setPalette(palette);
  verify(QIcon::Normal, QColor(25, 30, 40), 1);
  palette.setColor(QPalette::Active, QPalette::ButtonText, QColor(235, 240, 250));
  palette.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(100, 110, 120));
  palette.setColor(QPalette::Active, QPalette::HighlightedText, Qt::white);
  QApplication::setPalette(palette);
  verify(QIcon::Normal, QColor(235, 240, 250), 2);
  verify(QIcon::Disabled, QColor(100, 110, 120), 1);
  verify(QIcon::Selected, Qt::white, 1);
  QApplication::setPalette(original);
}
void presentation()
{
  QMainWindow window;
  window.setCentralWidget(new QWidget);
  auto *menu = window.menuBar()->addMenu("Commands");
  auto *action = menu->addAction("Shared action");
  action->setIconText("Shared");
  action->setCheckable(true);
  action->setShortcut(QKeySequence("Ctrl+K"));
  auto *about = menu->addAction("About");
  int triggered = 0;
  QObject::connect(action, &QAction::triggered, [&] { ++triggered; });
  std::vector<QAction *> many;
  for(int i = 0; i < 16; ++i) many.push_back(menu->addAction(QString("Command %1").arg(i)));
  auto *ribbon = new CommandRibbon(window, {}, {{"First", {{"Shared", {action}}}}, {"Second", {{"Many commands", many}}}}, about);
  auto *tabs = control<QTabBar>(*ribbon, "ribbonTabs");
  auto *body = control<QStackedWidget>(*ribbon, "ribbonPages");
  auto *classic = control<QAction>(*ribbon, "actionClassic_menus");
  auto *collapsed = control<QAction>(*ribbon, "actionCollapse_ribbon");
  window.resize(650, 400);
  window.show(); window.activateWindow();
  QApplication::processEvents();
  check(!window.menuBar()->isVisible() && body->isVisible(), "Default ribbon layout");
  auto *grip = control<QWidget>(*ribbon, "ribbonResizeHandle");
  const int initialHeight = body->height();
  const QPointF local = grip->rect().center(), global = grip->mapToGlobal(local.toPoint());
  QMouseEvent press(QEvent::MouseButtonPress, local, global, Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
  QApplication::sendEvent(grip, &press);
  QMouseEvent move(QEvent::MouseMove, local + QPointF(0, 40), global + QPointF(0, 40), Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
  QApplication::sendEvent(grip, &move);
  QMouseEvent release(QEvent::MouseButtonRelease, local, global + QPointF(0, 40), Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(grip, &release);
  QApplication::processEvents();
  check(body->height() == initialHeight + 40, "Dragging ribbon edge must resize it");
  key(*grip, Qt::Key_Up, Qt::NoModifier);
  check(body->height() == initialHeight + 32, "Keyboard ribbon resize failed");
  key(*grip, Qt::Key_Home, Qt::NoModifier);
  check(body->height() == initialHeight, "Default ribbon height was not restored");
  key(*grip, Qt::Key_Down, Qt::NoModifier);
  check(QSettings().value("interface/ribbonHeight").toInt() == body->height(), "Ribbon height not saved");

  auto *shared = button(*ribbon, action);
  shared->click();
  check(triggered == 1 && action->isChecked() && shared->isChecked(), "Single shared command dispatch");
  action->setChecked(false);
  action->setEnabled(false);
  check(!shared->isChecked() && !shared->isEnabled(), "Action state propagates to ribbon");
  action->setEnabled(true);
  tabs->setCurrentIndex(1);
  key(window, Qt::Key_K, Qt::ControlModifier);
  check(triggered == 2 && action->isChecked(), "Shortcut works on inactive tab with hidden menus");
  collapsed->setChecked(true);
  QApplication::processEvents();
  check(!body->isVisible() && tabs->isVisible() && !grip->isVisible(), "Ribbon collapse");
  check(QSettings().value("interface/ribbonCollapsed").toBool(), "Collapse preference saved");
  tabs->tabBarClicked(0);
  check(!collapsed->isChecked() && body->isVisible(), "Click tab expands ribbon");
  key(window, Qt::Key_F1, Qt::ControlModifier);
  check(collapsed->isChecked(), "Collapse shortcut");
  classic->setChecked(true);
  check(window.menuBar()->isVisible() && !tabs->isVisible() && !body->isVisible(), "Classic menu layout");
  check(control<QToolButton>(*ribbon, "ribbonInterface")->isVisible(), "Layout switch stays available");
  key(window, Qt::Key_K, Qt::ControlModifier);
  check(triggered == 3, "No duplicate shortcuts in classic layout");
  classic->setChecked(false);
  collapsed->setChecked(false);
  window.resize(380, 400);
  QApplication::processEvents();
  auto *scroll = qobject_cast<QScrollArea *>(body->currentWidget());
  check(scroll && scroll->horizontalScrollBar()->maximum() > 0, "Commands scroll in narrow windows");
  auto *interfaceButton = control<QToolButton>(*ribbon, "ribbonInterface");
  check(interfaceButton->mapTo(&window, interfaceButton->rect().topRight()).x() < window.width(), "Interface remains reachable at narrow width");
  classic->setChecked(true);
  collapsed->setChecked(true);
}
void workspaceLifecycle()
{
  struct View : WorkspaceView {
    QCustomPlot *plot = new QCustomPlot(this);
    QCustomPlot *map = new QCustomPlot(this);
    int refreshes = 0;
    bool deferred = false;
    View() {
      plot->setObjectName("section"); map->setObjectName("map");
      setCentralMaps({map}, {}, {});
    }
    void refreshData(const std::shared_ptr<MTSurveyData> &, const Responses &) override {
      ++refreshes;
      plot->xAxis->setRange(1., 10.); map->xAxis->setRange(1., 10.);
    }
    bool hasDeferredPlots() const override { return deferred; }
  };
  QWidget host;
  auto view = std::make_unique<View>();
  QPointer<QCustomPlot> map = view->map;
  map->setParent(&host);
  auto survey = std::make_shared<MTSurveyData>("Workspace test");
  view->refreshSurvey(survey, {});
  view->plot->xAxis->setRange(2., 3.); map->xAxis->setRange(4., 6.);
  view->refreshSurvey(survey, {});
  check(view->refreshes == 1, "Clean views must not recompute data");
  view->markSurveyDirty(); view->refreshSurvey(survey, {});
  check(view->refreshes == 2 && view->plot->xAxis->range() == QCPRange(2., 3.) &&
        map->xAxis->range() == QCPRange(4., 6.), "Refresh lost zoom in the view or reparented map");
  view->deferred = true;
  view->markSurveyDirty(); view->refreshSurvey(survey, {});
  view->plot->xAxis->setRange(10., 20.); view->restorePendingPlotRanges();
  check(view->plot->xAxis->range() == QCPRange(2., 3.), "Deferred refresh lost zoom");
  view->markSurveyDirty(); view->refreshSurvey(survey, {});
  SurveyCoordinates frame; frame.utm = true; frame.zone = 32;
  survey->set_coordinates(frame);
  view->markSurveyDirty(); view->refreshSurvey(survey, {});
  view->plot->xAxis->setRange(10., 20.); view->restorePendingPlotRanges();
  check(map->xAxis->range() == QCPRange(1., 10.) && view->plot->xAxis->range() == QCPRange(10., 20.),
        "Changing coordinates reapplied stale plot ranges");

  auto *list = new QListWidget(view.get()); list->setObjectName("stations");
  auto *site = new QListWidgetItem("Displayed name", list);
  site->setData(Qt::UserRole, "station-id"); site->setCheckState(Qt::Checked);
  auto *tree = new QTreeWidget(view.get()); tree->setObjectName("responses");
  auto *response = new QTreeWidgetItem(tree, {"Response"});
  response->setData(0, Qt::UserRole, "response-id"); response->setCheckState(0, Qt::Checked);
  const auto saved = view->saveViewState();
  site->setCheckState(Qt::Unchecked); site->setText("Renamed display label");
  response->setCheckState(0, Qt::Unchecked);
  int changes = 0;
  QObject::connect(list, &QListWidget::itemChanged, [&] { ++changes; });
  QObject::connect(tree, &QTreeWidget::itemChanged, [&] { ++changes; });
  view->restoreViewState(saved);
  check(site->checkState() == Qt::Checked && response->checkState(0) == Qt::Checked && changes == 0,
        "Restoring view selections lost stable identities or emitted editing signals");
  view.reset();
  check(map.isNull(), "Closing a view did not release its reparented map");
}

void sharedMapInteractions()
{
  QCustomPlot map;
  const auto original = QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables;
  map.setInteractions(original); map.setCursor(Qt::PointingHandCursor);
  {
    SharedMapInteraction session(&map);
    check(session.activate(), "Could not activate shared map");
    const QCursor drawing(Qt::CrossCursor);
    session.configure(QCP::iRangeDrag | QCP::iSelectPlottables, &drawing);
    check(map.interactions() == QCP::iRangeZoom && map.cursor().shape() == Qt::CrossCursor,
          "Drawing mode did not suspend map dragging and selection");
    session.activate(); session.configure();
    check(map.interactions() == original && map.cursor().shape() == Qt::PointingHandCursor,
          "Repeated activation overwrote the original map settings");
    session.configure(QCP::iSelectPlottables);
  }
  check(map.interactions() == original && map.cursor().shape() == Qt::PointingHandCursor,
        "Closing a map session did not restore interaction settings");
  auto *temporary = new QCustomPlot;
  SharedMapInteraction session(temporary); session.activate();
  delete temporary;
  check(!session.active() && !session.activate(), "Map session retained a destroyed map");
}

void integration()
{
  MainWindow window;
  auto *ribbon = control<QToolBar>(window, "commandRibbon");
  auto *classic = control<QAction>(*ribbon, "actionClassic_menus");
  auto *collapsed = control<QAction>(*ribbon, "actionCollapse_ribbon");
  check(classic->isChecked() && collapsed->isChecked(), "Preferences restored in new window");
  check(control<QStackedWidget>(*ribbon, "ribbonPages")->height() == QSettings().value("interface/ribbonHeight").toInt(), "Ribbon height not restored in new window");
  classic->setChecked(false);
  collapsed->setChecked(false);
  auto *tabs = control<QTabBar>(*ribbon, "ribbonTabs");
  check(tabs->count() == 8, "Permanent and TF context tabs");
  classic->setChecked(true); collapsed->setChecked(true);
  control<QAction>(window, "actionPlot_axis_ranges")->trigger();
  check(!classic->isChecked() && !collapsed->isChecked() && tabs->tabText(tabs->currentIndex()) == "TF appearance",
        "Classic view command must open the View ribbon");
  check(control<QComboBox>(window, "viewRangeAxis")->isVisible(), "Inline axis controls not visible");
  for(const auto *name: {"actionLoad_EDI", "actionNew_view", "actionRename_view", "actionProfile_analysis", "actionView_png"}) {
    auto *action = control<QAction>(window, name); button(*ribbon, action);
    check(window.actions().contains(action), "Window command registration");
  }
  check(control<QAction>(window, "actionLoad_EDI")->isEnabled(), "Import available without survey");
  check(!control<QAction>(window, "actionLoad_EDI")->icon().isNull(), "Import has a command icon");
  check(!control<QAction>(window, "actionSave_project")->isEnabled(), "Save needs survey");
  check(!control<QAction>(window, "actionImpedance_validity")->isEnabled(), "Analysis needs survey");
  check(!control<QAction>(window, "actionUndo")->isEnabled(), "Empty history disables undo");
  check(control<QAction>(window, "actionSave_as_PDF")->text() == "Save plot image…", "Image export correctly labelled");
  QApplication::processEvents();
  const auto artifacts = qEnvironmentVariable("EDITOOLS_TEST_ARTIFACT_DIR");
  if(!artifacts.isEmpty()) {
    window.showNormal(); window.resize(1280, 850);
    tabs->setCurrentIndex(0);
    QApplication::processEvents();
    window.grab().save(artifacts + "/ribbon-survey.png");
    tabs->setCurrentIndex(1);
    QApplication::processEvents();
    window.grab().save(artifacts + "/ribbon-edit.png");
  }
}
}
int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("EDIToolsTests");
  QCoreApplication::setApplicationName("CommandRibbon");
  QTemporaryDir settings;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, settings.path());
  try { iconRendering(); presentation(); workspaceLifecycle(); sharedMapInteractions(); integration(); }
  catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
  return 0;
}
