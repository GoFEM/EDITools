#include "include/CommandRibbon.h"
#include "mainwindow.h"
#include <QApplication>
#include <QKeyEvent>
#include <QMenuBar>
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
  check(!body->isVisible() && tabs->isVisible(), "Ribbon collapse");
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
  auto *interface = control<QToolButton>(*ribbon, "ribbonInterface");
  check(interface->mapTo(&window, interface->rect().topRight()).x() < window.width(), "Interface remains reachable at narrow width");
  classic->setChecked(true);
  collapsed->setChecked(true);
}
void integration()
{
  MainWindow window;
  auto *ribbon = control<QToolBar>(window, "commandRibbon");
  auto *classic = control<QAction>(*ribbon, "actionClassic_menus");
  auto *collapsed = control<QAction>(*ribbon, "actionCollapse_ribbon");
  check(classic->isChecked() && collapsed->isChecked(), "Preferences restored in new window");
  classic->setChecked(false);
  collapsed->setChecked(false);
  auto *tabs = control<QTabBar>(*ribbon, "ribbonTabs");
  check(tabs->count() == 5, "Five workflow tabs");
  for(auto *top : window.menuBar()->actions())
    for(auto *action : top->menu()->actions()) {
      if(action->isSeparator() || action->objectName() == "actionAbout_EDI_Tools") continue;
      button(*ribbon, action);
      check(window.actions().contains(action), "Command shortcut registered on main window");
    }
  check(control<QAction>(window, "actionLoad_EDI")->isEnabled(), "Import available without survey");
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
  try { presentation(); integration(); }
  catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
  return 0;
}
