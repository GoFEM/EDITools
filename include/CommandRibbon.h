#pragma once

#include <QToolBar>
#include <vector>

class QMainWindow;
class QTabBar;
class QStackedWidget;
class QToolButton;

// Presentation only: menus, shortcuts and ribbon buttons share the same actions.
class CommandRibbon : public QToolBar
{
public:
  struct Group { QString title; std::vector<QAction *> actions; };
  struct Page { QString title; std::vector<Group> groups; };

  CommandRibbon(QMainWindow &window, const std::vector<QAction *> &quickActions,
                const std::vector<Page> &pages, QAction *about);

private:
  void updateDisplay();
  QMainWindow &window;
  QTabBar *tabs;
  QStackedWidget *body;
  QToolButton *collapseButton;
  QAction *classic;
  QAction *collapsed;
};
