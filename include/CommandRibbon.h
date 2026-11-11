#pragma once

#include <QToolBar>
#include <QPointer>
#include <QHash>
#include <vector>

class QMainWindow;
class QTabBar;
class QStackedWidget;
class QToolButton;

// Presentation only: menus, shortcuts and ribbon buttons share the same actions.
class CommandRibbon : public QToolBar
{
public:
  struct Group { QString title; std::vector<QAction *> actions; QWidget *controls = nullptr; };
  struct Page { QString title; std::vector<Group> groups; QString id; };

  CommandRibbon(QMainWindow &window, const std::vector<QAction *> &quickActions,
                const std::vector<Page> &pages, QAction *about);
  void showPage(int index);
  void showPage(const QString &id);
  void setContext(const QString &key, const std::vector<Page> &pages, bool select = false);
  void releaseContext(const QString &key);

private:
  bool eventFilter(QObject *object, QEvent *event) override;
  void resizeBody(int height, bool save);
  QWidget *resizeHandle = nullptr;
  int defaultBodyHeight = 112;
  int dragStartY = 0, dragStartHeight = 0;
  bool dragging = false;
  void updateDisplay();
  void appendPage(const Page &page, const QString &context = {});
  QToolButton *button(QAction *action, QWidget *parent, bool quick = false);
  struct BorrowedControl { QPointer<QWidget> widget, parent; };
  struct PageRecord { QString id, context; QPointer<QWidget> widget; std::vector<BorrowedControl> controls; };
  std::vector<PageRecord> records;
  QString activeContext;
  QMainWindow &window;
  QTabBar *tabs;
  QStackedWidget *body;
  QToolButton *collapseButton;
  QAction *classic;
  QAction *collapsed;
};
