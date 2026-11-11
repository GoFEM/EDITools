#include "include/WorkspaceTabs.h"
#include <QAction>
#include <QLineEdit>
#include <QMenu>
#include <QTabBar>
#include <QUuid>
#include <QSignalBlocker>

WorkspaceTabs::WorkspaceTabs(QWidget *parent): QTabWidget(parent)
{
  setObjectName("workspaceTabs"); setMovable(true); setTabsClosable(true); setDocumentMode(true);
  tabBar()->setObjectName("workspaceTabBar"); tabBar()->setExpanding(false); tabBar()->setUsesScrollButtons(true);
  connect(this, &QTabWidget::tabCloseRequested, this, &WorkspaceTabs::closeTab);
  connect(tabBar(), &QTabBar::tabBarDoubleClicked, this, &WorkspaceTabs::renameTab);
  connect(tabBar(), &QTabBar::tabMoved, this, [this] {
    for(int i = 1; i < count(); ++i) if(type(widget(i)) == "tf") { tabBar()->moveTab(i, 0); break; }
  });
  tabBar()->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(tabBar(), &QWidget::customContextMenuRequested, this, [this](const QPoint &position) {
    const int index = tabBar()->tabAt(position); if(index < 0) return;
    setCurrentIndex(index);
    QMenu menu(this);
    menu.addAction(tr("Rename tab"), this, &WorkspaceTabs::renameCurrent);
    auto *duplicate = menu.addAction(tr("Duplicate view"), this, [this] { emit duplicateRequested(currentWidget()); });
    duplicate->setEnabled(type(currentWidget()) != "tf");
    auto *close = menu.addAction(tr("Close tab"), this, &WorkspaceTabs::closeCurrent);
    close->setEnabled(type(currentWidget()) != "tf");
    menu.addAction(tr("Reopen closed tab"), this, [this] { emit reopenRequested(); });
    menu.exec(tabBar()->mapToGlobal(position));
  });
  auto *close = new QAction(this); close->setShortcut(QKeySequence::Close); addAction(close);
  connect(close, &QAction::triggered, this, &WorkspaceTabs::closeCurrent);
  auto *reopen = new QAction(this); reopen->setShortcut(QKeySequence("Ctrl+Shift+T")); addAction(reopen);
  connect(reopen, &QAction::triggered, this, [this] { emit reopenRequested(); });
  auto *rename = new QAction(tabBar()); rename->setShortcut(Qt::Key_F2); rename->setShortcutContext(Qt::WidgetShortcut); tabBar()->addAction(rename);
  tabBar()->setFocusPolicy(Qt::StrongFocus); connect(rename, &QAction::triggered, this, &WorkspaceTabs::renameCurrent);
}
QString WorkspaceTabs::id(QWidget *view) const { return view ? view->property("workspaceId").toString() : QString(); }
QString WorkspaceTabs::type(QWidget *view) const { return view ? view->property("workspaceType").toString() : QString(); }
void WorkspaceTabs::addView(QWidget *view, const QString &kind, const QString &title, const QString &identifier)
{
  view->setProperty("workspaceType", kind);
  view->setProperty("workspaceId", identifier.isEmpty() ? QUuid::createUuid().toString(QUuid::WithoutBraces) : identifier);
  const int index = addTab(view, title); setTabToolTip(index, title + " · " + view->windowTitle());
  if(kind == "tf") { tabBar()->setTabButton(index, QTabBar::RightSide, nullptr); tabBar()->setTabButton(index, QTabBar::LeftSide, nullptr); }
  setCurrentIndex(index);
}
void WorkspaceTabs::renameCurrent() { renameTab(currentIndex()); }
void WorkspaceTabs::renameTab(int index)
{
  if(index < 0) return;
  auto *view = widget(index);
  auto *edit = new QLineEdit(tabText(index), tabBar()); edit->setObjectName("workspaceRename");
  edit->setGeometry(tabBar()->tabRect(index)); edit->show(); edit->selectAll(); edit->setFocus();
  connect(edit, &QLineEdit::editingFinished, this, [this, edit, view] {
    const int index = indexOf(view); const auto title = edit->text().trimmed();
    if(index >= 0 && !title.isEmpty()) { setTabText(index, title); setTabToolTip(index, title + " · " + view->windowTitle()); emit titleChanged(view, title); }
    edit->hide(); edit->deleteLater();
  });
}
void WorkspaceTabs::closeCurrent() { closeTab(currentIndex()); }
void WorkspaceTabs::closeTab(int index)
{
  auto *view = widget(index); if(!view || type(view) == "tf") return;
  emit closing(view); removeTab(index); delete view;
}
void WorkspaceTabs::clearViews()
{
  for(int i = count() - 1; i >= 0; --i) if(type(widget(i)) != "tf") closeTab(i);
}
