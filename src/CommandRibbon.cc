#include "include/CommandRibbon.h"

#include <QAction>
#include <QGridLayout>
#include <QLabel>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QStackedWidget>
#include <QTabBar>
#include <QToolButton>

CommandRibbon::CommandRibbon(QMainWindow &window, const std::vector<QAction *> &quickActions,
                             const std::vector<Page> &pages, QAction *about)
  : QToolBar(tr("Commands"), &window), window(window)
{
  setObjectName("commandRibbon");
  setMovable(false);
  setFloatable(false);
  setAllowedAreas(Qt::TopToolBarArea);
  setContextMenuPolicy(Qt::PreventContextMenu);
  toggleViewAction()->setEnabled(false);
  auto *panel = new QWidget(this);
  panel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
  auto *layout = new QVBoxLayout(panel);
  layout->setContentsMargins(4, 0, 4, 0);
  layout->setSpacing(0);
  auto *header = new QHBoxLayout;
  header->setSpacing(2);
  layout->addLayout(header);
  const auto button = [&window](QAction *action, QWidget *parent, bool quick = false) {
    window.addAction(action); // Keep window shortcuts active on every tab and in both layouts.
    auto *result = new QToolButton(parent);
    result->setDefaultAction(action);
    result->setToolButtonStyle(quick ? Qt::ToolButtonIconOnly : Qt::ToolButtonTextBesideIcon);
    result->setIconSize(QSize(16, 16));
    result->setAutoRaise(true);
    result->setAccessibleName(action->text().remove('&'));
    if(action->menu()) result->setPopupMode(QToolButton::InstantPopup);
    return result;
  };
  for(auto *action : quickActions) header->addWidget(button(action, panel, true));
  tabs = new QTabBar(panel);
  tabs->setObjectName("ribbonTabs");
  tabs->setAccessibleName(tr("Command categories"));
  tabs->setExpanding(false);
  tabs->setUsesScrollButtons(true);
  header->addWidget(tabs);
  header->addStretch();

  auto *menu = new QMenu(this);
  collapsed = menu->addAction(tr("Collapse ribbon"));
  collapsed->setObjectName("actionCollapse_ribbon");
  collapsed->setCheckable(true);
  collapsed->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_F1));
  collapsed->setToolTip(tr("Show only the tab row. Select a tab to expand the commands. Ctrl+F1 toggles this setting."));
  window.addAction(collapsed);
  collapseButton = button(collapsed, panel);
  collapseButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
  header->addWidget(collapseButton);
  classic = menu->addAction(tr("Classic menus"));
  classic->setObjectName("actionClassic_menus");
  classic->setCheckable(true);
  classic->setToolTip(tr("Use the traditional menu bar. Open, Save, Undo and Redo stay available below it."));
  menu->addSeparator();
  menu->addAction(about);
  auto *interface = new QToolButton(panel);
  interface->setObjectName("ribbonInterface");
  interface->setText(tr("Interface"));
  interface->setToolTip(tr("Choose classic menus or collapse the ribbon"));
  interface->setAutoRaise(true);
  interface->setMenu(menu);
  interface->setPopupMode(QToolButton::InstantPopup);
  header->addWidget(interface);

  body = new QStackedWidget(panel);
  body->setObjectName("ribbonPages");
  layout->addWidget(body);
  int bodyHeight = 0;
  for(const auto &page : pages) {
    tabs->addTab(page.title);
    auto *content = new QWidget;
    auto *row = new QHBoxLayout(content);
    row->setContentsMargins(0, 2, 0, 0);
    row->setSpacing(6);
    for(const auto &group : page.groups) {
      if(row->count()) {
        auto *line = new QFrame(content);
        line->setFrameShape(QFrame::VLine);
        row->addWidget(line);
      }
      auto *groupLayout = new QVBoxLayout;
      groupLayout->setSpacing(1);
      auto *commands = new QGridLayout;
      commands->setSpacing(1);
      int index = 0;
      int commandHeight = 0;
      for(auto *action : group.actions) {
        auto *tool = button(action, content);
        commandHeight = qMax(commandHeight, tool->sizeHint().height());
        commands->addWidget(tool, index % 2, index / 2, Qt::AlignLeft);
        ++index;
      }
      commands->setRowMinimumHeight(0, commandHeight);
      commands->setRowMinimumHeight(1, commandHeight);
      groupLayout->addLayout(commands, 1);
      auto *caption = new QLabel(group.title, content);
      caption->setAlignment(Qt::AlignCenter);
      auto font = caption->font();
      if(font.pointSizeF() > 0) font.setPointSizeF(font.pointSizeF() * 0.85);
      caption->setFont(font);
      groupLayout->addWidget(caption);
      row->addLayout(groupLayout);
    }
    row->addStretch();
    auto *scroll = new QScrollArea(body);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidgetResizable(true);
    scroll->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scroll->setWidget(content);
    bodyHeight = qMax(bodyHeight, content->sizeHint().height() + scroll->horizontalScrollBar()->sizeHint().height());
    body->addWidget(scroll);
  }
  body->setFixedHeight(bodyHeight);
  addWidget(panel);
  window.addToolBar(Qt::TopToolBarArea, this);

  QSettings settings;
  classic->setChecked(settings.value("interface/classicMenus", false).toBool());
  collapsed->setChecked(settings.value("interface/ribbonCollapsed", false).toBool());
  tabs->setCurrentIndex(qBound(0, settings.value("interface/ribbonTab", 0).toInt(), tabs->count() - 1));
  body->setCurrentIndex(tabs->currentIndex());
  connect(tabs, &QTabBar::currentChanged, this, [this](int index) {
    body->setCurrentIndex(index);
    collapsed->setChecked(false);
    QSettings().setValue("interface/ribbonTab", index);
  });
  connect(tabs, &QTabBar::tabBarClicked, this, [this](int index) { if(index >= 0) collapsed->setChecked(false); });
  connect(tabs, &QTabBar::tabBarDoubleClicked, this, [this](int index) { if(index >= 0) collapsed->toggle(); });
  connect(classic, &QAction::toggled, this, [this](bool checked) {
    QSettings().setValue("interface/classicMenus", checked);
    updateDisplay();
  });
  connect(collapsed, &QAction::toggled, this, [this](bool checked) {
    QSettings().setValue("interface/ribbonCollapsed", checked);
    updateDisplay();
  });
  updateDisplay();
}

void CommandRibbon::updateDisplay()
{
  window.menuBar()->setVisible(classic->isChecked());
  tabs->setVisible(!classic->isChecked());
  collapseButton->setVisible(!classic->isChecked());
  collapsed->setEnabled(!classic->isChecked());
  collapseButton->setArrowType(collapsed->isChecked() ? Qt::DownArrow : Qt::UpArrow);
  body->setVisible(!classic->isChecked() && !collapsed->isChecked());
}
