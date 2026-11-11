#include "include/WorkspaceView.h"
#include "qcustomplot.h"
#include <QAbstractButton>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QJsonArray>
#include <QLabel>
#include <QLineEdit>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTreeWidget>
#include <QListWidget>
#include <cmath>
#include <memory>
#include <QEvent>
#include <QRegularExpression>

WorkspaceView::WorkspaceView(QWidget *parent): QWidget(parent) {}
WorkspaceView::~WorkspaceView() { for(const auto &page: pages) delete page.widget.data(); }

void WorkspaceView::setEmbedded()
{
  embedded = true;
  if(standaloneTools) standaloneTools->hide();
}

std::vector<CommandRibbon::Page> WorkspaceView::toolPages() const
{
  std::vector<CommandRibbon::Page> result;
  for(const auto &page: pages) if(page.widget) result.push_back({page.title, {{QString(), {}, page.widget}}});
  return result;
}

void WorkspaceView::addToolGroup(const QString &pageName, const QString &title, QLayout *controls)
{
  QWidget *page = nullptr;
  for(const auto &candidate: pages) if(candidate.title == pageName) page = candidate.widget;
  if(!page) {
    page = new QWidget(this); page->hide();
    auto *row = new QHBoxLayout(page); row->setContentsMargins(0, 0, 0, 0); row->setSpacing(8);
    row->addStretch(); pages.push_back({pageName, page});
  }
  if(auto *parent = qobject_cast<QLayout *>(controls->parent())) parent->removeItem(controls);
  QList<QWidget *> widgets;
  const auto collect = [&](const auto &self, QLayout *layout) -> void {
    while(auto *item = layout->takeAt(0)) {
      if(item->widget()) widgets.push_back(item->widget());
      if(item->layout()) { self(self, item->layout()); delete item->layout(); }
      else delete item;
    }
  };
  collect(collect, controls); delete controls;
  auto *group = new QWidget(page); auto *column = new QVBoxLayout(group);
  column->setContentsMargins(0, 0, 0, 0); column->setSpacing(2);
  auto *grid = new QGridLayout; grid->setSpacing(3); column->addLayout(grid, 1);
  int position = 0;
  for(int i = 0; i < widgets.size(); ++i) {
    auto *unit = new QWidget(group); auto *pair = new QHBoxLayout(unit);
    pair->setContentsMargins(0, 0, 0, 0); pair->setSpacing(3);
    auto *widget = widgets[i];
    if(qobject_cast<QLabel *>(widget) && i + 1 < widgets.size() && !qobject_cast<QLabel *>(widgets[i + 1])) {
      pair->addWidget(widget); containers[widget] = unit; widget = widgets[++i];
    }
    if(auto *combo = qobject_cast<QComboBox *>(widget)) {
      combo->setMinimumWidth(0); combo->setMaximumWidth(220);
      combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon); combo->setMinimumContentsLength(10);
    }
    pair->addWidget(widget); widget->show(); containers[widget] = unit;
    if(!widget->objectName().isEmpty()) editors.push_back(widget);
    grid->addWidget(unit, position % 2, position / 2, Qt::AlignLeft); ++position;
  }
  auto *caption = new QLabel(title, group); caption->setAlignment(Qt::AlignCenter);
  auto font = caption->font(); if(font.pointSizeF() > 0) font.setPointSizeF(font.pointSizeF() * .85);
  caption->setFont(font); column->addWidget(caption);
  auto *row = qobject_cast<QHBoxLayout *>(page->layout());
  if(row->count() > 1) { auto *line = new QFrame(page); line->setFrameShape(QFrame::VLine); row->insertWidget(row->count() - 1, line); }
  row->insertWidget(row->count() - 1, group);
}

void WorkspaceView::setToolVisible(QWidget *control, bool visible)
{
  if(auto *container = containers.value(control).data()) container->setVisible(visible);
  else if(control) control->setVisible(visible);
}

void WorkspaceView::showEvent(QShowEvent *event)
{
  if(!embedded && !pages.empty() && !standaloneTools) {
    standaloneTools = new QTabWidget(this); standaloneTools->setMaximumHeight(120);
    for(const auto &page: pages) if(page.widget) standaloneTools->addTab(page.widget, page.title);
    if(auto *box = qobject_cast<QBoxLayout *>(layout())) box->insertWidget(0, standaloneTools);
  }
  QWidget::showEvent(event);
}

void WorkspaceView::setExportButton(const QString &format, QAbstractButton *button) { exports[format] = button; button->installEventFilter(this); }
QAbstractButton *WorkspaceView::exportButton(const QString &format) const { return exports.value(format); }
bool WorkspaceView::eventFilter(QObject *object, QEvent *event) {
  if(event->type() == QEvent::EnabledChange) emit toolsChanged();
  return QWidget::eventFilter(object, event);
}
QString WorkspaceView::displayTitle() const { const auto title = property("viewTitle").toString(); return title.isEmpty() ? windowTitle() : title; }
QString WorkspaceView::exportName(const QString &extension) const {
  auto name = displayTitle(); name.replace(QRegularExpression("[^\\w .-]"), "_"); return name + "." + extension;
}

QJsonObject WorkspaceView::savePlotRanges(QWidget *content)
{
  QJsonObject result;
  for(auto *plot: content->findChildren<QCustomPlot *>()) if(!plot->objectName().isEmpty()) {
    QJsonArray rects;
    for(auto *rect: plot->axisRects()) {
      const auto x = rect->axis(QCPAxis::atBottom)->range(), y = rect->axis(QCPAxis::atLeft)->range();
      rects.append(QJsonArray{x.lower, x.upper, y.lower, y.upper});
    }
    result[plot->objectName()] = rects;
  }
  return result;
}

void WorkspaceView::restorePlotRanges(QWidget *content, const QJsonObject &state)
{
  for(auto *plot: content->findChildren<QCustomPlot *>()) {
    const auto rects = state.value(plot->objectName()).toArray();
    for(int i = 0; i < rects.size() && i < plot->axisRectCount(); ++i) {
      const auto values = rects[i].toArray(); if(values.size() != 4) continue;
      for(int j = 0; j < 2; ++j) {
        const double lo = values[j * 2].toDouble(), hi = values[j * 2 + 1].toDouble();
        auto *axis = plot->axisRect(i)->axis(j == 0 ? QCPAxis::atBottom : QCPAxis::atLeft);
        if(std::isfinite(lo) && std::isfinite(hi) && lo < hi && QCPRange::validRange(lo, hi) &&
           (axis->scaleType() != QCPAxis::stLogarithmic || lo > 0)) axis->setRange(lo, hi);
      }
    }
    plot->replot(QCustomPlot::rpQueuedReplot);
  }
}

QJsonObject WorkspaceView::saveViewState() const
{
  QJsonObject values;
  for(const auto &entry: editors) if(auto *widget = entry.data()) {
    QJsonObject value;
    if(auto *c = qobject_cast<QComboBox *>(widget)) value = {{"index", c->currentIndex()}, {"text", c->currentText()},
      {"data", QJsonValue::fromVariant(c->currentData())}, {"hasData", c->currentData().isValid()}};
    else if(auto *s = qobject_cast<QDoubleSpinBox *>(widget)) value["value"] = s->value();
    else if(auto *s = qobject_cast<QSpinBox *>(widget)) value["value"] = s->value();
    else if(auto *c = qobject_cast<QCheckBox *>(widget)) value["checked"] = c->isChecked();
    else if(auto *e = qobject_cast<QLineEdit *>(widget)) value["text"] = e->text();
    if(!value.isEmpty()) values[widget->objectName()] = value;
  }
  QJsonObject splitters, tabs, selections;
  for(auto *list: findChildren<QListWidget *>()) if(!list->objectName().isEmpty()) {
    QJsonObject checked;
    for(int i = 0; i < list->count(); ++i) {
      const auto *item = list->item(i);
      const auto key = item->data(Qt::UserRole).isValid() ? item->data(Qt::UserRole).toString() : item->text();
      checked[key] = item->checkState() == Qt::Checked;
    }
    selections[list->objectName()] = checked;
  }
  for(auto *tree: findChildren<QTreeWidget *>()) if(!tree->objectName().isEmpty()) {
    QJsonObject checked;
    for(int i = 0; i < tree->topLevelItemCount(); ++i) {
      const auto *item = tree->topLevelItem(i);
      const auto key = item->data(0, Qt::UserRole).isValid() ? item->data(0, Qt::UserRole).toString() : item->text(0);
      checked[key] = item->checkState(0) == Qt::Checked;
    }
    selections[tree->objectName()] = checked;
  }
  for(auto *split: findChildren<QSplitter *>()) if(!split->objectName().isEmpty()) splitters[split->objectName()] = QString::fromLatin1(split->saveState().toBase64());
  for(auto *tab: findChildren<QTabWidget *>()) if(!tab->objectName().isEmpty()) tabs[tab->objectName()] = tab->currentIndex();
  return {{"controls", values}, {"plots", savePlotRanges(const_cast<WorkspaceView *>(this))}, {"splitters", splitters}, {"tabs", tabs}, {"selections", selections}};
}

void WorkspaceView::restoreViewState(const QJsonObject &state)
{
  const auto values = state.value("controls").toObject();
  std::vector<std::unique_ptr<QSignalBlocker>> blockers;
  for(const auto &entry: editors) if(entry) blockers.emplace_back(new QSignalBlocker(entry));
  for(const auto &entry: editors) if(auto *widget = entry.data()) {
    if(!values.contains(widget->objectName())) continue;
    const auto value = values.value(widget->objectName()).toObject();
    if(auto *c = qobject_cast<QComboBox *>(widget)) {
      int index = value.value("hasData").toBool() ? c->findData(value.value("data").toVariant()) : c->findText(value.value("text").toString());
      if(index < 0 && value.value("hasData").toBool()) {
        c->addItem(tr("Unavailable: %1").arg(value.value("text").toString()), value.value("data").toVariant()); index = c->count() - 1;
      }
      if(index >= 0) c->setCurrentIndex(index);
    } else if(auto *s = qobject_cast<QDoubleSpinBox *>(widget)) { if(value.value("value").isDouble()) s->setValue(value.value("value").toDouble()); }
    else if(auto *s = qobject_cast<QSpinBox *>(widget)) { if(value.value("value").isDouble()) s->setValue(value.value("value").toInt()); }
    else if(auto *c = qobject_cast<QCheckBox *>(widget)) c->setChecked(value.value("checked").toBool());
    else if(auto *e = qobject_cast<QLineEdit *>(widget)) e->setText(value.value("text").toString());
  }
  blockers.clear();
  const auto selections = state["selections"].toObject();
  for(auto *list: findChildren<QListWidget *>()) if(selections.contains(list->objectName())) {
    const QSignalBlocker block(list); const auto checked = selections[list->objectName()].toObject();
    for(int i = 0; i < list->count(); ++i) {
      auto *item = list->item(i); const auto key = item->data(Qt::UserRole).isValid() ? item->data(Qt::UserRole).toString() : item->text();
      if(checked.contains(key)) item->setCheckState(checked[key].toBool() ? Qt::Checked : Qt::Unchecked);
    }
  }
  for(auto *tree: findChildren<QTreeWidget *>()) if(selections.contains(tree->objectName())) {
    const QSignalBlocker block(tree); const auto checked = selections[tree->objectName()].toObject();
    for(int i = 0; i < tree->topLevelItemCount(); ++i) {
      auto *item = tree->topLevelItem(i); const auto key = item->data(0, Qt::UserRole).isValid() ? item->data(0, Qt::UserRole).toString() : item->text(0);
      if(checked.contains(key)) item->setCheckState(0, checked[key].toBool() ? Qt::Checked : Qt::Unchecked);
    }
  }
  if(refreshState) refreshState();
  const auto splits = state.value("splitters").toObject(), tabs = state.value("tabs").toObject();
  for(auto *split: findChildren<QSplitter *>()) if(splits.contains(split->objectName())) split->restoreState(QByteArray::fromBase64(splits.value(split->objectName()).toString().toLatin1()));
  for(auto *tab: findChildren<QTabWidget *>()) if(tabs.contains(tab->objectName())) tab->setCurrentIndex(tabs.value(tab->objectName()).toInt());
  restorePlotRanges(this, state.value("plots").toObject());
}
