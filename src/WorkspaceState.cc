#include "include/WorkspaceView.h"
#include "qcustomplot.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QJsonArray>
#include <QLineEdit>
#include <QListWidget>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTreeWidget>
#include <cmath>

namespace {
QList<QCustomPlot *> plotsFor(QWidget *content) {
  auto plots = content->findChildren<QCustomPlot *>();
  if (auto *view = dynamic_cast<WorkspaceView *>(content))
    for (auto *map : view->centralMaps())
      if (auto *plot = qobject_cast<QCustomPlot *>(map))
        if (!plots.contains(plot))
          plots.append(plot);
  return plots;
}

// Lists and station/response trees share the same root-row checkbox state format.
QList<QAbstractItemView *> selectionViews(const QWidget *content) {
  QList<QAbstractItemView *> result;
  for (auto *view : content->findChildren<QAbstractItemView *>())
    if (!view->objectName().isEmpty() &&
        (qobject_cast<QListWidget *>(view) || qobject_cast<QTreeWidget *>(view)))
      result.append(view);
  return result;
}
QString selectionKey(const QModelIndex &index) {
  const auto id = index.data(Qt::UserRole);
  return id.isValid() ? id.toString() : index.data(Qt::DisplayRole).toString();
}
QJsonObject saveSelections(const QWidget *content) {
  QJsonObject selections;
  for (auto *view : selectionViews(content)) {
    QJsonObject checked;
    const auto *model = view->model();
    for (int row = 0; row < model->rowCount(); ++row) {
      const auto index = model->index(row, 0);
      checked[selectionKey(index)] = index.data(Qt::CheckStateRole).toInt() == Qt::Checked;
    }
    selections[view->objectName()] = checked;
  }
  return selections;
}
void restoreSelections(QWidget *content, const QJsonObject &selections) {
  for (auto *view : selectionViews(content)) {
    if (!selections.contains(view->objectName()))
      continue;
    const QSignalBlocker block(view);
    const auto checked = selections[view->objectName()].toObject();
    auto *model = view->model();
    for (int row = 0; row < model->rowCount(); ++row) {
      const auto index = model->index(row, 0);
      const auto key = selectionKey(index);
      if (checked.contains(key))
        model->setData(index, checked[key].toBool() ? Qt::Checked : Qt::Unchecked, Qt::CheckStateRole);
    }
  }
}
} // namespace

QJsonObject WorkspaceView::savePlotRanges(QWidget *content) {
  QJsonObject result;
  for (auto *plot : plotsFor(content))
    if (!plot->objectName().isEmpty()) {
      QJsonArray rects;
      for (auto *rect : plot->axisRects()) {
        const auto x = rect->axis(QCPAxis::atBottom)->range(), y = rect->axis(QCPAxis::atLeft)->range();
        rects.append(QJsonArray{x.lower, x.upper, y.lower, y.upper});
      }
      result[plot->objectName()] = rects;
    }
  return result;
}

void WorkspaceView::restorePlotRanges(QWidget *content, const QJsonObject &state) {
  for (auto *plot : plotsFor(content)) {
    const auto rects = state.value(plot->objectName()).toArray();
    for (int i = 0; i < rects.size() && i < plot->axisRectCount(); ++i) {
      const auto values = rects[i].toArray();
      if (values.size() != 4)
        continue;
      for (int j = 0; j < 2; ++j) {
        const double lo = values[j * 2].toDouble(), hi = values[j * 2 + 1].toDouble();
        auto *axis = plot->axisRect(i)->axis(j == 0 ? QCPAxis::atBottom : QCPAxis::atLeft);
        if (std::isfinite(lo) && std::isfinite(hi) && lo < hi && QCPRange::validRange(lo, hi) &&
            (axis->scaleType() != QCPAxis::stLogarithmic || lo > 0))
          axis->setRange(lo, hi);
      }
    }
    plot->replot(QCustomPlot::rpQueuedReplot);
  }
}

QJsonObject WorkspaceView::saveViewState() const {
  QJsonObject values;
  for (const auto &entry : editors)
    if (auto *widget = entry.data()) {
      QJsonObject value;
      if (auto *c = qobject_cast<QComboBox *>(widget))
        value = {{"index", c->currentIndex()},
                 {"text", c->currentText()},
                 {"data", QJsonValue::fromVariant(c->currentData())},
                 {"hasData", c->currentData().isValid()}};
      else if (auto *s = qobject_cast<QDoubleSpinBox *>(widget))
        value["value"] = s->value();
      else if (auto *s = qobject_cast<QSpinBox *>(widget))
        value["value"] = s->value();
      else if (auto *c = qobject_cast<QCheckBox *>(widget))
        value["checked"] = c->isChecked();
      else if (auto *e = qobject_cast<QLineEdit *>(widget))
        value["text"] = e->text();
      if (!value.isEmpty())
        values[widget->objectName()] = value;
    }
  QJsonObject splitters, tabs;
  for (auto *split : findChildren<QSplitter *>())
    if (!split->objectName().isEmpty())
      splitters[split->objectName()] = QString::fromLatin1(split->saveState().toBase64());
  for (auto *tab : findChildren<QTabWidget *>())
    if (!tab->objectName().isEmpty())
      tabs[tab->objectName()] = tab->currentIndex();
  return {{"controls", values},
          {"plots", savePlotRanges(const_cast<WorkspaceView *>(this))},
          {"splitters", splitters},
          {"tabs", tabs},
          {"selections", saveSelections(this)}};
}

void WorkspaceView::restoreViewState(const QJsonObject &state) {
  const auto values = state.value("controls").toObject();
  std::vector<std::unique_ptr<QSignalBlocker>> blockers;
  for (const auto &entry : editors)
    if (entry)
      blockers.emplace_back(new QSignalBlocker(entry));
  for (const auto &entry : editors)
    if (auto *widget = entry.data()) {
      if (!values.contains(widget->objectName()))
        continue;
      const auto value = values.value(widget->objectName()).toObject();
      if (auto *c = qobject_cast<QComboBox *>(widget)) {
        int index = value.value("hasData").toBool() ? c->findData(value.value("data").toVariant())
                                                    : c->findText(value.value("text").toString());
        if (index < 0 && value.value("hasData").toBool()) {
          c->addItem(tr("Unavailable: %1").arg(value.value("text").toString()),
                     value.value("data").toVariant());
          index = c->count() - 1;
        }
        if (index >= 0)
          c->setCurrentIndex(index);
      } else if (auto *s = qobject_cast<QDoubleSpinBox *>(widget)) {
        if (value.value("value").isDouble())
          s->setValue(value.value("value").toDouble());
      } else if (auto *s = qobject_cast<QSpinBox *>(widget)) {
        if (value.value("value").isDouble())
          s->setValue(value.value("value").toInt());
      } else if (auto *c = qobject_cast<QCheckBox *>(widget))
        c->setChecked(value.value("checked").toBool());
      else if (auto *e = qobject_cast<QLineEdit *>(widget))
        e->setText(value.value("text").toString());
    }
  blockers.clear();
  restoreSelections(this, state["selections"].toObject());
  if (refreshState)
    refreshState();
  const auto splits = state.value("splitters").toObject(), tabs = state.value("tabs").toObject();
  for (auto *split : findChildren<QSplitter *>())
    if (splits.contains(split->objectName()))
      split->restoreState(QByteArray::fromBase64(splits.value(split->objectName()).toString().toLatin1()));
  for (auto *tab : findChildren<QTabWidget *>())
    if (tabs.contains(tab->objectName()))
      tab->setCurrentIndex(tabs.value(tab->objectName()).toInt());
  restorePlotRanges(this, state.value("plots").toObject());
}
