#include "ui/workspace/WorkspaceView.h"
#include "core/MTSurveyData.h"
#include <QAbstractButton>
#include <QComboBox>
#include <QEvent>
#include <QGridLayout>
#include <QLabel>
#include <QRegularExpression>
#include <QTabWidget>

WorkspaceView::WorkspaceView(QWidget *parent) : QWidget(parent) {}
WorkspaceView::~WorkspaceView() {
  // These widgets may be reparented to the main window while displayed there.
  for (const auto &map : sharedMaps)
    delete map.data();
  for (const auto &page : pages)
    delete page.widget.data();
}

void WorkspaceView::refreshSurvey(const std::shared_ptr<MTSurveyData> &survey, const Responses &responses) {
  if (!surveyDirty)
    return;
  const auto coordinates = survey ? QString::fromStdString(survey->coordinates().description()) : QString();
  const bool preserveRanges = surveyInitialized && surveyCoordinates == coordinates;
  const auto ranges = preserveRanges ? savePlotRanges(this) : QJsonObject();
  pendingPlotRanges = {};
  surveyDirty = false;
  refreshData(survey, responses);
  surveyCoordinates = coordinates;
  surveyInitialized = true;
  if (preserveRanges) {
    restorePlotRanges(this, ranges);
    if (hasDeferredPlots())
      pendingPlotRanges = ranges;
  }
}

void WorkspaceView::restorePendingPlotRanges() {
  if (pendingPlotRanges.isEmpty())
    return;
  restorePlotRanges(this, pendingPlotRanges);
  pendingPlotRanges = {};
}

void WorkspaceView::setCentralMaps(const QList<QWidget *> &maps, const QStringList &labels,
                                   const QString &toolPage) {
  for (auto *map : maps) {
    map->hide();
    map->setParent(this);
    map->setMinimumSize(0, 0);
    sharedMaps.append(map);
  }
  if (maps.size() > 1) {
    auto *row = new QHBoxLayout;
    auto *choice = new QComboBox(this);
    sharedMapChoice = choice;
    choice->setObjectName("sharedMapChoice");
    choice->addItems(labels);
    row->addWidget(new QLabel(tr("Map:"), this));
    row->addWidget(choice);
    connect(choice, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
            [this] { emit centralMapChanged(); });
    addToolGroup(toolPage, tr("Shared map"), row);
  }
}
QWidget *WorkspaceView::centralMap() const {
  const int index = sharedMapChoice ? sharedMapChoice->currentIndex() : 0;
  return index >= 0 && index < sharedMaps.size() ? sharedMaps[index].data() : nullptr;
}
QList<QWidget *> WorkspaceView::centralMaps() const {
  QList<QWidget *> result;
  for (const auto &map : sharedMaps)
    if (map)
      result.append(map);
  return result;
}

void WorkspaceView::setEmbedded() {
  embedded = true;
  if (standaloneTools)
    standaloneTools->hide();
}

std::vector<CommandRibbon::Page> WorkspaceView::toolPages() const {
  std::vector<CommandRibbon::Page> result;
  for (const auto &page : pages)
    if (page.widget)
      result.push_back({page.title, {{QString(), {}, page.widget}}});
  return result;
}

void WorkspaceView::addToolGroup(const QString &pageName, const QString &title, QLayout *controls) {
  QWidget *page = nullptr;
  for (const auto &candidate : pages)
    if (candidate.title == pageName)
      page = candidate.widget;
  if (!page) {
    page = new QWidget(this);
    page->hide();
    auto *row = new QHBoxLayout(page);
    row->setContentsMargins(0, 0, 0, 0);
    row->setSpacing(8);
    row->addStretch();
    pages.push_back({pageName, page});
  }
  if (auto *parent = qobject_cast<QLayout *>(controls->parent()))
    parent->removeItem(controls);
  QList<QWidget *> widgets;
  const auto collect = [&](const auto &self, QLayout *layout) -> void {
    while (auto *item = layout->takeAt(0)) {
      if (item->widget())
        widgets.push_back(item->widget());
      if (item->layout()) {
        self(self, item->layout());
        delete item->layout();
      } else
        delete item;
    }
  };
  collect(collect, controls);
  delete controls;
  auto *group = new QWidget(page);
  auto *column = new QVBoxLayout(group);
  column->setContentsMargins(0, 0, 0, 0);
  column->setSpacing(2);
  auto *grid = new QGridLayout;
  grid->setSpacing(3);
  column->addLayout(grid, 1);
  int position = 0;
  for (int i = 0; i < widgets.size(); ++i) {
    auto *unit = new QWidget(group);
    auto *pair = new QHBoxLayout(unit);
    pair->setContentsMargins(0, 0, 0, 0);
    pair->setSpacing(3);
    auto *widget = widgets[i];
    if (qobject_cast<QLabel *>(widget) && i + 1 < widgets.size() && !qobject_cast<QLabel *>(widgets[i + 1])) {
      pair->addWidget(widget);
      containers[widget] = unit;
      widget = widgets[++i];
    }
    if (auto *combo = qobject_cast<QComboBox *>(widget)) {
      combo->setMinimumWidth(0);
      combo->setMaximumWidth(220);
      combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
      combo->setMinimumContentsLength(10);
    }
    pair->addWidget(widget);
    widget->show();
    containers[widget] = unit;
    if (!widget->objectName().isEmpty())
      editors.push_back(widget);
    grid->addWidget(unit, position % 2, position / 2, Qt::AlignLeft);
    ++position;
  }
  auto *caption = new QLabel(title, group);
  caption->setAlignment(Qt::AlignCenter);
  auto font = caption->font();
  if (font.pointSizeF() > 0)
    font.setPointSizeF(font.pointSizeF() * .85);
  caption->setFont(font);
  column->addWidget(caption);
  auto *row = qobject_cast<QHBoxLayout *>(page->layout());
  if (row->count() > 1) {
    auto *line = new QFrame(page);
    line->setFrameShape(QFrame::VLine);
    row->insertWidget(row->count() - 1, line);
  }
  row->insertWidget(row->count() - 1, group);
}

void WorkspaceView::setToolVisible(QWidget *control, bool visible) {
  if (auto *container = containers.value(control).data())
    container->setVisible(visible);
  else if (control)
    control->setVisible(visible);
}

void WorkspaceView::showEvent(QShowEvent *event) {
  if (!embedded && !pages.empty() && !standaloneTools) {
    standaloneTools = new QTabWidget(this);
    standaloneTools->setMaximumHeight(120);
    for (const auto &page : pages)
      if (page.widget)
        standaloneTools->addTab(page.widget, page.title);
    if (auto *box = qobject_cast<QBoxLayout *>(layout()))
      box->insertWidget(0, standaloneTools);
  }
  QWidget::showEvent(event);
}

void WorkspaceView::setExportButton(const QString &format, QAbstractButton *button) {
  exports[format] = button;
  button->installEventFilter(this);
}
QAbstractButton *WorkspaceView::exportButton(const QString &format) const { return exports.value(format); }
bool WorkspaceView::eventFilter(QObject *object, QEvent *event) {
  if (event->type() == QEvent::EnabledChange)
    emit toolsChanged();
  return QWidget::eventFilter(object, event);
}
QString WorkspaceView::displayTitle() const {
  const auto title = property("viewTitle").toString();
  return title.isEmpty() ? windowTitle() : title;
}
QString WorkspaceView::exportName(const QString &extension) const {
  auto name = displayTitle();
  name.replace(QRegularExpression("[^\\w .-]"), "_");
  return name + "." + extension;
}
