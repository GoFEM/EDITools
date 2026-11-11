#pragma once

#include "CommandRibbon.h"
#include "SurveyCoordinates.h"
#include <QJsonObject>
#include <QPointer>
#include <QWidget>
#include <functional>
#include <map>
#include <memory>
#include <string>

class QLayout;
class QAbstractButton;
class QTabWidget;
class QComboBox;
class MTSurveyData;

// A view owns its content and its ribbon editors even when those editors are
// displayed by the main window. Each instance has independent display state.
class WorkspaceView : public QWidget {
  Q_OBJECT
public:
  explicit WorkspaceView(QWidget *parent = nullptr);
  ~WorkspaceView() override;
  using Responses = std::map<std::string, MTSurveyData>;
  // Refresh only dirty views, preserving zoom unless the coordinate frame changed.
  void markSurveyDirty() { surveyDirty = true; }
  void refreshSurvey(const std::shared_ptr<MTSurveyData> &survey, const Responses &responses);
  virtual void selectStation(const QString &) {}
  virtual void setSharedMapActive(bool, const SurveyCoordinates &) {}
  // Asynchronous views call this once their refreshed curves are ready.
  void restorePendingPlotRanges();
  void setEmbedded();
  void setCentralMaps(const QList<QWidget *> &maps, const QStringList &labels, const QString &toolPage);
  QWidget *centralMap() const;
  QList<QWidget *> centralMaps() const;
  virtual SurveyCoordinates centralMapCoordinates() const { return {}; }
  std::vector<CommandRibbon::Page> toolPages() const;
  void addToolGroup(const QString &page, const QString &title, QLayout *controls);
  void setToolVisible(QWidget *control, bool visible);
  void setStateRefresh(std::function<void()> refresh) { refreshState = std::move(refresh); }
  virtual QJsonObject saveViewState() const;
  virtual void restoreViewState(const QJsonObject &state);
  void setExportButton(const QString &format, QAbstractButton *button);
  QAbstractButton *exportButton(const QString &format) const;
  QString exportName(const QString &extension) const;
  QString displayTitle() const;
  template <class T> T *control(const QString &name) const {
    if (auto *found = findChild<T *>(name))
      return found;
    for (const auto &page : pages)
      if (page.widget)
        if (auto *found = page.widget->findChild<T *>(name))
          return found;
    for (const auto &map : sharedMaps)
      if (map) {
        if (map->objectName() == name)
          if (auto *found = qobject_cast<T *>(map.data()))
            return found;
        if (auto *found = map->findChild<T *>(name))
          return found;
      }
    return nullptr;
  }
  static QJsonObject savePlotRanges(QWidget *content);
  static void restorePlotRanges(QWidget *content, const QJsonObject &state);
signals:
  void stationActivated(const QString &name);
  void toolsChanged();
  void centralMapChanged();

protected:
  virtual void refreshData(const std::shared_ptr<MTSurveyData> &, const Responses &) {}
  virtual bool hasDeferredPlots() const { return false; }
  void showEvent(QShowEvent *event) override;
  bool eventFilter(QObject *object, QEvent *event) override;

private:
  bool surveyDirty = true;
  bool surveyInitialized = false;
  QString surveyCoordinates;
  QJsonObject pendingPlotRanges;
  QList<QPointer<QWidget>> sharedMaps;
  QPointer<QComboBox> sharedMapChoice;
  struct ToolPage {
    QString title;
    QPointer<QWidget> widget;
  };
  std::vector<ToolPage> pages;
  QList<QPointer<QWidget>> editors;
  QHash<QWidget *, QPointer<QWidget>> containers;
  QHash<QString, QPointer<QAbstractButton>> exports;
  std::function<void()> refreshState;
  QTabWidget *standaloneTools = nullptr;
  bool embedded = false;
};
