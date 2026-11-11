#pragma once

#include "CommandRibbon.h"
#include <QJsonObject>
#include <QPointer>
#include <QWidget>
#include <functional>

class QLayout;
class QAbstractButton;
class QTabWidget;

// A view owns its content and its ribbon editors even when those editors are
// displayed by the main window. Each instance has independent display state.
class WorkspaceView : public QWidget
{
  Q_OBJECT
public:
  explicit WorkspaceView(QWidget *parent = nullptr);
  ~WorkspaceView() override;
  void setEmbedded();
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
  template<class T> T *control(const QString &name) const {
    if(auto *found = findChild<T *>(name)) return found;
    for(const auto &page: pages) if(page.widget)
      if(auto *found = page.widget->findChild<T *>(name)) return found;
    return nullptr;
  }
  static QJsonObject savePlotRanges(QWidget *content);
  static void restorePlotRanges(QWidget *content, const QJsonObject &state);
signals:
  void stationActivated(const QString &name);
  void toolsChanged();
protected:
  void showEvent(QShowEvent *event) override;
  bool eventFilter(QObject *object, QEvent *event) override;
private:
  struct ToolPage { QString title; QPointer<QWidget> widget; };
  std::vector<ToolPage> pages;
  QList<QPointer<QWidget>> editors;
  QHash<QWidget *, QPointer<QWidget>> containers;
  QHash<QString, QPointer<QAbstractButton>> exports;
  std::function<void()> refreshState;
  QTabWidget *standaloneTools = nullptr;
  bool embedded = false;
};
