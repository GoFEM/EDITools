#pragma once
#include <QTabWidget>
#include <QJsonObject>

class WorkspaceTabs : public QTabWidget
{
  Q_OBJECT
public:
  explicit WorkspaceTabs(QWidget *parent = nullptr);
  void addView(QWidget *view, const QString &type, const QString &title, const QString &id = {});
  QString id(QWidget *view) const;
  QString type(QWidget *view) const;
  void renameCurrent();
  void closeCurrent();
  void clearViews();
signals:
  void closing(QWidget *view);
  void duplicateRequested(QWidget *view);
  void reopenRequested();
  void titleChanged(QWidget *view, const QString &title);
private:
  void renameTab(int index);
  void closeTab(int index);
};
