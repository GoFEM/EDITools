#pragma once
#include "include/AutoMask.h"
#include <QWidget>
#include <QStringList>
#include <functional>
class QCheckBox;
class QDoubleSpinBox;
class QSpinBox;
class WorkspaceView;

// Shared settings editor for station workflow and detailed batch review.
class AutoMaskSettings : public QWidget {
  Q_OBJECT
public:
  explicit AutoMaskSettings(const AutoMask::Options &options, QWidget *parent = nullptr);
  AutoMask::Options options() const;
  void setOptions(const AutoMask::Options &options);
  void moveToRibbon(WorkspaceView &view, const QString &page);
  static QStringList criterionNames();
  static QStringList criterionHelp();
signals:
  void changed();
private:
  std::function<void()> updateVisibility;
  AutoMask::Options initial;
  std::array<QCheckBox *, 4> components;
  std::array<QCheckBox *, AutoMask::criterionCount> criteria;
  std::array<QDoubleSpinBox *, AutoMask::criterionCount> thresholds;
  QSpinBox *votes, *neighbours, *density;
  QDoubleSpinBox *gap, *radius, *tolerance, *departure, *agreement, *rhoError, *phaseError,
    *tipperError, *tipperDeparture, *tipperAgreement;
};
bool editAutoMaskSettings(QWidget *parent, AutoMask::Options &options);
