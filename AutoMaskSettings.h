#pragma once
#include "include/AutoMask.h"
#include <QWidget>
#include <QStringList>
class QCheckBox;
class QDoubleSpinBox;
class QSpinBox;

// Shared settings editor for station workflow and detailed batch review.
class AutoMaskSettings : public QWidget {
  Q_OBJECT
public:
  explicit AutoMaskSettings(const AutoMask::Options &options, QWidget *parent = nullptr);
  AutoMask::Options options() const;
  static QStringList criterionNames();
  static QStringList criterionHelp();
signals:
  void changed();
private:
  AutoMask::Options initial;
  std::array<QCheckBox *, 4> components;
  std::array<QCheckBox *, AutoMask::criterionCount> criteria;
  std::array<QDoubleSpinBox *, AutoMask::criterionCount> thresholds;
  QSpinBox *votes, *neighbours, *density;
  QDoubleSpinBox *gap, *radius, *tolerance, *departure, *agreement, *rhoError, *phaseError,
    *tipperError, *tipperDeparture, *tipperAgreement;
};
bool editAutoMaskSettings(QWidget *parent, AutoMask::Options &options);
