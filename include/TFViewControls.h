#pragma once

#include "CommandRibbon.h"
#include "TFViewStyle.h"
#include <functional>

class QComboBox;
class QCheckBox;
class QLineEdit;
class QDoubleSpinBox;

// Inline ribbon editors, with no modal configuration windows.
class TFViewControls : public QObject
{
public:
  TFViewControls(QWidget *owner, std::function<TFViewRange(int)> range,
                 std::function<void(int, TFViewRange)> applyRange,
                 std::function<void()> appearanceChanged);
  std::vector<CommandRibbon::Group> groups() const;
  const TFViewStyle &style() const { return appearance; }
  void setStyle(TFViewStyle style);
  void refreshRange();

private:
  void refreshColor();
  TFViewStyle appearance;
  std::function<TFViewRange(int)> range;
  std::function<void(int, TFViewRange)> applyRange;
  std::function<void()> appearanceChanged;
  QWidget *ranges, *colors, *sizes, *tensor;
  QComboBox *axis, *colorTarget, *color, *colorMap;
  QCheckBox *automatic, *reverse;
  QLineEdit *lower, *upper;
  std::array<QDoubleSpinBox *, 4> sizeControls;
  QDoubleSpinBox *colorMin, *colorMax;
};
