#pragma once

#include <QStyledItemDelegate>

// Format coordinate cells without losing numeric values used for sorting.
class StationListDelegate : public QStyledItemDelegate
{
public:
  enum { PeriodCount = Qt::UserRole + 1, Coordinates, CoordinateReference, CoordinatePrecision };
  using QStyledItemDelegate::QStyledItemDelegate;
protected:
  void initStyleOption(QStyleOptionViewItem *option, const QModelIndex &index) const override;
};
