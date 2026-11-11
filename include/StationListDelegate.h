#pragma once

#include <QStyledItemDelegate>

// Keep DisplayRole as the station key so sorting, selection and edits use names.
class StationListDelegate : public QStyledItemDelegate
{
public:
  enum { PeriodCount = Qt::UserRole + 1, Coordinates, CoordinateReference };
  using QStyledItemDelegate::QStyledItemDelegate;
  QSize sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const override;
protected:
  void initStyleOption(QStyleOptionViewItem *option, const QModelIndex &index) const override;
};
