#include "include/StationListDelegate.h"

void StationListDelegate::initStyleOption(QStyleOptionViewItem *option, const QModelIndex &index) const
{
  QStyledItemDelegate::initStyleOption(option, index);
  if(index.column() >= 2 && index.data().metaType().id() == QMetaType::Double)
    option->text = QString::number(index.data().toDouble(), 'f', index.data(CoordinatePrecision).toInt());
}
