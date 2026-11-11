#include "include/StationListDelegate.h"

QSize StationListDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &) const
{
  return {240, option.fontMetrics.height() * 2 + 10};
}

void StationListDelegate::initStyleOption(QStyleOptionViewItem *option, const QModelIndex &index) const
{
  QStyledItemDelegate::initStyleOption(option, index);
  option->text = tr("%1 · %2 periods\n%3").arg(index.data().toString())
    .arg(index.data(PeriodCount).toULongLong()).arg(index.data(Coordinates).toString());
  option->text.replace(QLatin1Char('\n'), QChar::LineSeparator);
  // Let the platform style choose text/background colors together, including
  // inactive selections and high-contrast themes.
}
