#include "include/StationListDelegate.h"
#include <QApplication>
#include <QPainter>

QSize StationListDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &) const
{
  return {240, option.fontMetrics.height() * 2 + 10};
}

void StationListDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
  QStyleOptionViewItem item(option); initStyleOption(&item, index);
  const auto *style = item.widget ? item.widget->style() : QApplication::style();
  const QRect textRect = style->subElementRect(QStyle::SE_ItemViewItemText, &item, item.widget).adjusted(3, 3, -5, -3);
  item.text.clear(); style->drawControl(QStyle::CE_ItemViewItem, &item, painter, item.widget);
  painter->save();
  const auto foreground = item.state & QStyle::State_Selected ? QPalette::HighlightedText : QPalette::Text;
  painter->setPen(item.palette.color(foreground));
  QFont detailFont(item.font);
  if(detailFont.pointSizeF() > 0) detailFont.setPointSizeF(detailFont.pointSizeF() * .9);
  const QFontMetrics details(detailFont);
  const QString count = tr("%1 periods").arg(index.data(PeriodCount).toULongLong());
  const int countWidth = details.horizontalAdvance(count);
  QRect first(textRect.left(), textRect.top(), textRect.width(), option.fontMetrics.height());
  painter->setFont(detailFont);
  painter->drawText(first, Qt::AlignRight | Qt::AlignVCenter, count);
  QFont nameFont(item.font); nameFont.setBold(true); painter->setFont(nameFont);
  first.setRight(first.right() - countWidth - 12);
  painter->drawText(first, Qt::AlignLeft | Qt::AlignVCenter,
                    QFontMetrics(nameFont).elidedText(index.data().toString(), Qt::ElideRight, qMax(0, first.width())));
  painter->setFont(detailFont);
  const QRect second(textRect.left(), first.bottom() + 1, textRect.width(), details.height() + 2);
  painter->drawText(second, Qt::AlignLeft | Qt::AlignVCenter,
                    details.elidedText(index.data(Coordinates).toString(), Qt::ElideRight, second.width()));
  painter->restore();
}
