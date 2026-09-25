#ifndef EDITOOLS_HELP_BUTTON_H
#define EDITOOLS_HELP_BUTTON_H

#include <QHBoxLayout>
#include <QLabel>
#include <QToolButton>
#include <QToolTip>

namespace UiHelp {
// Hover, click, or focus and press Space to read the same contextual explanation.
inline QToolButton *button(QWidget *parent, const QString &topic, const QString &text, const QString &name = {})
{
  auto *help = new QToolButton(parent);
  help->setObjectName(name); help->setText("?"); help->setAutoRaise(true);
  help->setFixedSize(22, 22); help->setFocusPolicy(Qt::StrongFocus);
  help->setAccessibleName(QObject::tr("Help: %1").arg(topic)); help->setAccessibleDescription(text);
  help->setToolTip("<qt>" + text.toHtmlEscaped().replace('\n', "<br>") + "</qt>");
  QObject::connect(help, &QToolButton::clicked, help, [help] {
    QToolTip::showText(help->mapToGlobal(QPoint(0, help->height())), help->toolTip(), help);
  });
  return help;
}
inline QWidget *label(QWidget *parent, const QString &text, const QString &explanation, const QString &name = {})
{
  auto *widget = new QWidget(parent);
  widget->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
  auto *row = new QHBoxLayout(widget); row->setContentsMargins(0, 0, 0, 0); row->setSpacing(2);
  row->addWidget(new QLabel(text, widget)); row->addWidget(button(widget, text, explanation, name));
  return widget;
}
}
#endif
