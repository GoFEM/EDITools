#ifndef EDITOOLS_TOOL_TIPS_H
#define EDITOOLS_TOOL_TIPS_H

#include <QLabel>

namespace ToolTips {
inline void set(QWidget *control, const QString &text)
{
  control->setToolTip("<qt>" + text.toHtmlEscaped().replace('\n', "<br>") + "</qt>");
  control->setAccessibleDescription(text);
}

// Keep the same explanation on the input and its label.
inline QLabel *label(QWidget *control, const QString &text, const QString &explanation)
{
  auto *label = new QLabel(text, control->parentWidget());
  label->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
  label->setBuddy(control);
  set(label, explanation);
  set(control, explanation);
  return label;
}
}
#endif
