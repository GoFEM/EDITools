#include "include/CommandIcons.h"
#include <QApplication>
#include <QIconEngine>
#include <QPainter>
#include <QPalette>

namespace {
class CommandIconEngine final : public QIconEngine
{
public:
  explicit CommandIconEngine(const QString &path) : glyph(path) {}
  QIconEngine *clone() const override { return new CommandIconEngine(*this); }
  bool isNull() override { return glyph.isNull(); }
  QPixmap scaledPixmap(const QSize &size, QIcon::Mode mode, QIcon::State state, qreal scale) override
  {
    auto result = glyph.pixmap(size, scale, QIcon::Normal, state);
    if(result.isNull()) return result;
    const auto group = mode == QIcon::Disabled ? QPalette::Disabled : QPalette::Active;
    const auto role = mode == QIcon::Selected ? QPalette::HighlightedText : QPalette::ButtonText;
    QPainter tint(&result);
    tint.setCompositionMode(QPainter::CompositionMode_SourceIn);
    tint.fillRect(result.rect(), QApplication::palette().color(group, role));
    return result;
  }
  QPixmap pixmap(const QSize &size, QIcon::Mode mode, QIcon::State state) override
  {
    return scaledPixmap(size, mode, state, 1);
  }
  void paint(QPainter *painter, const QRect &rect, QIcon::Mode mode, QIcon::State state) override
  {
    painter->drawPixmap(rect, scaledPixmap(rect.size(), mode, state, painter->device()->devicePixelRatioF()));
  }
private:
  QIcon glyph;
};
}

QIcon commandIcon(const char *name)
{
  // Explicit initialization also retains the resource object when linking the UI library.
  static const bool initialized = [] { Q_INIT_RESOURCE(icons); return true; }();
  Q_UNUSED(initialized);
  return QIcon(new CommandIconEngine(QStringLiteral(":/Icons/fluent/%1.svg").arg(QLatin1String(name))));
}
