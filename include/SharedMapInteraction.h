#pragma once

#include "qcustomplot.h"
#include <QPointer>

// A view temporarily changes the survey map's interactions while its tab is active.
// Repeated activation preserves the original settings; leaving or destroying the
// view restores them. The map may be destroyed before the view during shutdown.
class SharedMapInteraction {
public:
  explicit SharedMapInteraction(QCustomPlot *plot) : plot(plot) {}
  ~SharedMapInteraction() { deactivate(); }
  SharedMapInteraction(const SharedMapInteraction &) = delete;
  SharedMapInteraction &operator=(const SharedMapInteraction &) = delete;

  bool active() const { return engaged && plot; }
  bool activate() {
    if (!plot)
      return false;
    if (!engaged) {
      interactions = plot->interactions();
      cursor = plot->cursor();
      engaged = true;
    }
    return true;
  }
  void configure(QCP::Interactions disabled = {}, const QCursor *replacement = nullptr) {
    if (!active())
      return;
    plot->setInteractions(interactions & ~disabled);
    plot->setCursor(replacement ? *replacement : cursor);
  }
  void deactivate() {
    if (active()) {
      configure();
      plot->replot(QCustomPlot::rpQueuedReplot);
    }
    engaged = false;
  }

private:
  QPointer<QCustomPlot> plot;
  bool engaged = false;
  QCP::Interactions interactions;
  QCursor cursor;
};
