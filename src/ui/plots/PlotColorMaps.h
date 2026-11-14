#ifndef EDITOOLS_PLOT_COLOR_MAPS_H
#define EDITOOLS_PLOT_COLOR_MAPS_H

#include "qcustomplot.h"

namespace PlotColorMaps {
const QStringList &names();
// Unknown names use the default Viridis palette.
QCPColorGradient gradient(const QString &name = "Viridis");
}
#endif
