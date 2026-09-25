#include "include/PlotColorMaps.h"

namespace PlotColorMaps {
const QStringList &names()
{
  static const QStringList values{"Viridis", "Thermal", "Jet", "Grayscale", "Polar"};
  return values;
}

QCPColorGradient gradient(const QString &name)
{
  if(name == "Thermal") return QCPColorGradient(QCPColorGradient::gpThermal);
  if(name == "Jet") return QCPColorGradient(QCPColorGradient::gpJet);
  if(name == "Grayscale") return QCPColorGradient(QCPColorGradient::gpGrayscale);
  if(name == "Polar") return QCPColorGradient(QCPColorGradient::gpPolar);
  QCPColorGradient colors;
  colors.setColorStopAt(0., QColor("#440154"));
  colors.setColorStopAt(.25, QColor("#3b528b"));
  colors.setColorStopAt(.5, QColor("#21918c"));
  colors.setColorStopAt(.75, QColor("#5ec962"));
  colors.setColorStopAt(1., QColor("#fde725"));
  return colors;
}
}
