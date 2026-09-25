#ifndef EDITOOLS_MAP_DISPLAY_SETTINGS_H
#define EDITOOLS_MAP_DISPLAY_SETTINGS_H

#include <QString>

// Shared appearance for interactive period maps and report maps.
struct MapDisplaySettings {
  double period = 1., tolerancePercent = .1, ellipseKm = .5, arrowKm = 1., colorMin = 0., colorMax = 90.;
  int convention = 0, colorBy = 0, mapLayers = 0;
  QString colormap = "Viridis";
  bool names = false, real = true, imaginary = false, reverse = false;
};
#endif
