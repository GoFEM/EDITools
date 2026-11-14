#ifndef EDITOOLS_TOPOGRAPHY_H
#define EDITOOLS_TOPOGRAPHY_H

#include "core/SurveyCoordinates.h"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <functional>
#include <memory>
#include <limits>

namespace Topography {
// Raster coordinates are absolute WGS84/UTM metres; heights are positive-up metres.
// NaNs represent missing coverage. The immutable grid can be shared by survey copies.
struct Grid {
  int width = 0, height = 0, epsg = 0;
  double west = 0., north = 0., spacing = 0.;
  std::vector<float> values;
  std::string source, datum, notice;
  void validate() const;
  double sample(double easting, double northing) const;
  std::array<double, 2> range() const;
  template <class Archive> void serialize(Archive &ar, unsigned) {
    ar &width &height &epsg &west &north &spacing &values &source &datum &notice;
    if (Archive::is_loading::value)
      validate();
  }
};
struct Style {
  bool visible = true, colors = true, shade = true, contours = false;
  double lower = 0., upper = 1000., opacity = .65, interval = 100.;
  std::string palette = "Terrain";
  template <class Archive> void serialize(Archive &ar, unsigned) {
    ar &visible &colors &shade &contours &lower &upper &opacity &interval &palette;
  }
};
struct State {
  std::shared_ptr<Grid> grid;
  Style style;
  template <class Archive> void serialize(Archive &ar, unsigned) { ar &grid &style; }
};
struct Request {
  int epsg = 0;
  double west = 0., south = 0., east = 0., north = 0., spacing = 30.;
  // Zero reads band units (unknown units must be explicitly supplied).
  double unitScale = 0.;
  std::string source = "User-provided DEM", datum = "Unknown";
  std::string notice;
  void validate() const;
};
using Progress = std::function<bool(double)>;
int utmEpsg(const SurveyCoordinates &coordinates);
Request surveyArea(const std::vector<std::array<double, 3>> &locations, const SurveyCoordinates &coordinates,
                   double marginMetres, double spacing);
// Transform conventional x/y coordinates (longitude/latitude for EPSG:4326).
void transform(int sourceEpsg, int targetEpsg, std::vector<double> &x, std::vector<double> &y);
std::array<double, 4> geographicBounds(const Request &request);
std::shared_ptr<Grid> load(const std::vector<std::string> &files, const Request &request,
                           const Progress &progress = {});
std::shared_ptr<Grid> reproject(const Grid &grid, const Request &request, const Progress &progress = {});
void writeGeoTiff(const std::string &path, const Grid &grid, const SurveyCoordinates &coordinates = {});
void writeXYZ(const std::string &path, const Grid &grid, const SurveyCoordinates &coordinates,
              bool positiveDown);
} // namespace Topography
#endif
