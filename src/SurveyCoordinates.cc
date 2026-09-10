#include "include/SurveyCoordinates.h"
#include <proj.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
void validate(const std::array<double, 3> &p)
{
  if(!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ||
     p[0] < -80 || p[0] > 84 || p[1] < -180 || p[1] > 180)
    throw std::runtime_error("UTM conversion needs finite latitude (-80 to 84 degrees), longitude (-180 to 180 degrees) and elevation.");
}
}

SurveyCoordinates SurveyCoordinates::suggested(const std::vector<std::array<double, 3>> &points)
{
  if(points.empty()) throw std::runtime_error("Load a survey before converting coordinates.");
  // Circular longitude mean also handles surveys straddling the date line.
  double latitude = 0, sine = 0, cosine = 0;
  const double radians = std::acos(-1.) / 180.;
  for(const auto &p: points) {
    validate(p);
    latitude += p[0] / points.size();
    sine += std::sin(p[1] * radians); cosine += std::cos(p[1] * radians);
  }
  const double longitude = std::atan2(sine, cosine) / radians;
  SurveyCoordinates result;
  result.utm = true;
  result.zone = std::max(1, std::min(60, int(std::floor((longitude + 180.) / 6.)) + 1));
  result.north = latitude >= 0.;
  return result;
}

std::vector<std::array<double, 3>> SurveyCoordinates::transform(
    const std::vector<std::array<double, 3>> &points) const
{
  if(!utm) return points;
  if(zone < 1 || zone > 60 || !std::isfinite(origin_easting) || !std::isfinite(origin_northing))
    throw std::runtime_error("Invalid UTM zone or origin.");
  // An explicit WGS84 projection does not need a runtime EPSG database or network.
  const std::string definition = "+proj=utm +zone=" + std::to_string(zone) + " +ellps=WGS84" + (north ? "" : " +south");
  std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)> context(proj_context_create(), proj_context_destroy);
  std::unique_ptr<PJ, decltype(&proj_destroy)> projection(proj_create(context.get(), definition.c_str()), proj_destroy);
  if(!projection) throw std::runtime_error("Could not initialize WGS84 / UTM projection.");
  const double radians = std::acos(-1.) / 180.;
  std::vector<std::array<double, 3>> result;
  result.reserve(points.size());
  for(const auto &p: points) {
    validate(p);
    proj_errno_reset(projection.get());
    const auto converted = proj_trans(projection.get(), PJ_FWD, proj_coord(p[1]*radians, p[0]*radians, p[2], 0));
    if(proj_errno(projection.get()) || !std::isfinite(converted.xy.x) || !std::isfinite(converted.xy.y))
      throw std::runtime_error("A station cannot be projected into the selected UTM zone.");
    // Existing maps and receiver writers expect North/East order. Elevation
    // is preserved; the native MT exporter changes elevation to positive-down z.
    result.push_back({{converted.xy.y-origin_northing, converted.xy.x-origin_easting, p[2]}});
  }
  return result;
}

SurveyCoordinates SurveyCoordinates::calculate(const std::vector<std::array<double, 3>> &points,
                                               int selected_zone, bool selected_north, bool center)
{
  if(points.empty()) throw std::runtime_error("Load a survey before converting coordinates.");
  SurveyCoordinates result;
  result.utm = true; result.zone = selected_zone; result.north = selected_north; result.centered = center;
  const auto projected = result.transform(points);
  if(center) {
    double min_e = projected.front()[1], max_e = min_e;
    double min_n = projected.front()[0], max_n = min_n;
    for(const auto &p: projected) {
      min_e = std::min(min_e, p[1]); max_e = std::max(max_e, p[1]);
      min_n = std::min(min_n, p[0]); max_n = std::max(max_n, p[0]);
    }
    result.origin_easting = min_e + (max_e-min_e)*.5;
    result.origin_northing = min_n + (max_n-min_n)*.5;
  }
  return result;
}

std::string SurveyCoordinates::description() const
{
  if(!utm) return "Latitude / longitude (degrees)";
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(17) << "WGS84 / UTM " << zone << (north ? 'N' : 'S')
       << " (EPSG:" << (north ? 32600 : 32700) + zone << ")\n"
       << "UTM origin (easting northing, m): " << origin_easting << ' ' << origin_northing
       << "\nHorizontal coordinates (m): x = northing - origin_northing; y = easting - origin_easting";
  return text.str();
}
