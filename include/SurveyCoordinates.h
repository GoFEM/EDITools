#ifndef SURVEY_COORDINATES_H
#define SURVEY_COORDINATES_H

#include <array>
#include <string>
#include <vector>

// Original station positions remain latitude, longitude, elevation. This
// survey-wide transform supplies map/export coordinates without losing them.
struct SurveyCoordinates
{
  bool utm = false;
  int zone = 31;
  bool north = true;
  bool centered = false;
  double origin_easting = 0, origin_northing = 0;

  static SurveyCoordinates suggested(const std::vector<std::array<double, 3>> &geographic);
  static SurveyCoordinates calculate(const std::vector<std::array<double, 3>> &geographic,
                                     int zone, bool north, bool centered);
  std::vector<std::array<double, 3>> transform(const std::vector<std::array<double, 3>> &geographic) const;
  std::string description() const;

  template<class Archive> void serialize(Archive &ar, const unsigned int)
  {
    ar & utm & zone & north & centered & origin_easting & origin_northing;
  }
};
#endif
