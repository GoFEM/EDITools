#pragma once

#include <array>
#include <string>

// Display settings belong to a project; changing them never changes observations.
struct TFViewStyle
{
  // XX, XY, YX, YY, masked tensor, real tipper, imaginary tipper, masked tipper.
  std::array<std::string, 8> colors{{"#579bc9", "#1d5fa2", "#b45309", "#d08738",
                                  "#9ca3af", "#111827", "#6b7280", "#bfc4cc"}};
  double markerSize = 5., lineWidth = 1., arrowScale = 90., ellipseSize = 26.;
  std::string colorMap = "Viridis";
  bool reverseColorMap = false;
  double colorMin = 0., colorMax = 90.;

  template<class Archive> void serialize(Archive &ar, const unsigned int)
  {
    ar & colors & markerSize & lineWidth & arrowScale & ellipseSize;
    ar & colorMap & reverseColorMap & colorMin & colorMax;
  }
};

struct TFViewRange
{
  bool automatic = true;
  double lower = .01, upper = 1000.;

  template<class Archive> void serialize(Archive &ar, const unsigned int)
  { ar & automatic & lower & upper; }
};
