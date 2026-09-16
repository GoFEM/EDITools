#ifndef MT_MAP_DATA_H
#define MT_MAP_DATA_H

#include "MTStationData.h"

namespace MTMapData {
struct PhaseTensor {
  double phiMin, phiMax; // Principal phases in degrees (atan of tensor axes).
  double skew, azimuth; // Degrees; azimuth clockwise from geographic north.
  double axisRatio;     // |Phi_min / Phi_max|, before taking atan.
};

// No interpolation or extrapolation: tolerance is a relative period difference.
int nearest_period(const dvector &frequencies, double period, double tolerance);
bool phase_tensor(const std::array<double, 4> &tensor, PhaseTensor &result);
bool phase_tensor(const MTStationData &station, unsigned frequency, PhaseTensor &result);
bool induction_vector(const MTStationData &station, unsigned frequency, bool imaginary,
                      std::array<double, 2> &northEast);
}
#endif
