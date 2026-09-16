#ifndef NATIVE_MT_H
#define NATIVE_MT_H

#include "MTResponseData.h"
#include <array>
#include <istream>
#include <map>
#include <string>
#include <vector>

// MT_DATA_EXPORT_SPEC.md. All values here are SI, with phases in degrees.
namespace NativeMT {
struct Mapping {
  RealDataType type;
  const char *observable;
  const char *component;
};
const std::vector<Mapping> &mappings();
const Mapping &mapping(RealDataType type);
using Observation = MTResponseData::Scalar;
using Receivers = std::map<std::string, std::array<double, 3>>;
struct Options {
  std::vector<RealDataType> types;
  std::vector<double> periods;
  Receivers receivers; // Already transformed Cartesian x=north, y=east, z=down.
  std::string coordinate_description;
  bool horizontal_axes_confirmed = false;
  bool conventions_confirmed = false;
  bool negative_time_convention = false;
  bool vertical_magnetic_up = false;
  bool distortion = false;
};
struct Export {
  std::vector<Observation> observations;
  Receivers receivers;
  std::vector<double> frequencies;
  std::vector<std::string> warnings;
  std::string data_text, receiver_text, frequency_text;
  std::size_t masked = 0, missing = 0, undefined_phase = 0;
};

bool same_frequency(double a, double b);
// Parsers preserve sparse scalar keys; no adjacency or rectangular-grid assumptions.
std::vector<Observation> read_observations(std::istream &input);
// Build plot-ready stations, deriving resistivity, phase and phase tensor from
// available impedances. Explicit scalar responses take precedence; gaps stay NaN.
std::map<std::string, MTStationData> read_responses(std::istream &input);
Receivers read_receivers(std::istream &input);
// Four columns: name easting northing elevation. Origin uses that same ordering.
Receivers read_projected_receivers(std::istream &input,
                                  const std::array<double, 3> &origin);
std::array<double, 3> model_coordinates(const std::array<double, 3> &projected,
                                      const std::array<double, 3> &origin);
Export prepare(const std::vector<const MTStationData *> &stations, const Options &options);
// Validate, then serialize. Useful for predicted-data round trips as well.
Export encode(const std::vector<Observation> &observations, const Receivers &receivers,
              const std::string &coordinate_description, bool distortion = false);
}
#endif
