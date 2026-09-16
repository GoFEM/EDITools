#ifndef MT_RESPONSE_DATA_H
#define MT_RESPONSE_DATA_H

#include "MTStationData.h"
#include <istream>

// Shared scalar representation for loaded MT responses, independent of file format.
namespace MTResponseData {
struct Scalar {
  double frequency;
  std::string receiver;
  RealDataType type;
  double value, error;
  template<class Archive> void serialize(Archive &ar, const unsigned int)
  {
    ar & frequency & receiver & type & value & error;
  }
};
enum class Format { Auto, Native, GoFEM };

// Only these readers depend on the input format. Auto inspects the file contents.
std::vector<Scalar> read(std::istream &input, Format format = Format::Auto);
// Build plot-ready stations from validated rows. Missing scalars stay NaN;
// explicit resistivity, phase and tensor values override derived values.
std::map<std::string, MTStationData> stations(const std::vector<Scalar> &rows);
}
#endif
