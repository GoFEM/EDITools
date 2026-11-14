#ifndef PERIOD_RESAMPLING_INFO_H
#define PERIOD_RESAMPLING_INFO_H

#include <string>
#include <vector>
#include <boost/serialization/vector.hpp>

namespace PeriodResampling {
enum class Status { Existing, Interpolated, Derived, Masked, Missing, Outside, Gap, InvalidError };
// Components 0..3: Zxx, Zxy, Zyx, Zyy; 4..5: Tzx, Tzy; 6..9: PTxx, PTxy, PTyx, PTyy.
struct Record {
  std::string station;
  unsigned component = 0;
  double period = 0., lower = 0., upper = 0.;
  Status status = Status::Missing;
  template<class Archive> void serialize(Archive &ar, const unsigned int) {
    ar & station & component & period & lower & upper & status;
  }
};
struct Info {
  std::string method;
  double maximumRatio = 2.;
  bool maskedBarriers = true;
  std::vector<double> periods;
  std::vector<Record> records;
  template<class Archive> void serialize(Archive &ar, const unsigned int) {
    ar & method & maximumRatio & maskedBarriers & periods & records;
  }
};
}
#endif
