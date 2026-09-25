#ifndef MT_STATION_DATA_H
#define MT_STATION_DATA_H

#include <string>
#include <vector>
#include <array>
#include <complex>
#include <map>
#include <set>
#include <limits>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "include/datum.h"

using dvector = std::vector<double>;
using cvector = std::vector<std::complex<double>>;

// Shared, read-only mapping from supported scalar types to tensor/vector columns.
extern const std::map<RealDataType, unsigned> type_to_column_table;


class MTStationData
{
  friend class boost::serialization::access;
  friend class EDIFileReader;
  friend struct NativeMTStationAccess;
  friend struct PeriodResamplingAccess;
  friend struct EDIPeriodMergeAccess;

private:
  std::string station_name;
  std::string file_name;

  std::vector<double> freqs;
  std::array<double, 3> location{{std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN()}};

  std::vector<cvector> Z;
  std::vector<dvector> Z_err, Z_err_floor;

  std::vector<cvector> T;
  std::vector<dvector> T_err, T_err_floor;

  std::vector<dvector> Rho, Phs, PT;
  std::vector<dvector> Rho_err, Phs_err, PT_err;

  std::vector<std::vector<bool>> Z_mask, T_mask, PT_mask;

  double error_floor = 0.;

  bool is_active = true;

public:
  MTStationData() = default;

  std::string name() const;
  void set_name(std::string new_name);

  std::string file_path() const;
  const std::array<double, 3> &position() const;
  void set_position(const std::array<double, 3> &position) { location = position; }

  bool active() const;
  void set_active(bool flag);

  const dvector &frequencies() const;
  // Returns false for masked or missing scalars. Error validation is left to
  // the caller so comparisons can choose observed or inversion errors.
  bool scalar_value(RealDataType type, unsigned frequency_index, double &value, double &error) const;

  void set_frequencies(const std::set<double> &fvalues);
  void set_data(double frequency,
                const std::vector<RealDataType> &types,
                const std::vector<double> &values,
                const std::vector<double> &errors);

  std::vector<std::vector<bool>> impedance_mask() const;
  std::vector<std::vector<bool>> tipper_mask() const;
  std::vector<std::vector<bool>> phase_tensor_mask() const;

  // Initialize all arrays and masks, replacing any previously loaded samples.
  void set_size(const unsigned n_frequencies, bool missing_values = false);

  void decimate();

  void get_apparent_resistivity(std::vector<dvector> &data,
                                std::vector<dvector> &error) const;

  void get_phase(std::vector<dvector> &data,
                 std::vector<dvector> &error) const;

  void get_tipper(std::vector<dvector> &data,
                  std::vector<dvector> &error) const;

  void get_phase_tensor(std::vector<dvector> &data,
                        std::vector<dvector> &error) const;

  // Calculate transfer functions using stored impedance tensor
  void calculate_apparent_resistivity();
  void calculate_phase();
  void calculate_phase_tensor();

  void set_error_floor(double val);

  // Optional UI linkage follows Phi = Re(Z)^-1 Im(Z): a complex Z component
  // relates to every PT entry, and any PT entry to all complex Z components.
  // Apply only to the opposite tensor, without recursive propagation.
  void mask_type(RealDataType type, bool on, bool linkTensorMasks = false);
  void set_data_mask(RealDataType type, double frequency, bool on, bool linkTensorMasks = false);
  static const std::vector<RealDataType> &linked_mask_types(RealDataType type);

  std::set<RealDataType> active_types() const;

  void write(std::ofstream &ofs,
             const std::vector<RealDataType> &types,
             const std::vector<double> &selected_periods,
             std::set<double> *written_frequencies = nullptr) const;

  // Calculates RMS with other stations
  double rms(const MTStationData &other) const;

private:
  template<class U>
  std::vector<U> decimate_vector(const std::vector<U> &vec);


  std::vector<bool> &mask_for(RealDataType type);
  void apply_error_floor();
  void propagate_rho_phase_error();
  void propagate_phase_tensor_error();

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & station_name;
      ar & file_name;
      ar & freqs;
      ar & location;
      ar & Z;
      ar & Z_err;
      ar & Z_err_floor;
      ar & T;
      ar & T_err;
      ar & T_err_floor;
      ar & Rho;
      ar & Phs;
      ar & PT;
      ar & Rho_err;
      ar & Phs_err;
      ar & PT_err;
      ar & T_mask;
      ar & Z_mask;
      ar & PT_mask;
      ar & is_active;
      ar & error_floor;
  }
};

template<class U>
std::vector<U> MTStationData::decimate_vector(const std::vector<U> &vec)
{
  std::vector<U> vec_dec;
  vec_dec.reserve(vec.size() / 2);

  for(size_t i = 0; i < vec.size(); i += 2)
    vec_dec.push_back(vec[i]);

  return vec_dec;
}

#endif // MT_STATION_DATA_H
