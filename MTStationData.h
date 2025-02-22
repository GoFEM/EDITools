#ifndef MT_STATION_DATA_H
#define MT_STATION_DATA_H

#include <string>
#include <vector>
#include <array>
#include <complex>
#include <map>
#include <set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "datum.h"

using dvector = std::vector<double>;
using cvector = std::vector<std::complex<double>>;

class weak_compare : public std::binary_function<double,double,bool>
{
public:
  weak_compare( double arg_ = 1e-6 ) : epsilon(arg_) {}
  bool operator()( const double &left, const double &right  ) const
  {
    return (fabs(left - right) > epsilon) && (left < right);
  }
  double epsilon;
};

static std::map<RealDataType, unsigned> type_to_column_table =
{
  {RhoZxx, 0}, {PhsZxx, 0},
  {RhoZxy, 1}, {PhsZxy, 1},
  {RhoZyx, 2}, {PhsZyx, 2},
  {RhoZyy, 3}, {PhsZyy, 3},
  {RealZxx, 0}, {ImagZxx, 0},
  {RealZxy, 1}, {ImagZxy, 1},
  {RealZyx, 2}, {ImagZyx, 2},
  {RealZyy, 3}, {ImagZyy, 3},
  {RealTzy, 1}, {ImagTzy, 1},
  {RealTzx, 0}, {ImagTzx, 0},
  {PTxx, 0}, {PTxy, 1},
  {PTyx, 2}, {PTyy, 3}
};

class MTStationData
{
  friend class boost::serialization::access;
  friend class EDIFileReader;

private:
  std::string station_name;
  std::string file_name;

  std::vector<double> freqs;
  std::array<double, 3> location;

  std::vector<cvector> Z;
  std::vector<dvector> Z_err, Z_err_floor;

  std::vector<cvector> T;
  std::vector<dvector> T_err, T_err_floor;

  std::vector<dvector> Rho, Phs, PT;
  std::vector<dvector> Rho_err, Phs_err, PT_err;

  std::vector<std::vector<bool>> Z_mask, T_mask, PT_mask;

  double error_floor;

  bool is_active;

public:
  MTStationData():
    error_floor(0)
  {}

  std::string name() const;
  void set_name(std::string new_name);

  std::string file_path() const;
  const std::array<double, 3> &position() const;

  bool active() const;
  void set_active(bool flag);

  const dvector &frequencies() const;

  void set_frequencies(const std::set<double> &fvalues);
  void set_data(double frequency,
                const std::vector<RealDataType> &types,
                const std::vector<double> &values,
                const std::vector<double> &errors);

  std::vector<std::vector<bool>> impedance_mask() const;
  std::vector<std::vector<bool>> tipper_mask() const;
  std::vector<std::vector<bool>> phase_tensor_mask() const;

  void set_size(const unsigned n_frequencies);

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

  void mask_type(RealDataType type, bool on);
  void set_data_mask(RealDataType type, double frequency, bool on);

  std::set<RealDataType> active_types() const;

  void write(std::ofstream &ofs,
             const std::vector<RealDataType> &types,
             const std::vector<double> &selected_periods) const;

  // Calculates RMS with other stations
  double rms(const MTStationData &other) const;

private:
  template<class U>
  std::vector<U> decimate_vector(const std::vector<U> &vec);

  std::complex<double> sgn(const std::complex<double> &val) const
  {
    return val / std::abs(val);
  }

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
