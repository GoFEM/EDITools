/*
 * The EDI Tools application.
 *
 * Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "include/MTStationData.h"

#include "include/datum.h"

#include <iomanip>
#include <fstream>
#include <algorithm>
#include <math.h>

#include <Eigen/Dense>

using dcomplex = std::complex<double>;

const double mu0 = 4.*M_PI*1e-7;

const std::map<RealDataType, unsigned> type_to_column_table =
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

std::string MTStationData::name() const
{
  return station_name;
}

void MTStationData::set_name(std::string new_name)
{
  station_name = new_name;
}

std::string MTStationData::file_path() const
{
  return file_name;
}

const std::array<double, 3> &MTStationData::position() const
{
  return location;
}

bool MTStationData::active() const
{
  return is_active;
}

void MTStationData::set_active(bool flag)
{
  is_active = flag;
}

const dvector &MTStationData::frequencies() const
{
  return freqs;
}

void MTStationData::set_frequencies(const std::set<double> &fvalues)
{
  if(freqs.size() != fvalues.size())
    throw std::runtime_error("Size mismatch.");

  freqs = std::vector<double>(fvalues.begin(), fvalues.end());
}

bool MTStationData::scalar_value(RealDataType type, unsigned f, double &value, double &error) const
{
  if(!is_active || f >= freqs.size()) return false;
  const auto column = type_to_column_table.find(type);
  if(column == type_to_column_table.end()) return false;
  const auto c = column->second;
  switch(type) {
  case RealZxx: case RealZxy: case RealZyx: case RealZyy:
  case ImagZxx: case ImagZxy: case ImagZyx: case ImagZyy:
    if(!Z_mask[c][f]) return false;
    value = (type == RealZxx || type == RealZxy || type == RealZyx || type == RealZyy) ?
      Z[c][f].real() : Z[c][f].imag();
    error = Z_err_floor[c][f];
    break;
  case RealTzx: case RealTzy: case ImagTzx: case ImagTzy:
    if(!T_mask[c][f]) return false;
    value = (type == RealTzx || type == RealTzy) ? T[c][f].real() : T[c][f].imag();
    error = T_err_floor[c][f];
    break;
  case PTxx: case PTxy: case PTyx: case PTyy:
    if(!PT_mask[c][f]) return false;
    value = PT[c][f]; error = PT_err[c][f];
    break;
  case RhoZxx: case RhoZxy: case RhoZyx: case RhoZyy:
    if(!Z_mask[c][f]) return false;
    value = Rho[c][f]; error = Rho_err[c][f];
    break;
  case PhsZxx: case PhsZxy: case PhsZyx: case PhsZyy:
    if(!Z_mask[c][f]) return false;
    value = Phs[c][f]; error = Phs_err[c][f];
    break;
  default: return false;
  }
  return std::isfinite(value);
}

void MTStationData::set_data(double frequency,
                             const std::vector<RealDataType> &types,
                             const std::vector<double> &values,
                             const std::vector<double> &errors)
{
  if(types.size() != values.size() || types.size() != errors.size())
    throw std::invalid_argument("Data types, values and errors must have the same length.");

  // Prefer the exact key so distinct, closely spaced response frequencies do
  // not overwrite one another. Keep the legacy tolerance for other callers.
  auto it = std::find(freqs.begin(), freqs.end(), frequency);
  if(it == freqs.end())
    it = std::find_if(freqs.begin(), freqs.end(), [frequency](double f) {
      return std::abs(f - frequency) < 1e-10;
    });
  if(it == freqs.end())
    return;
  const auto fidx = std::distance(freqs.begin(), it);

  for(unsigned i = 0; i < types.size(); ++i)
  {
    const auto type = types[i];
    const auto column = type_to_column_table.find(type);
    if(column == type_to_column_table.end()) continue;
    const unsigned c = column->second;
    switch(type) {
    case RealZxx: case RealZxy: case RealZyx: case RealZyy:
      Z[c][fidx].real(values[i]); Z_err[c][fidx] = errors[i]; break;
    case ImagZxx: case ImagZxy: case ImagZyx: case ImagZyy:
      Z[c][fidx].imag(values[i]); Z_err[c][fidx] = errors[i]; break;
    case RealTzx: case RealTzy:
      T[c][fidx].real(values[i]); T_err[c][fidx] = errors[i]; break;
    case ImagTzx: case ImagTzy:
      T[c][fidx].imag(values[i]); T_err[c][fidx] = errors[i]; break;
    case PTxx: case PTxy: case PTyx: case PTyy:
      PT[c][fidx] = values[i]; PT_err[c][fidx] = errors[i]; break;
    case RhoZxx: case RhoZxy: case RhoZyx: case RhoZyy:
      Rho[c][fidx] = values[i]; Rho_err[c][fidx] = errors[i]; break;
    case PhsZxx: case PhsZxy: case PhsZyx: case PhsZyy:
      Phs[c][fidx] = values[i]; Phs_err[c][fidx] = errors[i]; break;
    default: break;
    }
  }

  T_err_floor = T_err;
  Z_err_floor = Z_err;
}

std::vector<std::vector<bool>> MTStationData::impedance_mask() const
{
  if(is_active)
    return Z_mask;
  else
  {
    auto mask = Z_mask;
    for(auto &m: mask)
      std::fill(m.begin(), m.end(), false);

    return mask;
  }
}

std::vector<std::vector<bool> > MTStationData::tipper_mask() const
{
  std::vector<std::vector<bool>> mask(4);
  mask[0] = T_mask[0];
  mask[1] = T_mask[1];
  mask[2] = T_mask[0];
  mask[3] = T_mask[1];

  if(!is_active)
  {
    for(auto &m: mask)
      std::fill(m.begin(), m.end(), false);
  }

  return mask;
}

std::vector<std::vector<bool> > MTStationData::phase_tensor_mask() const
{
  if(is_active)
    return PT_mask;
  else
  {
    auto mask = PT_mask;
    for(auto &m: mask)
      std::fill(m.begin(), m.end(), false);

    return mask;
  }
}

void MTStationData::set_size(const unsigned n_frequencies, bool missing_values)
{
  freqs.assign(n_frequencies, 0.);

  const double initial = missing_values ? std::numeric_limits<double>::quiet_NaN() : 0.;
  Z.assign(4, cvector(n_frequencies, dcomplex(initial, initial)));
  T.assign(2, cvector(n_frequencies, dcomplex(initial, initial)));
  Rho.assign(4, dvector(n_frequencies, initial));
  Phs.assign(4, dvector(n_frequencies, initial));
  PT.assign(4, dvector(n_frequencies, initial));

  Z_err.assign(4, dvector(n_frequencies));
  T_err.assign(2, dvector(n_frequencies));
  Rho_err.assign(4, dvector(n_frequencies));
  Phs_err.assign(4, dvector(n_frequencies));
  PT_err.assign(4, dvector(n_frequencies));
  T_err_floor = T_err;
  Z_err_floor = Z_err;

  Z_mask.assign(4, std::vector<bool>(n_frequencies, true));
  T_mask.assign(2, std::vector<bool>(n_frequencies, true));
  PT_mask.assign(4, std::vector<bool>(n_frequencies, true));

  is_active = true;
}

void MTStationData::decimate()
{
  freqs = decimate_vector(freqs);

  for(size_t i = 0; i < T.size(); ++i)
  {
    T[i] = decimate_vector(T[i]);
    T_err[i] = decimate_vector(T_err[i]);
    T_err_floor[i] = decimate_vector(T_err_floor[i]);
    T_mask[i] = decimate_vector(T_mask[i]);
  }

  for(size_t i = 0; i < Z.size(); ++i)
  {
    Z[i] = decimate_vector(Z[i]);
    Z_err[i] = decimate_vector(Z_err[i]);
    Z_err_floor[i] = decimate_vector(Z_err_floor[i]);
    Z_mask[i] = decimate_vector(Z_mask[i]);

    Rho[i] = decimate_vector(Rho[i]);
    Phs[i] = decimate_vector(Phs[i]);
    Rho_err[i] = decimate_vector(Rho_err[i]);
    Phs_err[i] = decimate_vector(Phs_err[i]);

    PT[i] = decimate_vector(PT[i]);
    PT_err[i] = decimate_vector(PT_err[i]);
    PT_mask[i] = decimate_vector(PT_mask[i]);
  }
}

void MTStationData::get_apparent_resistivity(std::vector<dvector> &data,
                                             std::vector<dvector> &error) const
{
  data = Rho;
  error = Rho_err;
}

void MTStationData::get_phase(std::vector<dvector> &data,
                              std::vector<dvector> &error) const
{
  data = Phs;
  error = Phs_err;
}

void MTStationData::get_tipper(std::vector<dvector> &data,
                               std::vector<dvector> &error) const
{
  std::vector<dvector> T_real(4, dvector(freqs.size()));

  for(unsigned i = 0; i < freqs.size(); ++i)
  {
    T_real[0][i] = T[0][i].real();
    T_real[1][i] = T[1][i].real();
    T_real[2][i] = T[0][i].imag();
    T_real[3][i] = T[1][i].imag();
  }

  std::vector<dvector> T_err_real(4);
  T_err_real[0] = T_err_floor[0];
  T_err_real[1] = T_err_floor[1];
  T_err_real[2] = T_err_floor[0];
  T_err_real[3] = T_err_floor[1];

  data = T_real;
  error = T_err_real;
}

void MTStationData::get_phase_tensor(std::vector<dvector> &data,
                                     std::vector<dvector> &error) const
{
  data = PT;
  error = PT_err;
}

void MTStationData::calculate_apparent_resistivity()
{
  for(unsigned i = 0; i < Rho.size(); ++i)
    for(unsigned fidx = 0; fidx < freqs.size(); ++fidx)
      Rho[i][fidx] = 1. / (mu0 * 2. * M_PI * freqs[fidx]) * pow(std::abs(Z[i][fidx]), 2.);
}

void MTStationData::calculate_phase()
{
  for(unsigned i = 0; i < Rho.size(); ++i)
    for(unsigned fidx = 0; fidx < freqs.size(); ++fidx)
      Phs[i][fidx] = 180. / M_PI * atan2(Z[i][fidx].imag(), Z[i][fidx].real());
}

void MTStationData::calculate_phase_tensor()
{
  for(unsigned fidx = 0; fidx < freqs.size(); ++fidx)
  {
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> X, Y;

    X(0,0) = Z[0][fidx].real(); X(0,1) = Z[1][fidx].real();
    X(1,0) = Z[2][fidx].real(); X(1,1) = Z[3][fidx].real();

    Y(0,0) = Z[0][fidx].imag(); Y(0,1) = Z[1][fidx].imag();
    Y(1,0) = Z[2][fidx].imag(); Y(1,1) = Z[3][fidx].imag();

    auto PT_mat = X.inverse() * Y;

    PT[0][fidx] = PT_mat(0,0);
    PT[1][fidx] = PT_mat(0,1);
    PT[2][fidx] = PT_mat(1,0);
    PT[3][fidx] = PT_mat(1,1);
  }
}

void MTStationData::set_error_floor(double val)
{
  error_floor = val;
  apply_error_floor();
  propagate_rho_phase_error();
  propagate_phase_tensor_error();
}

void MTStationData::apply_error_floor()
{
  for(size_t i = 0; i < T.size(); ++i)
  {
    for(size_t j = 0; j < T[i].size(); ++j)
    {
      const double err = error_floor;
      //T_err_floor[i][j] = (err > T_err[i][j]) ? err : T_err[i][j];
      T_err_floor[i][j] = err;
    }
  }

  for(size_t j = 0; j < Z[0].size(); ++j)
  {
    {
      const double err = std::abs(Z[1][j]) * error_floor;
      Z_err_floor[0][j] = (err > Z_err[0][j]) ? err : Z_err[0][j];
      Z_err_floor[1][j] = (err > Z_err[1][j]) ? err : Z_err[1][j];
    }

    {
      const double err = std::abs(Z[2][j]) * error_floor;
      Z_err_floor[2][j] = (err > Z_err[2][j]) ? err : Z_err[2][j];
      Z_err_floor[3][j] = (err > Z_err[3][j]) ? err : Z_err[3][j];
    }
  }
}

const std::vector<RealDataType> &MTStationData::linked_mask_types(RealDataType type)
{
  static const std::vector<RealDataType> none;
  static const std::vector<RealDataType> impedance{RealZxx, RealZxy, RealZyx, RealZyy};
  static const std::vector<RealDataType> tensor{PTxx, PTxy, PTyx, PTyy};
  switch(type) {
  case PTxx: case PTxy: case PTyx: case PTyy: return impedance;
  case RealZxx: case RealZxy: case RealZyx: case RealZyy:
  case ImagZxx: case ImagZxy: case ImagZyx: case ImagZyy:
  case RhoZxx: case RhoZxy: case RhoZyx: case RhoZyy:
  case PhsZxx: case PhsZxy: case PhsZyx: case PhsZyy: return tensor;
  default: return none;
  }
}

std::vector<bool> &MTStationData::mask_for(RealDataType type)
{
  const auto it = type_to_column_table.find(type);
  if(it == type_to_column_table.end()) throw std::invalid_argument("Unsupported data type.");
  const unsigned column = it->second;
  switch(type) {
  case RealTzx: case RealTzy: case ImagTzx: case ImagTzy: return T_mask[column];
  case PTxx: case PTxy: case PTyx: case PTyy: return PT_mask[column];
  default: return Z_mask[column];
  }
}

void MTStationData::mask_type(RealDataType type, bool on, bool linkTensorMasks)
{
  auto &mask = mask_for(type);
  std::fill(mask.begin(), mask.end(), on);
  if(linkTensorMasks) for(auto related: linked_mask_types(type)) mask_type(related, on, false);
}

void MTStationData::set_data_mask(RealDataType type, double frequency, bool on, bool linkTensorMasks)
{
  auto &mask = mask_for(type);
  auto match = std::find(freqs.begin(), freqs.end(), frequency);
  if(match == freqs.end()) {
    double closest = 1e-3;
    for(auto it = freqs.begin(); it != freqs.end(); ++it) {
      const double distance = std::abs(frequency - *it) / *it;
      if(distance < closest) { closest = distance; match = it; }
    }
  }
  if(match == freqs.end()) return;
  const auto index = std::distance(freqs.begin(), match);
  mask[index] = on;
  if(linkTensorMasks) for(auto related: linked_mask_types(type)) set_data_mask(related, freqs[index], on, false);
}

std::set<RealDataType> MTStationData::active_types() const
{
  std::set<RealDataType> types;

  std::vector<RealDataType> types_v = {RealZxx, RealZxy, RealZyx, RealZyy};
  for(unsigned i = 0; i < Z_mask.size(); ++i)
  {
    bool flag = false;
    for(unsigned j = 0; j < Z_mask[i].size(); ++j)
    {
      flag |= Z_mask[i][j];
    }

    if(flag)
      types.insert(types_v[i]);
  }

  types_v = {RealTzx, RealTzy};
  for(unsigned i = 0; i < T_mask.size(); ++i)
  {
    bool flag = false;
    for(unsigned j = 0; j < T_mask[i].size(); ++j)
    {
      flag |= T_mask[i][j];
    }

    if(flag)
      types.insert(types_v[i]);
  }

  types_v = {PTxx, PTxy, PTyx, PTyy};
  for(unsigned i = 0; i < PT_mask.size(); ++i)
  {
    bool flag = false;
    for(unsigned j = 0; j < PT_mask[i].size(); ++j)
    {
      flag |= PT_mask[i][j];
    }

    if(flag)
      types.insert(types_v[i]);
  }

  return types;
}

void MTStationData::propagate_rho_phase_error()
{
  for(size_t i = 0; i < Z.size(); ++i)
  {
    for(size_t j = 0; j < Z[i].size(); ++j)
    {
      const dcomplex z = Z[i][j];
      const double ze = Z_err_floor[i][j];
      const double omega = 2. * M_PI * freqs[j];
      const double absz = std::abs(z);
      Rho_err[i][j] = std::abs(sqrt(pow((2*absz)/(mu0*omega), 2.) * ze*ze));
      Phs_err[i][j] = 180. / M_PI * ze / absz;
    }
  }

//  for(size_t j = 0; j < Z[0].size(); ++j)
//  {
//    std::cerr << freqs[j] << "\t";
//    for(size_t i = 0; i < Z.size(); ++i)
//    {
//      std::cerr << Rho[i][j] << "\t" << Rho_err[i][j] << "\t";
//    }
//    for(size_t i = 0; i < Z.size(); ++i)
//    {
//      std::cerr << Phs[i][j] << "\t" << Phs_err[i][j] << "\t";
//    }
//    std::cerr << "\n";
//  }
}

void MTStationData::propagate_phase_tensor_error()
{
  double dPTdX[2][2][2][2],
         dPTdY[2][2][2][2];

  dcomplex Zt[2][2];
  double PTt[2][2], PTt_err[2][2];

  for(size_t f = 0; f < freqs.size(); ++f)
  {
    PTt[0][0] = PT[0][f];
    PTt[0][1] = PT[1][f];
    PTt[1][0] = PT[2][f];
    PTt[1][1] = PT[3][f];

    Zt[0][0] = Z[0][f];
    Zt[0][1] = Z[1][f];
    Zt[1][0] = Z[2][f];
    Zt[1][1] = Z[3][f];

    const double detX = Zt[0][0].real()*Zt[1][1].real() - Zt[1][0].real()*Zt[0][1].real();

    {
      dPTdX[0][0][0][0] =(-PTt[0][0] * Zt[1][1].real()) / detX;
      dPTdX[0][0][0][1] =( PTt[0][0] * Zt[1][0].real() - Zt[1][0].imag()) / detX;
      dPTdX[0][0][1][0] =( PTt[0][0] * Zt[0][1].real()) / detX;
      dPTdX[0][0][1][1] =(-PTt[0][0] * Zt[0][0].real() + Zt[0][0].imag()) / detX;

      dPTdY[0][0][0][0] = Zt[1][1].real() / detX;
      dPTdY[0][0][0][1] = 0;
      dPTdY[0][0][1][0] =-Zt[0][1].real() / detX;
      dPTdY[0][0][1][1] = 0;
    }

    // dPTxy
    {
      dPTdX[0][1][0][0] = (-PTt[0][1] * Zt[1][1].real()) / detX;
      dPTdX[0][1][0][1] = ( PTt[0][1] * Zt[1][0].real() - Zt[1][1].imag()) / detX;
      dPTdX[0][1][1][0] = ( PTt[0][1] * Zt[0][1].real()) / detX;
      dPTdX[0][1][1][1] = (-PTt[0][1] * Zt[0][0].real() + Zt[0][1].imag()) / detX;

      dPTdY[0][1][0][0] = 0;
      dPTdY[0][1][0][1] = Zt[1][1].real() / detX;
      dPTdY[0][1][1][0] = 0;
      dPTdY[0][1][1][1] =-Zt[0][1].real() / detX;
    }

    // dPTyx
    {
      dPTdX[1][0][0][0] = (-PTt[1][0] * Zt[1][1].real() + Zt[1][0].imag()) / detX;
      dPTdX[1][0][0][1] = ( PTt[1][0] * Zt[1][0].real()) / detX;
      dPTdX[1][0][1][0] = ( PTt[1][0] * Zt[0][1].real() - Zt[0][0].imag()) / detX;
      dPTdX[1][0][1][1] = (-PTt[1][0] * Zt[0][0].real()) / detX;

      dPTdY[1][0][0][0] =-Zt[1][0].real() / detX;
      dPTdY[1][0][0][1] = 0;
      dPTdY[1][0][1][0] = Zt[0][0].real() / detX;
      dPTdY[1][0][1][1] = 0;
    }

    // dPTyy
    {
      dPTdX[1][1][0][0] = (-PTt[1][1] * Zt[1][1].real() + Zt[1][1].imag()) /  detX;
      dPTdX[1][1][0][1] = ( PTt[1][1] * Zt[1][0].real()) / detX;
      dPTdX[1][1][1][0] = ( PTt[1][1] * Zt[0][1].real() - Zt[0][1].imag()) / detX;
      dPTdX[1][1][1][1] = (-PTt[1][1] * Zt[0][0].real()) / detX;

      dPTdY[1][1][0][0] = 0;
      dPTdY[1][1][0][1] =-Zt[1][0].real() / detX;
      dPTdY[1][1][1][0] = 0;
      dPTdY[1][1][1][1] = Zt[0][0].real() / detX;
    }

    for(int k = 0; k < 2; ++k)
      for(int l = 0; l < 2; ++l)
      {
        double propagated_error = 0;
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            propagated_error += pow(dPTdY[k][l][i][j], 2.)*pow(Z_err_floor[j+2*i][f], 2.) +
                                pow(dPTdX[k][l][i][j], 2.)*pow(Z_err_floor[j+2*i][f], 2.);

        PTt_err[k][l] = sqrt(propagated_error);
      }

    PT_err[0][f] = PTt_err[0][0];
    PT_err[1][f] = PTt_err[0][1];
    PT_err[2][f] = PTt_err[1][0];
    PT_err[3][f] = PTt_err[1][1];

//    std::cerr << freqs[f] << "\t";
//    for(size_t i = 0; i < Z.size(); ++i)
//    {
//      std::cerr << PT[i][f] << "\t" << PT_err[i][f] << "\t";
//    }
//    std::cerr << "\n";
  }
}

void MTStationData::write(std::ofstream &ofs, const std::vector<RealDataType> &types,
                          const std::vector<double> &selected_periods,
                          std::set<double> *written_frequencies) const
{
  ofs << std::setprecision(std::numeric_limits<double>::max_digits10);

  for(unsigned f = 0; f < freqs.size(); ++f)
  {
    const double period = 1. / freqs[f];

    bool is_selected = false;
    for(double selected_period: selected_periods)
    {
      if(fabs(selected_period - period) / period < 1e-5)
      {
        is_selected = true;
        break;
      }
    }

    if(!is_selected)
      continue;

    std::set<RealDataType> processed_types;
    const auto start = ofs.tellp();

    for(RealDataType type: types)
    {
      if(processed_types.count(type))
        continue;

      switch (type) {
      case RealZxx:
      case ImagZxx:
        if(Z_mask[0][f] && !std::isnan(std::abs(Z[0][f])))
        {
          ofs << Datum::convert_type_to_string(RealZxx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[0][f].real() << "\t" << Z_err_floor[0][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagZxx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[0][f].imag() << "\t" << Z_err_floor[0][f] << "\n";
        }
        break;
      case RealZxy:
      case ImagZxy:
        if(Z_mask[1][f] && !std::isnan(std::abs(Z[1][f])))
        {
          ofs << Datum::convert_type_to_string(RealZxy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[1][f].real() << "\t" << Z_err_floor[1][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagZxy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[1][f].imag() << "\t" << Z_err_floor[1][f] << "\n";
        }
        break;
      case RealZyx:
      case ImagZyx:
        if(Z_mask[2][f] && !std::isnan(std::abs(Z[2][f])))
        {
          ofs << Datum::convert_type_to_string(RealZyx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[2][f].real() << "\t" << Z_err_floor[2][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagZyx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[2][f].imag() << "\t" << Z_err_floor[2][f] << "\n";
        }
        break;
      case RealZyy:
      case ImagZyy:
        if(Z_mask[3][f] && !std::isnan(std::abs(Z[3][f])))
        {
          ofs << Datum::convert_type_to_string(RealZyy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[3][f].real() << "\t" << Z_err_floor[3][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagZyy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Z[3][f].imag() << "\t" << Z_err_floor[3][f] << "\n";
        }
        break;
      case RealTzx:
      case ImagTzx:
        if(T_mask[0][f] && std::abs(T[0][f]) > 0.
              && !std::isnan(std::abs(T[0][f])))
        {
          ofs << Datum::convert_type_to_string(RealTzx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << T[0][f].real() << "\t" << T_err_floor[0][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagTzx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << T[0][f].imag() << "\t" << T_err_floor[0][f] << "\n";
        }
        break;
      case RealTzy:
      case ImagTzy:
        if(T_mask[1][f] && std::abs(T[1][f]) > 0.
             && !std::isnan(std::abs(T[1][f])))
        {
          ofs << Datum::convert_type_to_string(RealTzy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << T[1][f].real() << "\t" << T_err_floor[1][f] << "\n";
          ofs << Datum::convert_type_to_string(ImagTzy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << T[1][f].imag() << "\t" << T_err_floor[1][f] << "\n";
        }
        break;
      case PTxx:
        if(PT_mask[0][f] && !std::isnan(std::abs(PT[0][f])))
        {
          ofs << Datum::convert_type_to_string(PTxx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << PT[0][f] << "\t" << PT_err[0][f] << "\n";
        }
        break;
      case PTxy:
        if(PT_mask[1][f] && !std::isnan(std::abs(PT[1][f])))
        {
          ofs << Datum::convert_type_to_string(PTxy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << PT[1][f] << "\t" << PT_err[1][f] << "\n";
        }
        break;
      case PTyx:
        if(PT_mask[2][f] && !std::isnan(std::abs(PT[2][f])))
        {
          ofs << Datum::convert_type_to_string(PTyx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << PT[2][f] << "\t" << PT_err[2][f] << "\n";
        }
        break;
      case PTyy:
        if(PT_mask[3][f] && !std::isnan(std::abs(PT[3][f])))
        {
          ofs << Datum::convert_type_to_string(PTyy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << PT[3][f] << "\t" << PT_err[3][f] << "\n";
        }
        break;
      case RhoZxx:
        if(Z_mask[0][f] && !std::isnan(std::abs(Rho[0][f])))
        {
          ofs << Datum::convert_type_to_string(RhoZxx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Rho[0][f] << "\t" << Rho_err[0][f] << "\n";
        }
        break;
      case RhoZxy:
        if(Z_mask[1][f] && !std::isnan(std::abs(Rho[1][f])))
        {
          ofs << Datum::convert_type_to_string(RhoZxy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Rho[1][f] << "\t" << Rho_err[1][f] << "\n";
        }
        break;
      case RhoZyx:
        if(Z_mask[2][f] && !std::isnan(std::abs(Rho[2][f])))
        {
          ofs << Datum::convert_type_to_string(RhoZyx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Rho[2][f] << "\t" << Rho_err[2][f] << "\n";
        }
        break;
      case RhoZyy:
        if(Z_mask[3][f] && !std::isnan(std::abs(Rho[3][f])))
        {
          ofs << Datum::convert_type_to_string(RhoZyy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Rho[3][f] << "\t" << Rho_err[3][f] << "\n";
        }
        break;
      case PhsZxx:
        if(Z_mask[0][f] && !std::isnan(std::abs(Phs[3][f])))
        {
          ofs << Datum::convert_type_to_string(PhsZxx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Phs[0][f] << "\t" << Phs_err[0][f] << "\n";
        }
        break;
      case PhsZxy:
        if(Z_mask[1][f] && !std::isnan(std::abs(Phs[1][f])))
        {
          ofs << Datum::convert_type_to_string(PhsZxy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Phs[1][f] << "\t" << Phs_err[1][f] << "\n";
        }
        break;
      case PhsZyx:
        if(Z_mask[2][f] && !std::isnan(std::abs(Phs[2][f])))
        {
          ofs << Datum::convert_type_to_string(PhsZyx) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Phs[2][f] << "\t" << Phs_err[2][f] << "\n";
        }
        break;
      case PhsZyy:
        if(Z_mask[3][f] && !std::isnan(std::abs(Phs[3][f])))
        {
          ofs << Datum::convert_type_to_string(PhsZyy) << "\t"
              << freqs[f] << "\tPlane_wave\t" << station_name << "\t"
              << Phs[3][f] << "\t" << Phs_err[3][f] << "\n";
        }
        break;
      default:
        break;
      }

      processed_types.insert(type);
      RealDataType complement_type = Datum::get_complex_complement(type);
      if(complement_type != InvalidType)
        processed_types.insert(complement_type);
    }
    if(written_frequencies && ofs.tellp() > start) written_frequencies->insert(freqs[f]);
  }
}

double MTStationData::rms(const MTStationData &other) const
{
  unsigned n = 0;
  double residual = 0;

  auto add = [&](double observed, double predicted, double error,
                 double response_error, bool active, bool phase = false) {
    if(!active || !std::isfinite(observed) || !std::isfinite(predicted) ||
       !std::isfinite(error) || error <= 0 ||
       !std::isfinite(response_error) || response_error <= 0)
      return;
    const double difference = phase ? std::remainder(observed - predicted, 360.) : observed - predicted;
    residual += std::pow(difference / error, 2.);
    ++n;
  };

  if(is_active && other.is_active)
    for(size_t j = 0; j < freqs.size(); ++j) {
      auto match = std::find(other.freqs.begin(), other.freqs.end(), freqs[j]);
      if(match == other.freqs.end()) {
        // EDI and exported response frequencies may differ in printed precision.
        double closest = 1e-3;
        for(auto it = other.freqs.begin(); it != other.freqs.end(); ++it) {
          const double distance = std::abs(*it - freqs[j]) / freqs[j];
          if(distance < closest) { closest = distance; match = it; }
        }
      }
      if(match == other.freqs.end()) continue;
      const auto k = std::distance(other.freqs.begin(), match);
      for(size_t i = 0; i < Z.size(); ++i) {
        const bool impedance_active = Z_mask[i][j] && other.Z_mask[i][k];
        add(Z[i][j].real(), other.Z[i][k].real(), Z_err_floor[i][j], other.Z_err_floor[i][k], impedance_active);
        add(Z[i][j].imag(), other.Z[i][k].imag(), Z_err_floor[i][j], other.Z_err_floor[i][k], impedance_active);
        add(Rho[i][j], other.Rho[i][k], Rho_err[i][j], other.Rho_err[i][k], impedance_active);
        add(Phs[i][j], other.Phs[i][k], Phs_err[i][j], other.Phs_err[i][k], impedance_active, true);
        add(PT[i][j], other.PT[i][k], PT_err[i][j], other.PT_err[i][k], PT_mask[i][j] && other.PT_mask[i][k]);
        if(i < T.size()) {
          const bool tipper_active = T_mask[i][j] && other.T_mask[i][k];
          add(T[i][j].real(), other.T[i][k].real(), T_err_floor[i][j], other.T_err_floor[i][k], tipper_active);
          add(T[i][j].imag(), other.T[i][k].imag(), T_err_floor[i][j], other.T_err_floor[i][k], tipper_active);
        }
      }
    }

  return n ? std::sqrt(residual / n) : std::numeric_limits<double>::quiet_NaN();
}
