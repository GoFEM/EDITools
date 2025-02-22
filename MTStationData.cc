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

#include "MTStationData.h"

#include "datum.h"

#include <iomanip>
#include <fstream>
#include <math.h>

#include <Eigen/Dense>

using dcomplex = std::complex<double>;

const double mu0 = 4.*M_PI*1e-7;

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

void MTStationData::set_data(double frequency,
                             const std::vector<RealDataType> &types,
                             const std::vector<double> &values,
                             const std::vector<double> &errors)
{
  unsigned fidx = std::numeric_limits<unsigned>::max();
  for(unsigned i = 0; i < freqs.size(); ++i)
    if(fabs(freqs[i] - frequency) < 1e-10)
    {
      fidx = i;
      break;
    }

  if(fidx > freqs.size())
    return;

  for(unsigned i = 0; i < types.size(); ++i)
  {
    switch (types[i])
    {
    case RealZxx:
      Z[0][fidx].real(values[i]);
      Z_err[0][fidx] = errors[i];
      break;
    case RealZxy:
      Z[1][fidx].real(values[i]);
      Z_err[1][fidx] = errors[i];
      break;
    case RealZyx:
      Z[2][fidx].real(values[i]);
      Z_err[2][fidx] = errors[i];
      break;
    case RealZyy:
      Z[3][fidx].real(values[i]);
      Z_err[3][fidx] = errors[i];
      break;
    case ImagZxx:
      Z[0][fidx].imag(values[i]);
      Z_err[0][fidx] = errors[i];
      break;
    case ImagZxy:
      Z[1][fidx].imag(values[i]);
      Z_err[1][fidx] = errors[i];
      break;
    case ImagZyx:
      Z[2][fidx].imag(values[i]);
      Z_err[2][fidx] = errors[i];
      break;
    case ImagZyy:
      Z[3][fidx].imag(values[i]);
      Z_err[3][fidx] = errors[i];
      break;
    case PTxx:
      PT[0][fidx] = values[i];
      PT_err[0][fidx] = errors[i];
      break;
    case PTxy:
      PT[1][fidx] = values[i];
      PT_err[1][fidx] = errors[i];
      break;
    case PTyx:
      PT[2][fidx] = values[i];
      PT_err[2][fidx] = errors[i];
      break;
    case PTyy:
      PT[3][fidx] = values[i];
      PT_err[3][fidx] = errors[i];
      break;
    case RealTzx:
      T[0][fidx].real(values[i]);
      T_err[0][fidx] = errors[i];
      break;
    case RealTzy:
      T[1][fidx].real(values[i]);
      T_err[1][fidx] = errors[i];
      break;
    case ImagTzx:
      T[0][fidx].imag(values[i]);
      T_err[0][fidx] = errors[i];
      break;
    case ImagTzy:
      T[1][fidx].imag(values[i]);
      T_err[1][fidx] = errors[i];
      break;
    default:
      break;
    }
  }

  T_err_floor = T_err;
  Z_err_floor = Z_err;
  PT_err = PT_err;
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

void MTStationData::set_size(const unsigned n_frequencies)
{
  freqs.resize(n_frequencies);

  Z.resize(4, cvector(n_frequencies));
  T.resize(2, cvector(n_frequencies));
  Rho.resize(4, dvector(n_frequencies));
  Phs.resize(4, dvector(n_frequencies));
  PT.resize(4, dvector(n_frequencies));

  Z_err.resize(4, dvector(n_frequencies));
  T_err.resize(2, dvector(n_frequencies));
  Rho_err.resize(4, dvector(n_frequencies));
  Phs_err.resize(4, dvector(n_frequencies));
  PT_err.resize(4, dvector(n_frequencies));
  T_err_floor = T_err;
  Z_err_floor = Z_err;

  Z_mask.resize(4, std::vector<bool>(n_frequencies, true));
  T_mask.resize(2, std::vector<bool>(n_frequencies, true));
  PT_mask.resize(4, std::vector<bool>(n_frequencies, true));

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

void MTStationData::mask_type(RealDataType type, bool on)
{
  auto it = type_to_column_table.find(type);
  if(it == type_to_column_table.end())
    throw std::runtime_error("Unsupported type.");

  unsigned column = it->second;

  std::vector<bool>* v;
  switch (type) {
  case RealTzx:
  case RealTzy:
    v = &T_mask[column];
    break;
  case PTxx:
  case PTxy:
  case PTyx:
  case PTyy:
    v = &PT_mask[column];
    break;
  default:
    v = &Z_mask[column];
    break;
  }

  std::fill(v->begin(), v->end(), on);
}

void MTStationData::set_data_mask(RealDataType type, double frequency, bool on)
{
  auto it = type_to_column_table.find(type);
  if(it == type_to_column_table.end())
    throw std::runtime_error("Unsupported type.");

  unsigned column = it->second;

  std::vector<bool>* v;
  switch (type) {
  case RealTzx:
  case RealTzy:
  case ImagTzx:
  case ImagTzy:
    v = &T_mask[column];
    break;
  case PTxx:
  case PTxy:
  case PTyx:
  case PTyy:
    v = &PT_mask[column];
    break;
  default:
    v = &Z_mask[column];
    break;
  }

  for(unsigned i = 0; i < freqs.size(); ++i)
    if(fabs(frequency - freqs[i]) / freqs[i] < 1e-3)
    {
      (*v)[i] = on;
      break;
    }
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
                          const std::vector<double> &selected_periods) const
{
  ofs << std::setprecision(10);

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
  }
}

double MTStationData::rms(const MTStationData &other) const
{
  unsigned n = 0;

  double residual = 0;
  for(size_t i = 0; i < Z.size(); ++i)
  {
    for(size_t j = 0; j < Z[i].size(); ++j)
    {
      if(Z_err_floor[i][j] != 0. && other.Z_err_floor[i][j] != 0.)
      {
        residual += pow((Z[i][j] - other.Z[i][j]).real() / Z_err_floor[i][j], 2.) +
                    pow((Z[i][j] - other.Z[i][j]).imag() / Z_err_floor[i][j], 2.);
        n += 2;
      }

      if(PT_err[i][j] != 0. && other.PT_err[i][j] != 0.)
      {
        residual += pow((PT[i][j] - other.PT[i][j]) / PT_err[i][j], 2.);
        ++n;
      }

      if(i < T.size())
      {
        if(T_err_floor[i][j] != 0. && other.T_err_floor[i][j] != 0.)
        {
          residual += pow((T[i][j] - other.T[i][j]).real() / T_err_floor[i][j], 2.) +
                      pow((T[i][j] - other.T[i][j]).imag() / T_err_floor[i][j], 2.);
          n += 2;
        }
      }
    }
  }

  return sqrt(1.0 / n * residual);
}
