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

#ifndef EDI_FILE_DATA_H
#define EDI_FILE_DATA_H

#include <string>
#include <map>
#include <vector>
#include <set>

#include <Eigen/Dense>

#include "include/MTStationData.h"

using StringMap = std::map<std::string, std::string>;
using Matrix2cd = Eigen::Matrix<std::complex<double>, 2, 2, Eigen::RowMajor>;
using Matrix2d = Eigen::Matrix<double, 2, 2, Eigen::RowMajor>;

struct DEFINEMEAS_DATA
{
  StringMap options;
  std::vector<StringMap> MEAS;
};

struct SPECTRA_DATA
{
  StringMap options;
  Eigen::Matrix<double, 7, 7, Eigen::RowMajor> data;
};

struct MTSECT_DATA
{
  StringMap options;
  std::vector<double> frequencies;
  std::vector<std::vector<double>> Zdata;
  std::vector<std::vector<double>> Zvar;
  std::vector<std::vector<double>> Tdata;
  std::vector<std::vector<double>> Tvar;
};


struct SPECTRASECT_DATA
{
  StringMap options;
  std::vector<std::string> channel_ids;
  std::vector<SPECTRA_DATA> spectra_data;
};

enum BlockID
{
  EDI_HEAD, EDI_INFO, EDI_DEFINEMEAS, EDI_MTSECT, EDI_SPECTRASECT, EDI_END, EDI_INVALID
};

const std::map<std::string, BlockID> block_names =
      {{"HEAD", EDI_HEAD},
       {"INFO", EDI_INFO},
       {"DEFINEMEAS", EDI_DEFINEMEAS},
       {"MTSECT", EDI_MTSECT},
       {"SPECTRASECT", EDI_SPECTRASECT},
       {"END", EDI_END}};

class EDIFileReader
{
public:
  EDIFileReader(const std::string edi_file_name);

  const MTStationData &get_mt_data() const;

private:
  void read_edi();
  BlockID get_block_id(const std::string &line) const;

  void read_head_block(std::istringstream &ss);
  void read_info_block(std::istringstream &ss);
  void read_definemeas_section(std::istringstream &ss);
  void read_spectrasect_section(std::istringstream &ss);
  void read_mtsect_section(std::istringstream &ss);

  std::string trim(const std::string& str, const std::string& whitespace);

  template<typename T>
  T get_option_value(const StringMap &options,
                     const std::string &option_name) const;

  void calculate_data_from_spectra();
  void calculate_data_from_mtsect();
  void fill_stations_location();

private:
  std::string m_edi_file_name;

  std::string edi_info;
  StringMap head_options;
  DEFINEMEAS_DATA definemeas;
  SPECTRASECT_DATA spectrasect;
  MTSECT_DATA mtsect;

  double empty_value;

  bool is_data_spectra;

  MTStationData station_data;
};

#endif // EDI_FILE_DATA_H
