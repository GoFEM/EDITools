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

#ifndef MT_SURVEY_DATA_H
#define MT_SURVEY_DATA_H

#include <string>
#include <vector>
#include <map>

#include <boost/serialization/map.hpp>

#include "EDIFileReader.h"

class MTSurveyData
{
public:
  MTSurveyData();
  MTSurveyData(const std::string survey_name);

  std::string get_survey_name() const;

  bool is_station_present(const std::string &name) const;

  // Returns stations names which have duplicates (and therefore were ignored)
  std::vector<std::string> load_from_edi(std::vector<std::string> &file_list);
  void load_from_gofem(std::string file_path);

  std::vector<std::string> get_stations_names() const;

  MTStationData &get_station_data(const std::string &name);

  std::vector<std::array<double, 3>> get_stations_locations() const;
  std::vector<std::array<double, 3>> get_stations_locations(const std::vector<std::string> &names) const;

  std::string closest_station_name(const double &lat, const double &lon) const;

  void set_active_flag(const std::string &name, const bool flag);
  void set_active_flag(const std::string &name, RealDataType type, const bool flag);
  bool is_active(const std::string &name) const;

  // Decimates data w.r.t. periods, i.e. keep every second period
  void decimate();

  std::vector<double> get_unique_periods() const;

  void write_gofem(const std::string &file,
                   const std::vector<RealDataType> &types,
                   const std::vector<double> &periods) const;

  void set_error_floor(double error_floor);

  unsigned n_stations() const;

  std::set<RealDataType> get_active_types(const std::string &name) const;

  void remove_station(const std::string &name);
  void rename_station(const std::string &name, const std::string &new_name);

private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & m_survey_name;
      ar & m_stations_data;
  }

private:
  std::string m_survey_name;
  std::map<std::string, MTStationData> m_stations_data;
};

#endif // MT_SURVEY_DATA_H
