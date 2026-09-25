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
#include <memory>
#include <boost/serialization/shared_ptr.hpp>

#include <boost/serialization/map.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/vector.hpp>

#include "MTStationData.h"
#include "SurveyCoordinates.h"
#include "MTResponseData.h"
#include "PeriodResamplingInfo.h"

class MTSurveyData
{
public:
  MTSurveyData() = default;
  explicit MTSurveyData(const std::string &survey_name);

  std::string get_survey_name() const;
  const PeriodResampling::Info &resampling_info() const { return m_resampling_info; }
  const MTSurveyData *resampling_source() const { return m_resampling_source.get(); }

  bool is_station_present(const std::string &name) const;

  // Returns stations names which have duplicates (and therefore were ignored)
  std::vector<std::string> load_from_edi(const std::vector<std::string> &file_list);
  void load_responses(const std::string &file_path, MTResponseData::Format format = MTResponseData::Format::Auto);
  void load_from_gofem(const std::string &file_path);
  void load_from_native_responses(const std::string &file_path);

  std::vector<std::string> get_stations_names() const;

  MTStationData &get_station_data(const std::string &name);
  const MTStationData &get_station_data(const std::string &name) const;
  const std::vector<MTResponseData::Scalar> &response_observations() const { return m_response_observations; }

  std::vector<std::array<double, 3>> get_stations_locations() const;
  std::vector<std::array<double, 3>> get_stations_locations(const std::vector<std::string> &names) const;
  std::vector<std::array<double, 3>> geographic_locations() const;
  const SurveyCoordinates &coordinates() const { return m_coordinates; }
  void set_coordinates(const SurveyCoordinates &coordinates);
  void ensure_utm_coordinates();

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
  friend struct PeriodResamplingAccess;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & m_survey_name;
      ar & m_stations_data;
      if(version >= 1) ar & m_coordinates;
      else if(Archive::is_loading::value) m_coordinates = SurveyCoordinates{};
      if(version >= 2) ar & m_response_observations;
      else if(Archive::is_loading::value) m_response_observations.clear();
      if(version >= 3) ar & m_resampling_info & m_resampling_source;
      else if(Archive::is_loading::value) { m_resampling_info = {}; m_resampling_source.reset(); }
  }

private:
  std::string m_survey_name;
  std::map<std::string, MTStationData> m_stations_data;
  SurveyCoordinates m_coordinates;
  std::vector<MTResponseData::Scalar> m_response_observations;
  PeriodResampling::Info m_resampling_info;
  std::shared_ptr<MTSurveyData> m_resampling_source;
};

BOOST_CLASS_VERSION(MTSurveyData, 3)

#endif // MT_SURVEY_DATA_H
