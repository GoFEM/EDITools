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

#include "include/MTSurveyData.h"

#include "include/datum.h"

#include <fstream>
#include <iomanip>
#include <locale>

MTSurveyData::MTSurveyData()
{}

MTSurveyData::MTSurveyData(const std::string survey_name):
  m_survey_name(survey_name)
{}

std::string MTSurveyData::get_survey_name() const
{
  return m_survey_name;
}

bool MTSurveyData::is_station_present(const std::string &name) const
{
  auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    return false;
  else
    return true;
}

std::vector<std::string> MTSurveyData::load_from_edi(std::vector<std::string> &file_list)
{
  std::vector<std::string> duplicates;

  for(auto& file: file_list)
  {
    EDIFileReader edi_file_reader(file);
    const MTStationData &data = edi_file_reader.get_mt_data();

    const auto it = m_stations_data.find(data.name());
    if(it == m_stations_data.end())
    {
      m_stations_data.insert({data.name(), data});
      m_stations_data[data.name()].set_error_floor(0.);
    }
    else
      duplicates.push_back(data.name() + " (" + file + ")");
  }

  return duplicates;
}

void MTSurveyData::load_responses(const std::string &file_path, MTResponseData::Format format)
{
  std::ifstream input(file_path);
  if(!input.is_open()) throw std::runtime_error("Cannot open file " + file_path);
  auto rows = MTResponseData::read(input, format);
  auto stations = MTResponseData::stations(rows);
  m_stations_data = std::move(stations);
  m_response_observations = std::move(rows);
  m_survey_name = file_path;
  m_coordinates = SurveyCoordinates{};
}

void MTSurveyData::load_from_gofem(std::string file_path)
{
  load_responses(file_path, MTResponseData::Format::GoFEM);
}

void MTSurveyData::load_from_native_responses(const std::string &file_path)
{
  load_responses(file_path, MTResponseData::Format::Native);
}

std::vector<std::string> MTSurveyData::get_stations_names() const
{
  std::vector<std::string> names;
  for(const auto &station: m_stations_data)
    names.push_back(station.first);

  return names;
}

MTStationData &MTSurveyData::get_station_data(const std::string &name)
{
  auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  return it->second;
}

const MTStationData &MTSurveyData::get_station_data(const std::string &name) const
{
  return m_stations_data.at(name);
}

std::vector<std::array<double, 3>> MTSurveyData::get_stations_locations() const
{
  return m_coordinates.transform(geographic_locations());
}

std::vector<std::array<double, 3>> MTSurveyData::geographic_locations() const
{
  std::vector<std::array<double, 3>> locations;

  for(auto &station: m_stations_data)
  {
    const MTStationData &data = station.second;
    locations.push_back(data.position());
  }

  return locations;
}

void MTSurveyData::set_coordinates(const SurveyCoordinates &coordinates)
{
  coordinates.transform(geographic_locations()); // Validate before changing the survey.
  m_coordinates = coordinates;
}

void MTSurveyData::ensure_utm_coordinates()
{
  if(!m_coordinates.utm) set_coordinates(SurveyCoordinates::suggested(geographic_locations()));
}

std::vector<std::array<double, 3> > MTSurveyData::get_stations_locations(const std::vector<std::string> &names) const
{
  std::vector<std::array<double, 3>> locations;

  for(auto &name: names)
  {
    const auto it = m_stations_data.find(name);
    if(it == m_stations_data.end())
      throw std::runtime_error("Station " + name + " not found");

    const MTStationData &data = it->second;

    locations.push_back(data.position());
  }

  return m_coordinates.transform(locations);
}

std::string MTSurveyData::closest_station_name(const double &lat,
                                               const double &lon) const
{
  std::map<double, std::string> dist_map;
  const auto locations = get_stations_locations();
  unsigned index = 0;

  for(const auto &station: m_stations_data)
  {
    const MTStationData &data = station.second;
    const double lat_i = locations[index][0];
    const double lon_i = locations[index++][1];
    const double distance = sqrt((lat - lat_i)*(lat-lat_i) + (lon-lon_i)*(lon-lon_i));
    dist_map.insert(std::make_pair(distance, data.name()));
  }

  return dist_map.begin()->second;
}

void MTSurveyData::set_active_flag(const std::string &name, const bool flag)
{
  auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  MTStationData &data = it->second;
  data.set_active(flag);
}

void MTSurveyData::set_active_flag(const std::string &name, RealDataType type, const bool flag)
{
  auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  MTStationData &data = it->second;
  data.mask_type(type, flag);
}

bool MTSurveyData::is_active(const std::string &name) const
{
  const auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  const MTStationData &data = it->second;

  return data.active();
}

void MTSurveyData::decimate()
{
  for(auto &station: m_stations_data)
  {
    MTStationData &data = station.second;

    data.decimate();
  }
}

std::vector<double> MTSurveyData::get_unique_periods() const
{
  // Preserve distinct samples even at short periods; absolute-tolerance ordering
  // can collapse neighboring points on a common logarithmic grid.
  std::set<double> periods_set;

  for(auto &station: m_stations_data)
  {
    const MTStationData &data = station.second;

    for(auto f: data.frequencies())
      periods_set.insert(1./f);
  }

  std::vector<double> periods;
  std::copy(periods_set.begin(), periods_set.end(), std::back_inserter(periods));

  return periods;
}

void MTSurveyData::write_gofem(const std::string &file,
                               const std::vector<RealDataType> &types,
                               const std::vector<double> &periods) const
{
  std::set<double> written_frequencies;
  const auto geographic = geographic_locations();
  const auto output_coordinates = m_coordinates.utm ? m_coordinates : SurveyCoordinates::suggested(geographic);
  const auto locations = output_coordinates.transform(geographic);
  // Write data file
  {
    std::ofstream ofs(file);

    if(!ofs.is_open())
      throw std::runtime_error("Cannot open file " + file);
    ofs.imbue(std::locale::classic());

    ofs << "# DataType Frequency SourceName ReceiverName Value Error" << std::endl;

    for(auto &station: m_stations_data)
    {
      const MTStationData &data = station.second;

      if(!data.active())
        continue;

      data.write(ofs, types, periods, &written_frequencies);
    }
    if(!ofs) throw std::runtime_error("Cannot write data file " + file);
  }

  // Write receiver file
  {
    std::ofstream ofs(file + ".recvs");

    if(!ofs.is_open())
      throw std::runtime_error("Cannot open file " + file);
    ofs.imbue(std::locale::classic());
    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);

    ofs << "# Type Name Electrodes x y z" << std::endl;
    {
      std::istringstream provenance(output_coordinates.description());
      std::string line;
      while(std::getline(provenance, line)) ofs << "# " << line << '\n';
    }
    unsigned index = 0;

    for(auto &station: m_stations_data)
    {
      const MTStationData &data = station.second;
      const auto &position = locations[index++];

      if(!data.active())
        continue;

      ofs << "Dipole\t" << data.name() << "\t1\t"
          << position[0] << "\t"
          << position[1] << "\t"
          << position[2] << std::endl;
    }
  }
  std::ofstream frequencies(file + ".freqs");
  frequencies.imbue(std::locale::classic());
  frequencies << std::setprecision(std::numeric_limits<double>::max_digits10);
  for(double f: written_frequencies) frequencies << f << '\n';
  if(!frequencies) throw std::runtime_error("Cannot write frequency file " + file + ".freqs");
}

void MTSurveyData::set_error_floor(double error_floor)
{
  for(auto &station: m_stations_data)
  {
    MTStationData &data = station.second;
    data.set_error_floor(error_floor);
  }
}

unsigned MTSurveyData::n_stations() const
{
  return m_stations_data.size();
}

std::set<RealDataType> MTSurveyData::get_active_types(const std::string &name) const
{
  const auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  const MTStationData &data = it->second;

  if(!data.active())
    return std::set<RealDataType>();

  return data.active_types();
}

void MTSurveyData::remove_station(const std::string &name)
{
  m_stations_data.erase(m_stations_data.find(name));
}

void MTSurveyData::rename_station(const std::string &name, const std::string &new_name)
{
  auto it = m_stations_data.find(name);
  if(it == m_stations_data.end())
    throw std::runtime_error("Station " + name + " not found");

  MTStationData &data = it->second;
  data.set_name(new_name);

  std::swap(m_stations_data[new_name], data);
  m_stations_data.erase(it);
}
