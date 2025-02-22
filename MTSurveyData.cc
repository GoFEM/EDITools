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

#include "MTSurveyData.h"

#include "datum.h"

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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

void MTSurveyData::load_from_gofem(std::string file_path)
{
  std::ifstream ifs(file_path);

  if(!ifs.is_open())
    throw std::runtime_error("Cannot open file " + file_path);

  std::vector<std::vector<std::string>> columns(6);

  std::string line;
  while (!ifs.eof())
  {
    std::getline (ifs, line);
    // Trim string
    boost::trim_copy(line);

    // Skip empty lines and comments
    if (line.length() < 1 || (line[0] == '#' && line[1] == '!'))
      continue;

    std::vector<std::string> strs;
    boost::split(strs, line, boost::is_any_of("\t "), boost::token_compress_on);

    if(strs.size() != columns.size())
      throw std::runtime_error("Wrong string format: " + line);

    for(unsigned i = 0; i < columns.size(); ++i)
      columns[i].push_back(strs[i]);
  }

  // Extract unique station names
  std::set<std::string> station_names(columns[3].begin(), columns[3].end());

  auto rho_func = [](RealDataType type)
                    { return (type == RhoZxx | type == RhoZxy |
                              type == RhoZyx | type == RhoZyy); };

  auto phase_func = [](RealDataType type)
                    { return (type == PhsZxx | type == PhsZxy |
                              type == PhsZyx | type == PhsZyy); };

  auto pt_func = [](RealDataType type)
                    { return (type == PTxx | type == PTxy |
                              type == PTyx | type == PTyy); };

  // Fill stations with data
  for(auto sname: station_names)
  {
    MTStationData sdata;
    sdata.set_name(sname);

    std::vector<unsigned> sindices;
    std::set<double> frequencies;
    for(unsigned i = 0; i < columns[3].size(); ++i)
    {
      if(sname == columns[3][i])
      {
        sindices.push_back(i);

        double frequency = boost::lexical_cast<double>(columns[1][i]);
        frequencies.insert(frequency);
      }
    }

    sdata.set_size(frequencies.size());

    sdata.set_frequencies(frequencies);

    for(double frequency: frequencies)
    {
      std::vector<RealDataType> types;
      std::vector<double> values, errors;

      for(unsigned i = 0; i < sindices.size(); ++i)
      {
        const double f = boost::lexical_cast<double>(columns[1][sindices[i]]);
        if(fabs(frequency - f) < 1e-10)
        {
          RealDataType type = Datum::convert_string_to_type(columns[0][sindices[i]]);
          const double value = boost::lexical_cast<double>(columns[4][sindices[i]]);
          const double error = boost::lexical_cast<double>(columns[5][sindices[i]]);

          types.push_back(type);
          values.push_back(value);
          errors.push_back(error);
        }
      }

      sdata.set_data(frequency, types, values, errors);

      if (std::none_of(types.cbegin(), types.cend(), rho_func))
        sdata.calculate_apparent_resistivity();

      if (std::none_of(types.cbegin(), types.cend(), phase_func))
        sdata.calculate_phase();

      if (std::none_of(types.cbegin(), types.cend(), pt_func))
        sdata.calculate_phase_tensor();
    }

    m_stations_data.insert(std::make_pair(sname, sdata));
  }
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

std::vector<std::array<double, 3>> MTSurveyData::get_stations_locations() const
{
  std::vector<std::array<double, 3>> locations;

  for(auto &station: m_stations_data)
  {
    const MTStationData &data = station.second;
    locations.push_back(data.position());
  }

  return locations;
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

  return locations;
}

std::string MTSurveyData::closest_station_name(const double &lat,
                                               const double &lon) const
{
  std::map<double, std::string> dist_map;

  for(const auto &station: m_stations_data)
  {
    const MTStationData &data = station.second;
    const double &lat_i = data.position()[0];
    const double &lon_i = data.position()[1];
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
  std::set<double, weak_compare> periods_set;

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
  // Write data file
  {
    std::ofstream ofs(file);

    if(!ofs.is_open())
      throw std::runtime_error("Cannot open file " + file);

    ofs << "# DataType Frequency SourceName ReceiverName Value Error" << std::endl;

    for(auto &station: m_stations_data)
    {
      const MTStationData &data = station.second;

      if(!data.active())
        continue;

      data.write(ofs, types, periods);
    }
  }

  // Write receiver file
  {
    std::ofstream ofs(file + ".recvs");

    if(!ofs.is_open())
      throw std::runtime_error("Cannot open file " + file);

    ofs << "# Type Name Electrodes x y z" << std::endl;

    for(auto &station: m_stations_data)
    {
      const MTStationData &data = station.second;

      if(!data.active())
        continue;

      ofs << "Dipole\t" << data.name() << "\t1\t"
          << data.position()[0] << "\t"
          << data.position()[1] << "\t"
          << data.position()[2] << std::endl;
    }
  }

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
