#include "include/MTResponseData.h"
#include "include/NativeMT.h"
#include <cmath>
#include <iterator>
#include <locale>
#include <sstream>
#include <tuple>

namespace {
std::vector<MTResponseData::Scalar> read_gofem(std::istream &input)
{
  std::vector<MTResponseData::Scalar> rows;
  std::set<std::tuple<double, std::string, RealDataType>> keys;
  std::string line;
  unsigned lineno = 0;
  while(std::getline(input, line)) {
    ++lineno;
    try {
      std::istringstream stream(line);
      stream.imbue(std::locale::classic());
      std::string type, source, extra;
      if(!(stream >> type) || type.front() == '#') continue;
      MTResponseData::Scalar row;
      if(!(stream >> row.frequency >> source >> row.receiver >> row.value >> row.error) ||
         (stream >> extra && extra.front() != '#'))
        throw std::runtime_error("Expected type frequency source receiver value error.");
      row.type = Datum::convert_string_to_type(type);
      if(!type_to_column_table.count(row.type))
        throw std::runtime_error("Unsupported MT response type: " + type);
      if(!std::isfinite(row.frequency) || row.frequency <= 0. || !std::isfinite(row.value) ||
         !std::isfinite(row.error) || row.error < 0.)
        throw std::runtime_error("Frequency must be positive; values and nonnegative errors must be finite.");
      if((row.type == RhoZxx || row.type == RhoZxy || row.type == RhoZyx || row.type == RhoZyy) && row.value < 0.)
        throw std::runtime_error("Apparent resistivity must be nonnegative.");
      if(!keys.emplace(row.frequency, row.receiver, row.type).second)
        throw std::runtime_error("Duplicate response for " + row.receiver + ", " + type + " at this frequency.");
      rows.push_back(row);
    } catch(const std::exception &error) {
      throw std::runtime_error("GoFEM line " + std::to_string(lineno) + ": " + error.what());
    }
  }
  if(input.bad()) throw std::runtime_error("Failed reading GoFEM responses.");
  if(rows.empty()) throw std::runtime_error("Empty GoFEM response file.");

  return rows;
}
}

namespace MTResponseData {
std::vector<Scalar> read(std::istream &input, Format format)
{
  if(format == Format::Auto) {
    // Inspect content, not extensions, so mixed selections and extensionless
    // solver output use the same import action. Replay also supports non-seekable streams.
    const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if(input.bad()) throw std::runtime_error("Failed reading responses.");
    std::istringstream probe(text);
    std::string line, first;
    while(std::getline(probe, line)) {
      std::istringstream fields(line);
      if(fields >> first && first.front() != '#') break;
      first.clear();
    }
    if(first.empty()) throw std::runtime_error("Empty response file.");
    const auto detected = first.find_first_of("0123456789+-.") == 0 ? Format::Native : Format::GoFEM;
    std::istringstream replay(text);
    return read(replay, detected);
  }
  return format == Format::Native ? NativeMT::read_observations(input) : read_gofem(input);
}

std::map<std::string, MTStationData> stations(const std::vector<Scalar> &rows)
{
  std::map<std::string, std::map<double, std::vector<Scalar>>> grouped;
  for(const auto &row: rows)
    grouped[row.receiver][row.frequency].push_back(row);

  std::map<std::string, MTStationData> stations;
  for(const auto &receiver: grouped) {
    auto &station = stations[receiver.first];
    std::set<double> frequencies;
    for(const auto &entry: receiver.second) frequencies.insert(entry.first);
    station.set_name(receiver.first);
    station.set_size(frequencies.size(), true);
    station.set_frequencies(frequencies);

    for(const auto &entry: receiver.second) {
      std::vector<RealDataType> types;
      std::vector<double> values, errors;
      for(const auto &row: entry.second) {
        types.push_back(row.type);
        values.push_back(row.value);
        errors.push_back(row.error);
      }
      station.set_data(entry.first, types, values, errors);
    }

    station.calculate_apparent_resistivity();
    station.calculate_phase();
    station.calculate_phase_tensor();

    // Derived curves fill only what can be computed from complete impedances.
    // Restore explicitly supplied scalars after deriving the other curves.
    for(const auto &entry: receiver.second)
      for(const auto &row: entry.second)
        if((row.type >= RhoZxx && row.type <= PhsZyy) || (row.type >= PTxx && row.type <= PTyy))
          station.set_data(row.frequency, {row.type}, {row.value}, {row.error});
  }
  return stations;
}
}
