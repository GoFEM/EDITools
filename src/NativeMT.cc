#include "include/NativeMT.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <tuple>

namespace {
void require(bool condition, const std::string &message)
{
  if(!condition) throw std::runtime_error(message);
}
void valid_name(const std::string &name)
{
  require(!name.empty() && name.find_first_not_of(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-") == std::string::npos,
    "Receiver name must use only A-Z, a-z, 0-9, underscore, dot or hyphen: " + name);
}
double number(const std::string &token)
{
  require(!token.empty() && token.find_first_not_of("0123456789+-.eE") == std::string::npos,
          "Invalid numeric token: " + token);
  std::istringstream stream(token);
  stream.imbue(std::locale::classic());
  double value;
  require(bool(stream >> value) && stream.peek() == std::char_traits<char>::eof()
          && std::isfinite(value), "Invalid finite number: " + token);
  return value;
}
std::vector<std::string> fields(const std::string &line)
{
  std::istringstream stream(line);
  stream.imbue(std::locale::classic());
  std::vector<std::string> result;
  std::string token;
  while(stream >> token) result.push_back(token);
  if(!result.empty() && result.front().front() == '#') result.clear();
  return result;
}
std::string numeric(double value)
{
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return stream.str();
}
std::string context(const NativeMT::Observation &row)
{
  const auto &m = NativeMT::mapping(row.type);
  return row.receiver + " at " + numeric(row.frequency) + " Hz, " + m.observable + " " + m.component;
}
bool rho(RealDataType type) { return type == RhoZxx || type == RhoZxy || type == RhoZyx || type == RhoZyy; }
bool phase(RealDataType type) { return type == PhsZxx || type == PhsZxy || type == PhsZyx || type == PhsZyy; }
void validate(const NativeMT::Observation &row)
{
  NativeMT::mapping(row.type);
  valid_name(row.receiver);
  require(std::isfinite(row.frequency) && row.frequency > 0, "Frequency must be finite and positive.");
  require(std::isfinite(row.value), "Non-finite value: " + context(row));
  require(std::isfinite(row.error) && row.error > 0,
          "Standard error must be finite and strictly positive: " + context(row));
  require(!rho(row.type) || row.value >= 0, "Negative apparent resistivity: " + context(row));
}
void validate_keys(const std::vector<NativeMT::Observation> &rows)
{
  std::set<std::tuple<double, std::string, RealDataType>> keys;
  std::set<double> frequencies;
  for(const auto &row: rows) {
    validate(row);
    require(keys.emplace(row.frequency, row.receiver, row.type).second,
            "Duplicate observation: " + context(row));
    frequencies.insert(row.frequency);
  }
  for(auto it = frequencies.begin(); it != frequencies.end(); ++it) {
    auto next = std::next(it);
    if(next != frequencies.end())
      require(!NativeMT::same_frequency(*it, *next), "Near-duplicate frequencies " + numeric(*it) +
              " and " + numeric(*next) + " Hz. Resolve them explicitly before export.");
  }
}
}

namespace NativeMT {
const std::vector<Mapping> &mappings()
{
  static const std::vector<Mapping> table = {
    {RealZxx, "impedance_xx", "real"}, {ImagZxx, "impedance_xx", "imag"},
    {RealZxy, "impedance_xy", "real"}, {ImagZxy, "impedance_xy", "imag"},
    {RealZyx, "impedance_yx", "real"}, {ImagZyx, "impedance_yx", "imag"},
    {RealZyy, "impedance_yy", "real"}, {ImagZyy, "impedance_yy", "imag"},
    {RealTzx, "induction_x", "real"}, {ImagTzx, "induction_x", "imag"},
    {RealTzy, "induction_y", "real"}, {ImagTzy, "induction_y", "imag"},
    {PTxx, "phase_tensor_xx", "value"}, {PTxy, "phase_tensor_xy", "value"},
    {PTyx, "phase_tensor_yx", "value"}, {PTyy, "phase_tensor_yy", "value"},
    {RhoZxx, "apparent_resistivity_xx", "value"}, {RhoZxy, "apparent_resistivity_xy", "value"},
    {RhoZyx, "apparent_resistivity_yx", "value"}, {RhoZyy, "apparent_resistivity_yy", "value"},
    {PhsZxx, "phase_xx", "value"}, {PhsZxy, "phase_xy", "value"},
    {PhsZyx, "phase_yx", "value"}, {PhsZyy, "phase_yy", "value"}
  };
  return table;
}
const Mapping &mapping(RealDataType type)
{
  for(const auto &m: mappings()) if(m.type == type) return m;
  throw std::runtime_error("Unsupported native MT data type: " + std::to_string(int(type)) +
                           ". Log-resistivity is not supported.");
}
bool same_frequency(double a, double b)
{
  return std::abs(a-b) <= 1e-12 * std::max({1., std::abs(a), std::abs(b)});
}
std::vector<Observation> read_observations(std::istream &input)
{
  std::vector<Observation> rows;
  std::string line;
  std::size_t lineno = 0;
  while(std::getline(input, line)) {
    ++lineno;
    try {
      auto f = fields(line);
      if(f.empty()) continue;
      require(f.size() == 6, "Expected exactly six fields (no inline comments).");
      auto it = std::find_if(mappings().begin(), mappings().end(), [&](const Mapping &m) {
        return f[2] == m.observable && f[3] == m.component;
      });
      require(it != mappings().end(), "Unknown observable/component: " + f[2] + " " + f[3]);
      Observation row{number(f[0]), f[1], it->type, number(f[4]), number(f[5])};
      validate(row);
      rows.push_back(row);
    } catch(const std::exception &e) {
      throw std::runtime_error("Observation line " + std::to_string(lineno) + ": " + e.what());
    }
  }
  require(!input.bad(), "Failed reading observations.");
  require(!rows.empty(), "Empty observation file.");
  validate_keys(rows);
  return rows;
}
Receivers read_receivers(std::istream &input)
{
  Receivers receivers;
  std::string line;
  std::size_t lineno = 0;
  while(std::getline(input, line)) {
    ++lineno;
    try {
      auto f = fields(line);
      if(f.empty()) continue;
      require(f.size() == 6 && f[0] == "Dipole" && f[2] == "1",
              "Expected Dipole name 1 x y z. Intersite export is not supported.");
      valid_name(f[1]);
      require(receivers.emplace(f[1], std::array<double, 3>{{number(f[3]), number(f[4]), number(f[5])}}).second,
              "Duplicate receiver: " + f[1]);
    } catch(const std::exception &e) {
      throw std::runtime_error("Receiver line " + std::to_string(lineno) + ": " + e.what());
    }
  }
  require(!input.bad(), "Failed reading receivers.");
  require(!receivers.empty(), "Empty receiver file.");
  return receivers;
}
std::array<double, 3> model_coordinates(const std::array<double, 3> &p,
                                      const std::array<double, 3> &origin)
{
  for(unsigned i = 0; i < 3; ++i)
    require(std::isfinite(p[i]) && std::isfinite(origin[i]), "Coordinates and origin must be finite metres.");
  std::array<double, 3> result{{p[1]-origin[1], p[0]-origin[0], origin[2]-p[2]}};
  for(double v: result) require(std::isfinite(v), "Coordinate conversion overflow.");
  return result;
}
Receivers read_projected_receivers(std::istream &input, const std::array<double, 3> &origin)
{
  Receivers receivers;
  std::string line;
  std::size_t lineno = 0;
  while(std::getline(input, line)) {
    ++lineno;
    try {
      const auto f = fields(line);
      if(f.empty()) continue;
      require(f.size() == 4, "Expected name easting northing elevation, in metres.");
      valid_name(f[0]);
      auto position = model_coordinates({{number(f[1]), number(f[2]), number(f[3])}}, origin);
      require(receivers.emplace(f[0], position).second, "Duplicate receiver: " + f[0]);
    } catch(const std::exception &e) {
      throw std::runtime_error("Projected receiver line " + std::to_string(lineno) + ": " + e.what());
    }
  }
  require(!input.bad(), "Failed reading projected receivers.");
  require(!receivers.empty(), "Empty projected receiver file.");
  return receivers;
}

Export encode(const std::vector<Observation> &rows, const Receivers &receivers,
              const std::string &description, bool distortion)
{
  require(!rows.empty(), "Empty export: no selected valid observations remain.");
  validate_keys(rows);
  Export result;
  result.observations = rows;
  std::set<double> frequencies;
  std::set<std::tuple<std::string, double, RealDataType>> keys;
  for(const auto &row: rows) {
    const auto it = receivers.find(row.receiver);
    require(it != receivers.end(), "Missing model coordinates for receiver: " + row.receiver);
    for(double coordinate: it->second)
      require(std::isfinite(coordinate), "Non-finite model coordinates for receiver: " + row.receiver);
    result.receivers.insert(*it);
    frequencies.insert(row.frequency);
    keys.emplace(row.receiver, row.frequency, row.type);
  }
  result.frequencies.assign(frequencies.begin(), frequencies.end());
  if(distortion) {
    const std::map<RealDataType, RealDataType> companions = {
      {RealZxx, RealZyx}, {ImagZxx, ImagZyx}, {RealZyx, RealZxx}, {ImagZyx, ImagZxx},
      {RealZxy, RealZyy}, {ImagZxy, ImagZyy}, {RealZyy, RealZxy}, {ImagZyy, ImagZxy}};
    for(const auto &row: rows) {
      require(!rho(row.type) && !phase(row.type),
              "Distortion estimate = true does not support apparent resistivity or phase. Disable distortion compatibility.");
      const auto it = companions.find(row.type);
      if(it != companions.end())
        require(keys.count(std::make_tuple(row.receiver, row.frequency, it->second)),
                "Distortion compatibility: missing " + std::string(mapping(it->second).observable) + " " +
                mapping(it->second).component + " companion for " + context(row));
    }
  }
  bool correlated = false, derived = false;
  for(const auto &row: rows) {
    if(rho(row.type) || phase(row.type)) {
      derived = true;
      const unsigned component = type_to_column_table.at(row.type);
      const RealDataType real = static_cast<RealDataType>(int(RealZxx) + 2*component);
      correlated |= keys.count(std::make_tuple(row.receiver, row.frequency, real)) ||
                    keys.count(std::make_tuple(row.receiver, row.frequency, static_cast<RealDataType>(int(real)+1)));
    }
  }
  if(correlated) result.warnings.push_back("Impedance and derived resistivity/phase share measurements; the file cannot represent their error correlations.");
  if(derived) result.warnings.push_back("Use Distortion estimate = false with apparent resistivity or phase observations.");

  // Export ordering follows the table, rather than the numeric data-type enum.
  std::map<RealDataType, unsigned> order;
  for(unsigned i = 0; i < mappings().size(); ++i) order[mappings()[i].type] = i;
  std::sort(result.observations.begin(), result.observations.end(), [&](const Observation &a, const Observation &b) {
    return std::make_tuple(a.receiver, a.frequency, order.at(a.type)) <
           std::make_tuple(b.receiver, b.frequency, order.at(b.type));
  });
  std::ostringstream data, recvs, freqs;
  for(auto stream: {&data, &recvs, &freqs}) {
    stream->imbue(std::locale::classic());
    *stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  }
  data << "# frequency receiver observable component value error\n"
       << "# SI units; x=North y=East z=down; exp(+i omega t); phase in degrees\n";
  recvs << "# kind name point_count x y z\n";
  // Prefix every provenance line: descriptions can never become data records.
  std::istringstream provenance(description);
  std::string line;
  while(std::getline(provenance, line)) {
    for(char &ch: line) if(static_cast<unsigned char>(ch) < 32 || static_cast<unsigned char>(ch) > 126) ch = ' ';
    recvs << "# " << line << '\n';
  }
  for(const auto &row: result.observations) {
    const auto &m = mapping(row.type);
    data << row.frequency << ' ' << row.receiver << ' ' << m.observable << ' ' << m.component << ' '
         << row.value << ' ' << row.error << '\n';
  }
  for(const auto &receiver: result.receivers)
    recvs << "Dipole " << receiver.first << " 1 " << receiver.second[0] << ' ' << receiver.second[1]
          << ' ' << receiver.second[2] << '\n';
  for(double f: result.frequencies) freqs << f << '\n';
  result.data_text = data.str(); result.receiver_text = recvs.str(); result.frequency_text = freqs.str();
  return result;
}
}

// Access only; native export never mutates the survey or its plot conventions.
struct NativeMTStationAccess {
  enum Status { Available, Masked, Missing, UndefinedPhase };
  static Status get(const MTStationData &s, unsigned f, RealDataType type, double &value, double &error)
  {
    const unsigned c = type_to_column_table.at(type);
    const auto &m = NativeMT::mapping(type);
    const std::string observable(m.observable);
    if(observable.find("impedance_") == 0) {
      if(!s.Z_mask.at(c).at(f)) return Masked;
      value = std::string(m.component) == "real" ? s.Z.at(c).at(f).real() : s.Z.at(c).at(f).imag();
      error = s.Z_err_floor.at(c).at(f);
    } else if(observable.find("induction_") == 0) {
      if(!s.T_mask.at(c).at(f)) return Masked;
      value = std::string(m.component) == "real" ? s.T.at(c).at(f).real() : s.T.at(c).at(f).imag();
      // The legacy plot floor replaces tipper errors. Native export retains
      // measured uncertainties when larger than that absolute floor.
      error = s.T_err.at(c).at(f);
      if(!std::isfinite(s.T_err_floor.at(c).at(f)))
        error = s.T_err_floor.at(c).at(f);
      else if(std::isfinite(error) && error >= 0)
        error = std::max(error, s.T_err_floor.at(c).at(f));
    } else if(observable.find("phase_tensor_") == 0) {
      if(!s.PT_mask.at(c).at(f)) return Masked;
      value = s.PT.at(c).at(f); error = s.PT_err.at(c).at(f);
    } else {
      if(!s.Z_mask.at(c).at(f)) return Masked;
      if(rho(type)) {
        value = s.Rho.at(c).at(f); error = s.Rho_err.at(c).at(f);
      } else {
        const auto z = s.Z.at(c).at(f);
        if(std::isnan(z.real()) || std::isnan(z.imag())) return Missing;
        require(std::isfinite(z.real()) && std::isfinite(z.imag()), "Non-finite impedance used to derive phase: " + s.name());
        if(z == std::complex<double>(0., 0.)) return UndefinedPhase;
        value = std::atan2(z.imag(), z.real()) * 180. / std::acos(-1.);
        error = s.Phs_err.at(c).at(f);
      }
    }
    // The EDI reader represents its explicit EMPTY sentinel as NaN.
    return std::isnan(value) ? Missing : Available;
  }
};

namespace NativeMT {
Export prepare(const std::vector<const MTStationData *> &stations, const Options &options)
{
  require(options.horizontal_axes_confirmed,
          "Confirm responses are aligned with the model North/East axes, including sensor azimuth and grid-north corrections.");
  require(options.conventions_confirmed, "Specify verified source time and vertical magnetic conventions.");
  require(!options.coordinate_description.empty(), "Describe the coordinate system, origin and vertical datum.");
  require(!options.types.empty(), "Select at least one data component.");
  require(!options.periods.empty(), "Select at least one period.");
  std::set<RealDataType> types;
  for(auto type: options.types) {
    mapping(type);
    require(types.insert(type).second, "Duplicate selected data type.");
  }
  // Use exact original periods for selection; never merge nearby stations' data.
  std::map<double, double> selected;
  for(double period: options.periods) {
    require(std::isfinite(period) && period > 0, "Selected periods must be finite and positive.");
    const double frequency = 1. / period;
    require(std::isfinite(frequency) && frequency > 0, "Selected period cannot be represented as a frequency.");
    require(selected.emplace(period, frequency).second, "Duplicate selected period.");
  }
  std::vector<Observation> rows;
  std::size_t masked = 0, missing = 0, undefined = 0;
  for(const auto *station: stations) {
    require(station != nullptr, "Invalid station.");
    if(!station->active()) continue;
    const auto &frequencies = station->frequencies();
    for(unsigned f = 0; f < frequencies.size(); ++f) {
      require(std::isfinite(frequencies[f]) && frequencies[f] > 0, "Invalid frequency at station " + station->name());
      const auto selection = selected.find(1. / frequencies[f]);
      if(selection == selected.end()) continue;
      for(auto type: types) {
        double value = 0, error = 0;
        auto status = NativeMTStationAccess::get(*station, f, type, value, error);
        if(status == NativeMTStationAccess::Masked) { ++masked; continue; }
        if(status == NativeMTStationAccess::Missing) { ++missing; continue; }
        if(status == NativeMTStationAccess::UndefinedPhase) { ++undefined; continue; }
        const auto &m = mapping(type);
        if(options.negative_time_convention &&
           (std::string(m.component) == "imag" || phase(type) || std::string(m.observable).find("phase_tensor_") == 0))
          value = -value;
        if(options.vertical_magnetic_up && std::string(m.observable).find("induction_") == 0)
          value = -value;
        rows.push_back({selection->second, station->name(), type, value, error});
      }
    }
  }
  Export result = encode(rows, options.receivers, options.coordinate_description, options.distortion);
  result.masked = masked; result.missing = missing; result.undefined_phase = undefined;
  return result;
}
}
