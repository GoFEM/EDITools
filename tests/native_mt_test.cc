#include "include/NativeMT.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>

using namespace NativeMT;
namespace {
void check(bool ok, const std::string &message)
{
  if(!ok) throw std::runtime_error(message);
}
void near(double a, double b)
{
  check(std::abs(a-b) <= 1e-13 * std::max({1., std::abs(a), std::abs(b)}), "Unexpected numeric value");
}
void rejects(const std::function<void()> &action, const std::string &part)
{
  try { action(); }
  catch(const std::exception &e) {
    check(std::string(e.what()).find(part) != std::string::npos,
          "Wrong diagnostic: " + std::string(e.what()) + "; wanted " + part);
    return;
  }
  throw std::runtime_error("Expected rejection: " + part);
}
const Receivers receivers{{"S01", {{1.5, 1.5, 1.5}}}, {"S02", {{2., 2., 2.}}}};
std::vector<Observation> parse(const std::string &text)
{
  std::istringstream input(text);
  return read_observations(input);
}
std::vector<Observation> example()
{
  std::ifstream input(NATIVE_MT_TEST_DATA "/native_mt_example.data");
  return read_observations(input);
}
MTStationData station(const std::string &name = "S01", const std::set<double> &frequencies = {1.})
{
  MTStationData s;
  s.set_name(name); s.set_size(frequencies.size()); s.set_frequencies(frequencies);
  for(double f: frequencies) {
    s.set_data(f, {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy,
                   RealTzx, ImagTzx, RealTzy, ImagTzy},
                  {.001, .001, .002, .002, -.002, -.002, .001, .001, 0., .1, -.2, 0.},
                  {.0001, .0001, .0001, .0001, .0001, .0001, .0001, .0001, .01, .01, .02, .02});
  }
  s.calculate_apparent_resistivity(); s.calculate_phase(); s.calculate_phase_tensor();
  s.set_error_floor(0.);
  return s;
}
Options options()
{
  Options o;
  for(const auto &m: mappings()) o.types.push_back(m.type);
  o.periods = {1.}; o.receivers = receivers;
  o.coordinate_description = "Test model frame: origin E=0 N=0 elevation=0, local datum";
  o.horizontal_axes_confirmed = true; o.conventions_confirmed = true;
  return o;
}
const Observation &find(const Export &result, RealDataType type)
{
  for(const auto &row: result.observations) if(row.type == type) return row;
  throw std::runtime_error("Missing expected type");
}

void format_roundtrip()
{
  auto rows = example();
  check(rows.size() == 24, "Expected 24 example scalars");
  for(unsigned i = 0; i < 24; ++i) check(rows[i].type == mappings()[i].type, "Incorrect mapping order");
  auto result = encode(rows, receivers, "Example");
  auto reread = parse(result.data_text);
  check(reread.size() == rows.size(), "Round trip dropped rows");
  for(unsigned i = 0; i < rows.size(); ++i) {
    check(rows[i].type == reread[i].type && rows[i].receiver == reread[i].receiver &&
          rows[i].frequency == reread[i].frequency && rows[i].value == reread[i].value &&
          rows[i].error == reread[i].error, "Round trip changed key, units or precision");
  }
  check(result.receivers.size() == 1 && result.frequencies == std::vector<double>{1.}, "Wrong companion coverage");
  check(result.data_text.find("Plane_wave") == std::string::npos, "Legacy source token leaked");
  // Predicted rows can be shuffled and have unpaired complex scalars.
  std::reverse(rows.begin(), rows.end());
  rows.erase(rows.begin() + 5);
  for(auto &row: rows) row.value *= .5;
  auto predicted = parse(encode(rows, receivers, "Prediction").data_text);
  check(predicted.size() == 23, "Prediction pairing lost sparse rows");
  for(const auto &row: rows) {
    auto it = std::find_if(predicted.begin(), predicted.end(), [&](const Observation &p) {return p.type == row.type;});
    check(it != predicted.end() && it->value == row.value && it->error == row.error, "Prediction reverse mapping failed");
  }
  auto phase_rows = parse("1 S01 phase_yx value 225 2\n1 S01 phase_xy value 0 3\n");
  check(parse(encode(phase_rows, receivers, "Phases").data_text)[1].value == 225, "Unwrapped observed phase changed");
}
void invalid_format()
{
  for(const std::string &line: {
    "1 S01 impedance_xy real 1 0", "1 S01 impedance_xy real 1 -1",
    "0 S01 impedance_xy real 1 1", "-1 S01 impedance_xy real 1 1",
    "1 S01 impedance_xy imaginary 1 1", "1 S01 impedance_xy value 1 1",
    "1 S01 phase_xy real 1 1", "1 S01 log10RhoZxy value 1 1",
    "1 S01 apparent_resistivity_xy value -1 1", "1 S01 impedance_xy real nan 1",
    "1 S01 impedance_xy real 1 inf", "1D0 S01 impedance_xy real 1 1",
    "1,0 S01 impedance_xy real 1 1", "1 S01 impedance_xy real 1 1 extra",
    "1 S01 impedance_xy real 1 1 # comment", "1 S#01 impedance_xy real 1 1",
    "1 S01 impedance_xy real 1", "\xEF\xBB\xBF" "1 S01 impedance_xy real 1 1"})
    rejects([&] {parse(line);}, "line 1");
  rejects([] {parse(" # only a comment\n");}, "Empty");
  rejects([] {parse("1 S01 phase_xy value 1 2\n1 S01 phase_xy value 1 2");}, "Duplicate");
  auto rows = example();
  rows[0].receiver = "unknown";
  rejects([&] {encode(rows, receivers, "");}, "Missing model coordinates");
  rows = example(); rows[0].type = log10RhoZxx;
  rejects([&] {encode(rows, receivers, "");}, "Unsupported");
  rows = example(); rows[0].value = std::numeric_limits<double>::infinity();
  rejects([&] {encode(rows, receivers, "");}, "Non-finite");
  rows[0].value = std::numeric_limits<double>::quiet_NaN();
  rejects([&] {encode(rows, receivers, "");}, "Non-finite");
  auto invalid_receivers = receivers;
  invalid_receivers["S01"][0] = std::numeric_limits<double>::infinity();
  rejects([&] {encode(example(), invalid_receivers, "");}, "Non-finite model coordinates");
  rejects([&] {encode({}, receivers, "");}, "Empty export");
}
struct DecimalComma: std::numpunct<char> {
  char do_decimal_point() const override {return ',';}
  char do_thousands_sep() const override {return '.';}
  std::string do_grouping() const override {return "\3";}
};
void precision_and_comments()
{
  const auto saved = std::locale();
  std::locale::global(std::locale(saved, new DecimalComma));
  try {
    auto rows = parse(" \t# comment\r\n\r\n1.2345678901234567 S01 phase_xy value -135.12345678901235 1.2345678901234567\r\n");
    auto result = encode(rows, receivers, "first line\nsecond line");
    auto reread = parse(result.data_text);
    check(reread[0].frequency == rows[0].frequency && reread[0].value == rows[0].value && reread[0].error == rows[0].error,
          "Lost double precision");
    check(result.data_text.find(',') == std::string::npos, "Locale-dependent decimal");
    check(result.frequency_text.find("1.2345678901234567") != std::string::npos, "Frequency text differs");
    check(result.receiver_text.find("# second line") != std::string::npos, "Uncommented provenance");
  } catch(...) {std::locale::global(saved); throw;}
  std::locale::global(saved);
}
void coordinates()
{
  auto converted = model_coordinates({{501234., 6123456., 320.}}, {{500000., 6123000., 400.}});
  check(converted == std::array<double,3>{{456., 1234., 80.}}, "Wrong N/E order or depth sign");
  std::istringstream projected("# name E N elevation\r\nS01 501234 6123456 320\r\n");
  check(read_projected_receivers(projected, {{500000., 6123000., 400.}}).at("S01") == converted, "Projected reader");
  std::istringstream direct("  # receivers\r\nDipole S01 1 456 1234 80\r\n");
  check(read_receivers(direct).at("S01") == converted, "Receiver reader");
  for(const auto &text: {"Dipole S01 1 0 0 nan", "Dipole S01 2 0 0 0", "Intersite S01 2 0 0 0 1 1 1",
                         "Dipole S01 1 0 0 0 extra", "Dipole S01 1 0 0 0\nDipole S01 1 1 1 1"})
    rejects([&] {std::istringstream stream(text); read_receivers(stream);}, "Receiver line");
}
void frequencies_and_distortion()
{
  auto rows = example();
  rows.push_back({2., "S02", RealZxy, 0., .1});
  rows.push_back({.1, "S02", ImagTzx, -.1, .1});
  auto result = encode(rows, receivers, "Sparse");
  check(result.frequencies == std::vector<double>({.1,1.,2.}), "Sparse frequency union");
  rows.back().frequency = 1. + 5e-13;
  rejects([&] {encode(rows, receivers, "");}, "Near-duplicate frequencies");
  auto full = example(); full.resize(16);
  encode(full, receivers, "", true);
  full.erase(full.begin());
  rejects([&] {encode(full, receivers, "", true);}, "missing impedance_xx real");
  rejects([&] {encode(example(), receivers, "", true);}, "does not support apparent resistivity or phase");
  encode({{1., "S01", PTxy, 0., .1}}, receivers, "", true);
  encode({{1., "S01", RhoZxy, 0., .1}}, receivers, "");
}
void station_export()
{
  auto s = station(); auto o = options();
  auto result = prepare({&s}, o);
  check(result.observations.size() == 24, "Station export missing scalar types");
  near(find(result, RealZxy).value, .002); near(find(result, RealZxy).error, .0001);
  near(find(result, PhsZxy).value, 45.); near(find(result, PhsZyx).value, -135.);
  near(find(result, PhsZxy).error, 180. / std::acos(-1.) * .0001 / std::abs(std::complex<double>(.002,.002)));
  near(find(result, RealTzx).value, 0.); near(find(result, RealTzx).error, .01);
  near(find(result, RhoZxy).value, 1.0132118364233778);
  o.negative_time_convention = true; o.vertical_magnetic_up = true;
  auto changed = prepare({&s}, o);
  near(find(changed, ImagZxy).value, -.002); near(find(changed, PhsZyx).value, 135.);
  near(find(changed, PTxx).value, -find(result, PTxx).value);
  near(find(changed, RealTzy).value, .2); near(find(changed, ImagTzx).value, .1);
  near(find(changed, RhoZxy).value, find(result, RhoZxy).value);
  for(const auto &row: changed.observations) near(row.error, find(result, row.type).error);
  o = options(); o.types = {RealZxy, PhsZxy, RhoZxy};
  s.set_error_floor(.1);
  result = prepare({&s}, o);
  near(find(result, RealZxy).error, std::abs(std::complex<double>(.002,.002)) * .1);
  near(find(result, PhsZxy).error, 180. / std::acos(-1.) * .1);
  o.horizontal_axes_confirmed = false;
  rejects([&] {prepare({&s}, o);}, "Confirm responses");
  o.horizontal_axes_confirmed = true; o.conventions_confirmed = false;
  rejects([&] {prepare({&s}, o);}, "Specify verified source");
}
void masks_and_missing()
{
  auto s1 = station("S01", {1., 2.}); auto s2 = station("S02", {2., 4.});
  auto o = options(); o.types = {RealZxy, ImagTzx}; o.periods = {1., .5, .25};
  s1.set_data_mask(RealZxy, 1., false);
  auto result = prepare({&s2, &s1}, o);
  check(result.observations.size() == 7 && result.masked == 1, "Mask or sparse station selection ignored");
  check(result.frequencies == std::vector<double>({1.,2.,4.}), "Station frequency coverage");
  s2.set_active(false); o.periods = {.5}; o.types = {RealZxy};
  result = prepare({&s1, &s2}, o);
  check(result.observations.size() == 1 && result.frequencies == std::vector<double>{2.}, "Disabled station / selection ignored");
  s1.set_data(2., {RealZxy}, {std::numeric_limits<double>::quiet_NaN()}, {.1});
  o.types = {RealZxy, RealTzx};
  result = prepare({&s1}, o);
  check(result.missing == 1 && result.observations.size() == 1, "EDI missing sentinel not omitted");
  s1.set_data(2., {RealZxy}, {std::numeric_limits<double>::infinity()}, {.1});
  rejects([&] {prepare({&s1}, o);}, "Non-finite value");
  s1.set_data(2., {RealZxy}, {0.}, {0.});
  rejects([&] {prepare({&s1}, o);}, "Standard error");
  s1.set_data(2., {RealZxy, ImagZxy}, {0., 0.}, {.1, .1});
  o.types = {RealZxy, ImagZxy, PhsZxy};
  result = prepare({&s1}, o);
  check(result.undefined_phase == 1 && result.observations.size() == 2, "Undefined phase or valid zero handling");
  s1.set_data(2., {RealZxy}, {1.}, {.1});
  result = prepare({&s1}, o);
  near(find(result, PhsZxy).value, 0.);
  o.types = {PhsZxy}; s1.set_data_mask(RealZxy, 2., false);
  rejects([&] {prepare({&s1}, o);}, "Empty export");
}

void halfspace_and_selection_validation()
{
  auto s = station();
  s.set_data(1., {RealZxx, ImagZxx, RealZyy, ImagZyy}, {0.,0.,0.,0.}, {.0001,.0001,.0001,.0001});
  s.calculate_apparent_resistivity(); s.calculate_phase(); s.calculate_phase_tensor(); s.set_error_floor(0.);
  auto o = options();
  o.types = {RealZxx, ImagZxx, PhsZxx, PhsZxy, PhsZyx, PTxx, PTxy, PTyx, PTyy};
  auto result = prepare({&s}, o);
  check(result.undefined_phase == 1, "Halfspace diagonal phase must be omitted");
  near(find(result, RealZxx).value, 0.); near(find(result, ImagZxx).value, 0.);
  near(find(result, PhsZxy).value, 45.); near(find(result, PhsZyx).value, -135.);
  near(find(result, PTxx).value, 1.); near(find(result, PTyy).value, 1.);
  near(find(result, PTxy).value, 0.); near(find(result, PTyx).value, 0.);
  o.types = {RealZxy, RealZxy};
  rejects([&] {prepare({&s}, o);}, "Duplicate selected data type");
  o.types = {log10RhoZxy};
  rejects([&] {prepare({&s}, o);}, "Log-resistivity");
  o.types = {RealZxy}; o.periods = {1.,1.};
  rejects([&] {prepare({&s}, o);}, "Duplicate selected period");
  o.periods = {-1.};
  rejects([&] {prepare({&s}, o);}, "Selected periods");
  o.periods = {1.};
  rejects([&] {prepare({&s, &s}, o);}, "Duplicate observation");
  auto other = station("S02", {1. + 5e-13});
  o.periods.push_back(1. / other.frequencies()[0]);
  rejects([&] {prepare({&s, &other}, o);}, "Near-duplicate frequencies");
}

void resistivity_phase_contract()
{
  const auto rows = parse(
    "# frequency_Hz receiver observable component value standard_error\n"
    "1.0 S01 apparent_resistivity_xy value 100 5\n"
    "1.0 S01 apparent_resistivity_yx value 100 5\n"
    "1.0 S01 phase_xy value 45 2\n"
    "1.0 S01 phase_yx value -135 2\n");
  const auto exported = parse(encode(rows, receivers, "Contract example").data_text);
  check(exported.size() == 4, "Resistivity/phase export added companion impedance rows");
  for(unsigned i = 0; i < rows.size(); ++i)
    check(exported[i].type == rows[i].type && exported[i].value == rows[i].value &&
          exported[i].error == rows[i].error, "Linear values, degree phases or absolute errors changed");

  auto s = station();
  auto o = options();
  for(const auto &selection: std::vector<std::vector<RealDataType>>{
        {RhoZxy, RhoZyx}, {PhsZxy, PhsZyx}, {RhoZxy, RhoZyx, PhsZxy, PhsZyx}}) {
    o.types = selection;
    const auto result = prepare({&s}, o);
    check(result.observations.size() == selection.size(), "Standalone selection requires impedance");
    for(const auto &row: result.observations)
      check(std::find(selection.begin(), selection.end(), row.type) != selection.end(),
            "Unexpected companion row");
  }
  auto result = prepare({&s}, o);
  const double sigma_rho = 2. * 1.0132118364233778 * .0001 / std::hypot(.002, .002);
  near(find(result, RhoZxy).error, sigma_rho);
  s.set_data_mask(RealZxy, 1., false);
  result = prepare({&s}, o);
  check(result.observations.size() == 2 && result.masked == 2, "Derived data did not respect impedance masks");
  for(const auto &row: result.observations)
    check(row.type == RhoZyx || row.type == PhsZyx, "Masked xy observation exported");
}

void response_import()
{
  auto rows = example();
  std::reverse(rows.begin(), rows.end());
  std::istringstream all(encode(rows, receivers, "Responses").data_text);
  auto stations = read_responses(all);
  check(stations.size() == 1, "Response receiver grouping");
  auto &s = stations.at("S01");
  std::vector<dvector> values, errors;
  for(const auto &row: rows) {
    const std::string observable = mapping(row.type).observable;
    const auto c = type_to_column_table.at(row.type);
    if(observable.find("apparent_resistivity_") == 0) s.get_apparent_resistivity(values, errors);
    else if(observable.find("phase_tensor_") == 0) s.get_phase_tensor(values, errors);
    else if(observable.find("phase_") == 0) s.get_phase(values, errors);
    else continue;
    near(values[c][0], row.value);
    near(errors[c][0], row.error);
  }

  std::istringstream sparse(
    " # frequency_hz receiver observable component value error\r\n"
    "100 S02 phase_yx value -135 2\r\n"
    "2 S01 induction_x real 0 0.01\r\n"
    "1 S01 impedance_xy imag 0.002 0.0001\r\n"
    "4 S01 impedance_xy real 0.004 0.0001\r\n"
    "1 S01 phase_tensor_xy value 0.75 0.1\r\n"
    "1 S01 impedance_xy real 0.002 0.0001\r\n");
  stations = read_responses(sparse);
  check(stations.size() == 2 && stations.at("S01").frequencies() == dvector({1., 2., 4.}),
        "Sparse, shuffled response keys were not preserved");
  stations.at("S01").get_apparent_resistivity(values, errors);
  near(values[1][0], 1.0132118364233778);
  check(std::isnan(values[0][0]) && std::isnan(values[1][1]) && std::isnan(values[1][2]),
        "Missing impedance components were filled with zero");
  stations.at("S01").get_phase(values, errors);
  near(values[1][0], 45.);
  check(std::isnan(values[1][2]), "Incomplete impedance yielded a phase");
  stations.at("S01").get_phase_tensor(values, errors);
  near(values[1][0], .75);
  check(std::isnan(values[0][0]), "Incomplete tensor yielded a phase tensor");
  stations.at("S01").get_tipper(values, errors);
  near(values[0][1], 0.);
  check(std::isnan(values[2][1]) && std::isnan(values[1][1]), "Missing tipper scalar became a zero");
  stations.at("S02").get_phase(values, errors);
  near(values[2][0], -135.);

  std::istringstream close_frequencies(
    "1e-8 S01 apparent_resistivity_xy value 10 1\n"
    "1.005e-8 S01 apparent_resistivity_xy value 20 1\n");
  auto close = read_responses(close_frequencies);
  close.at("S01").get_apparent_resistivity(values, errors);
  near(values[1][0], 10.); near(values[1][1], 20.);

  auto dense = example();
  dense.resize(12); // Impedance and tipper only; the other curves must be derived.
  std::istringstream derived(encode(dense, receivers, "Derived").data_text);
  auto calculated = read_responses(derived);
  calculated.at("S01").get_phase_tensor(values, errors);
  for(const auto &component: values)
    check(std::isfinite(component[0]), "Full impedance did not produce a phase tensor");
}

void response_rms()
{
  auto observed = station("S01", {1., 2., 4.});
  std::istringstream partial(
    "2 S01 impedance_xy imag 0.0022 0.0001\n"
    "2 S01 impedance_xy real 0.0021 0.0001\n"
    "8 S01 impedance_xy real 9 0.0001\n");
  auto responses = read_responses(partial);
  near(observed.rms(responses.at("S01")), std::sqrt(2.5));
  observed.set_data_mask(RealZxy, 2., false);
  check(std::isnan(observed.rms(responses.at("S01"))), "RMS included masked/unmatched observations");

  std::istringstream rounded("1.0000001 S01 impedance_xy real 0.0021 0.0001\n");
  responses = read_responses(rounded);
  near(observed.rms(responses.at("S01")), 1.);

  std::istringstream phase_only("1 S01 phase_xy value 405 2\n");
  responses = read_responses(phase_only);
  near(observed.rms(responses.at("S01")), 0.);
  observed.set_active(false);
  check(std::isnan(observed.rms(responses.at("S01"))), "RMS included disabled station");
}

void station_initialization()
{
  auto s = station();
  double value, error;
  rejects([&] { s.set_data(1., {RealZxy, ImagZxy}, {.3}, {.01, .01}); }, "same length");
  check(s.scalar_value(RealZxy, 0, value, error), "Rejected update changed station"); near(value, .002);
  s.set_active(false); s.set_data_mask(RealTzx, 1., false);
  s.set_size(3, true); s.set_frequencies({1., 2., 3.});
  check(s.active() && s.tipper_mask()[0].size() == 3 && s.tipper_mask()[0][2], "Reinitialization left stale masks or sizes");
  check(!s.scalar_value(RealZxy, 0, value, error), "Reinitialization retained old values");
  s.set_data(3., {RealZxy, ImagZxy}, {.3, .4}, {.01, .01});
  check(s.scalar_value(ImagZxy, 2, value, error), "Reinitialized rows were not resized"); near(value, .4);
}
}

int main(int argc, char **argv)
{
  try {
    format_roundtrip(); invalid_format(); precision_and_comments(); coordinates();
    frequencies_and_distortion(); station_export(); masks_and_missing(); halfspace_and_selection_validation();
    resistivity_phase_contract();
    response_import(); response_rms(); station_initialization();
    if(argc == 2) {
      const auto output = encode(example(), receivers, "Specification example in local Cartesian model coordinates");
      std::ofstream(std::string(argv[1]) + ".data") << output.data_text;
      std::ofstream(std::string(argv[1]) + ".recvs") << output.receiver_text;
      std::ofstream(std::string(argv[1]) + ".freqs") << output.frequency_text;
    }
    std::cout << "Native MT format, validation, conventions, masks and round-trip checks passed.\n";
  } catch(const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}
