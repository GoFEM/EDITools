#include "include/MTSpectra.h"
#include "include/EDIFileReader.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
using Complex = std::complex<double>;
using Matrix = std::array<std::array<Complex, 7>, 7>;
const std::array<Complex, 4> impedance{{{1., .5}, {2., 3.}, {-4., 2.}, {.7, -.2}}};
const std::array<Complex, 2> tipper{{{.2, .05}, {-.1, .07}}};
void check(bool value, const char *message) { if(!value) throw std::runtime_error(message); }
void near(double a, double b, double tolerance = 1e-12) {
  check(std::isfinite(a) && std::abs(a - b) <= tolerance * std::max(1., std::abs(b)), "Spectral numerical mismatch");
}
void near(Complex a, Complex b) { near(a.real(), b.real()); near(a.imag(), b.imag()); }
void rejects(const std::function<void()> &run) {
  try { run(); } catch(const std::exception &) { return; }
  throw std::runtime_error("Invalid spectral input was accepted");
}
MTSpectra::PackedMatrix synthetic() {
  // Independent unit-power sources. Remote channels contain independent noise,
  // so the exact transfer functions and residual variances are known separately.
  Matrix mixing{};
  mixing[0][0] = 2.; mixing[1][1] = 3.;
  for(unsigned c = 0; c < 2; ++c) {
    mixing[2][c] = tipper[c] * mixing[c][c];
    mixing[3][c] = impedance[c] * mixing[c][c];
    mixing[4][c] = impedance[c + 2] * mixing[c][c];
  }
  mixing[2][2] = .5; mixing[3][3] = .6; mixing[4][4] = .8;
  mixing[5][0] = 2.; mixing[5][5] = .4;
  mixing[6][1] = 3.; mixing[6][6] = .7;
  Matrix power{};
  for(unsigned i = 0; i < 7; ++i) for(unsigned j = 0; j < 7; ++j)
    for(unsigned k = 0; k < 7; ++k) power[i][j] += mixing[i][k] * std::conj(mixing[j][k]);
  MTSpectra::PackedMatrix packed{};
  for(unsigned i = 0; i < 7; ++i) for(unsigned j = 0; j < 7; ++j)
    packed[i][j] = i < j ? -power[i][j].imag() : power[i][j].real();
  return packed;
}
void analytical() {
  const auto packed = synthetic();
  for(double scale: {1e-200, 1., 1e200}) {
    auto scaled = packed; for(auto &row: scaled) for(auto &value: row) value *= scale;
    const auto result = MTSpectra::estimate(scaled, 200.);
    for(unsigned c = 0; c < 4; ++c) {
      near(result.impedance[c], impedance[c]);
      const double weight = c % 2 ? 9.49 / 81. : 4.16 / 16.;
      near(result.impedanceError[c], std::sqrt((c < 2 ? .36 : .64) * weight / 200.));
    }
    for(unsigned c = 0; c < 2; ++c) {
      near(result.tipper[c], tipper[c]);
      near(result.tipperError[c], std::sqrt(.25 * (c ? 9.49 / 81. : 4.16 / 16.) / 200.));
    }
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto missing = packed;
  for(unsigned i = 0; i < 7; ++i) missing[2][i] = missing[i][2] = nan;
  auto result = MTSpectra::estimate(missing, 200.);
  for(unsigned c = 0; c < 4; ++c) near(result.impedance[c], impedance[c]);
  check(std::isnan(result.tipper[0].real()), "Missing Hz became zero tipper");
  missing = packed;
  for(unsigned i = 0; i < 7; ++i) missing[3][i] = missing[i][3] = nan;
  result = MTSpectra::estimate(missing, 200.);
  check(std::isnan(result.impedance[0].real()), "Missing Ex became zero impedance");
  near(result.impedance[2], impedance[2]); near(result.tipper[0], tipper[0]);
  for(double averages: {0., -1., nan}) rejects([&] { MTSpectra::estimate(packed, averages); });
  auto singular = packed;
  for(unsigned remote: {5u, 6u}) singular[remote][1] = singular[1][remote] = 0.;
  rejects([&] { MTSpectra::estimate(singular, 200.); });
  singular = packed; singular[5][5] = 0.;
  rejects([&] { MTSpectra::estimate(singular, 200.); });
}
void reference() {
  std::ifstream file(SPECTRA_TEST_DATA "/mt_spectra_reference.txt");
  check(bool(file), "Missing spectral regression cases");
  unsigned cases = 0;
  for(std::string line; std::getline(file, line);) {
    if(line.empty() || line[0] == '#') continue;
    std::istringstream row(line); double averages; row >> averages;
    MTSpectra::PackedMatrix packed;
    for(auto &values: packed) for(auto &value: values) row >> value;
    const auto result = MTSpectra::estimate(packed, averages);
    std::vector<double> actual;
    for(auto z: result.impedance) { actual.push_back(z.real()); actual.push_back(z.imag()); }
    for(auto z: result.tipper) { actual.push_back(z.real()); actual.push_back(z.imag()); }
    actual.insert(actual.end(), result.impedanceError.begin(), result.impedanceError.end());
    actual.insert(actual.end(), result.tipperError.begin(), result.tipperError.end());
    for(double value: actual) { double expected; check(bool(row >> expected), "Incomplete spectral reference case"); near(value, expected, 2e-5); }
    ++cases;
  }
  check(cases == 6, "Missing spectral reference cases");
}
std::string edi(const MTSpectra::PackedMatrix &packed, unsigned channels = 7, unsigned values = 49, double frequency = 2.) {
  std::ostringstream text; text << std::setprecision(17);
  text << ">HEAD\nDATAID=SPECTRA\nLAT=55\nLONG=9\nELEV=100\nEMPTY=1e32\n>=SPECTRASECT\nNCHAN=" << channels
       << "\n// " << channels << '\n';
  for(unsigned c = 0; c < channels; ++c) text << c + 1 << ' ';
  text << "\n>SPECTRA FREQ= " << frequency << " AVGT= 200 // 49\n";
  for(unsigned n = 0; n < values; ++n) text << packed[n / 7][n % 7] << (n % 7 == 6 ? '\n' : ' ');
  text << "\n>END\n";
  return text.str();
}
void reader() {
  const std::string path = SPECTRA_TEST_OUTPUT "/native-spectra-test.edi";
  auto load = [&](const std::string &text) {
    { std::ofstream file(path); file << text; check(bool(file), "Cannot write EDI fixture"); }
    return EDIFileReader(path).get_mt_data();
  };
  const auto packed = synthetic(); const auto text = edi(packed);
  const auto station = load(text);
  const double factor = 4. * std::acos(-1.) * 1e-4;
  const auto expected = MTSpectra::estimate(packed, 200.);
  const std::array<RealDataType, 4> real{{RealZxx, RealZxy, RealZyx, RealZyy}}, imag{{ImagZxx, ImagZxy, ImagZyx, ImagZyy}};
  for(unsigned c = 0; c < 4; ++c) {
    double value, error;
    check(station.scalar_value(real[c], 0, value, error), "Spectral impedance unavailable");
    near(value, factor * impedance[c].real()); near(error, factor * expected.impedanceError[c]);
    check(station.scalar_value(imag[c], 0, value, error), "Spectral imaginary impedance unavailable");
    near(value, factor * impedance[c].imag());
  }
  double value, error;
  check(station.scalar_value(RealTzx, 0, value, error), "Spectral tipper unavailable");
  near(value, tipper[0].real()); near(error, expected.tipperError[0]);
  check(station.scalar_value(RhoZxy, 0, value, error), "Derived resistivity unavailable");
  near(value, .2 * std::norm(impedance[1]) / 2.);
  check(station.scalar_value(PhsZxy, 0, value, error), "Derived phase unavailable");
  near(value, std::arg(impedance[1]) * 180. / std::acos(-1.));
  check(station.scalar_value(PTxx, 0, value, error), "Derived phase tensor unavailable");
  // Without EMPTY, zero cross-powers remain valid and the default sentinel applies.
  auto noEmpty = text;
  noEmpty.erase(noEmpty.find("EMPTY=1e32\n"), std::string("EMPTY=1e32\n").size());
  const auto defaultEmpty = load(noEmpty);
  check(defaultEmpty.scalar_value(RealZxy, 0, value, error), "Zero cross-power treated as missing");
  near(value, factor * impedance[1].real());
  check(defaultEmpty.scalar_value(RealTzx, 0, value, error), "Tipper missing without EMPTY header");
  near(value, tipper[0].real());
  // Standard scientific exponent spelling is accepted in spectral values.
  auto exponent = text; const auto pos = exponent.find('\n', exponent.find(">SPECTRA FREQ")) + 1;
  exponent.replace(pos, 1, "4D+00"); load(exponent);
  rejects([&] { load(edi(packed, 5)); });
  rejects([&] { load(edi(packed, 7, 47)); });
  rejects([&] { load(edi(packed, 7, 49, 0.)); });
}
}
int main() {
  try { analytical(); reference(); reader(); std::cout << "Native spectra, reference values and EDI integration passed.\n"; }
  catch(const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
