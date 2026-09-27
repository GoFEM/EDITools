#ifndef EDITOOLS_IMPEDANCE_ANALYSIS_H
#define EDITOOLS_IMPEDANCE_ANALYSIS_H

#include "MTStationData.h"
#include <cmath>

namespace ImpedanceAnalysis {
enum class Quantity { XX, XY, YX, YY, SSQ, Berdichevsky, Determinant };
const char *name(Quantity quantity);
struct Sample {
  double period = 0.;
  std::complex<double> impedance{NAN, NAN}; // SI E/H; YX is sign-reversed.
  double error = NAN; // One-standard-deviation error in either complex part.
  double rhoError = NAN, phaseError = NAN; // Main-data propagated errors: Ω m and degrees.
  bool usable() const;
  double resistivity() const;
  double phase() const; // Degrees, without folding anomalous phases.
};
// Keeps missing/masked periods as NaNs so plots do not bridge missing samples.
std::vector<Sample> curve(const MTStationData &station, Quantity quantity);

struct RegionalPoint {
  double period = 0., resistivity = NAN, phase = NAN;
  double logRhoSpread = NAN; // Sample standard deviation of ln(rho), not error of the mean.
  unsigned stations = 0;
};
// Equal station weights. Bounded period groups; ambiguous repeated station
// samples are excluded from that group. No interpolation or extrapolation.
std::vector<RegionalPoint> regionalMean(const std::vector<std::vector<Sample>> &curves,
                                      double relativeTolerance = 0., unsigned minimumStations = 1);

struct RhoPlusOptions {
  unsigned intervalsPerDecade = 12;
};
struct RhoPlusResult {
  std::vector<Sample> observed, fitted, response;
  std::vector<double> rhoResidual, phaseResidual;
  double chiSquare = NAN, normalizedRms = NAN, refinementChange = NAN;
  unsigned excluded = 0, invalidErrors = 0;
  bool converged = false, gridStable = false;
};
// Finite-grid rho+ approximation: bounded spectral density in log admittance,
// exact interval integrals, a free real offset, and both spectral tails.
// Fits ln(rho) and phase jointly, then repeats with twice the grid density.
// Diagnostic of 1-D consistency, not proof of general 3-D data validity.
RhoPlusResult rhoPlus(const std::vector<Sample> &samples, const RhoPlusOptions &options = {});
// Indices into observed. Either normalized residual must strictly exceed the
// threshold; provisional fits cannot supply masking candidates.
std::vector<unsigned> rhoPlusOutliers(const RhoPlusResult &result, double threshold);
}
#endif
