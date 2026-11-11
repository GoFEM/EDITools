#include "include/ImpedanceAnalysis.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

using namespace ImpedanceAnalysis;
namespace {
const double pi = std::acos(-1.), mu = 4e-7*pi;
using Complex = std::complex<double>;
void check(bool good, const char *text) { if(!good) throw std::runtime_error(text); }
void near(double a, double b, double tolerance = 1e-8) { check(std::isfinite(a) && std::abs(a-b) <= tolerance*std::max(1., std::abs(b)), "Numerical mismatch"); }
void propagate(Sample &s) {
  const double relative=s.error/std::abs(s.impedance);
  s.rhoError=2.*s.resistivity()*relative; s.phaseError=relative*180./pi;
}
Sample halfspace(double period, double rho) {
  Sample s; s.period = period; s.impedance = std::sqrt(Complex(0., mu*2*pi/period*rho)); s.error = .02*std::abs(s.impedance); propagate(s); return s;
}
MTStationData station(const std::array<Complex, 4> &z) {
  MTStationData s; s.set_size(1, true); s.set_frequencies({1.});
  s.set_data(1., {RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy, ImagZyy},
    {z[0].real(), z[0].imag(), z[1].real(), z[1].imag(), z[2].real(), z[2].imag(), z[3].real(), z[3].imag()}, std::vector<double>(8, .01));
  return s;
}
void invariants() {
  const auto z = halfspace(1., 100.).impedance;
  auto s = station({{0., z, -z, 0.}});
  for(auto quantity: {Quantity::XY, Quantity::YX, Quantity::SSQ, Quantity::Berdichevsky, Quantity::Determinant}) {
    const auto values = curve(s, quantity); near(values[0].resistivity(), 100.); near(values[0].phase(), 45.);
  }
  near(curve(s, Quantity::SSQ)[0].error, .01/std::sqrt(2.));
  s.set_data_mask(RealZxx, 1., false);
  check(!curve(s, Quantity::SSQ)[0].usable(), "SSQ used masked tensor entry");
  check(!curve(s, Quantity::Determinant)[0].usable(), "Determinant used masked entry");
  check(curve(s, Quantity::Berdichevsky)[0].usable(), "Off-diagonal average requires diagonals");
  s.set_active(false); check(!curve(s, Quantity::XY)[0].usable(), "Disabled station used");
  // Rotate a nonsymmetric full tensor; SSQ uses complex squares, not magnitudes.
  const std::array<Complex,4> tensor{{{.2,.1},{2.,3.},{-4.,2.},{.7,-.2}}};
  const double c = std::cos(.71), t = std::sin(.71);
  std::array<Complex,4> rotated{{c*c*tensor[0]+c*t*(tensor[1]+tensor[2])+t*t*tensor[3],
    -c*t*tensor[0]+c*c*tensor[1]-t*t*tensor[2]+c*t*tensor[3],
    -c*t*tensor[0]-t*t*tensor[1]+c*c*tensor[2]+c*t*tensor[3],
    t*t*tensor[0]-c*t*(tensor[1]+tensor[2])+c*c*tensor[3]}};
  for(auto q: {Quantity::SSQ, Quantity::Berdichevsky, Quantity::Determinant}) {
    auto a=curve(station(tensor),q)[0], b=curve(station(rotated),q)[0];
    near(a.resistivity(),b.resistivity()); near(a.phase(),b.phase());
  }
}
void rotations() {
  const std::array<Complex,4> z{{{.2,.1},{2.,3.},{-4.,2.},{.7,-.2}}};
  auto original = station(z);
  original.calculate_apparent_resistivity(); original.calculate_phase(); original.set_error_floor(.05);
  const auto before = sample(original, Quantity::XY, 0);
  const double radians = 30. * pi / 180., c = std::cos(radians), s = std::sin(radians);
  const Complex expected = -c*s*z[0] + c*c*z[1] - s*s*z[2] + c*s*z[3];
  auto rotated = sample(original, Quantity::XY, 0, 30.);
  near(rotated.impedance.real(), expected.real()); near(rotated.impedance.imag(), expected.imag());
  double variance = 0.; const double weights[]{-c*s,c*c,-s*s,c*s};
  const RealDataType types[]{RealZxx,RealZxy,RealZyx,RealZyy};
  for(unsigned k=0;k<4;++k){double value,error;original.scalar_value(types[k],0,value,error);variance+=weights[k]*weights[k]*error*error;}
  near(rotated.error,std::sqrt(variance));
  near(rotated.rhoError,2.*rotated.resistivity()*rotated.error/std::abs(rotated.impedance));
  near(rotated.phaseError,180./pi*rotated.error/std::abs(rotated.impedance));
  near(sample(original,Quantity::XY,0,90.).impedance.real(),-z[2].real());
  near(sample(original,Quantity::YX,0,90.).impedance.imag(),z[1].imag());
  for(auto q:{Quantity::SSQ,Quantity::Berdichevsky,Quantity::Determinant}){
    const auto a=sample(original,q,0), b=sample(original,q,0,37.);
    near(a.resistivity(),b.resistivity());near(a.phase(),b.phase());near(a.error,b.error);
  }
  near(sample(original,Quantity::XY,0,180.).rhoError,before.rhoError);
  near(sample(original,Quantity::XY,0).error,before.error);
  original.set_data_mask(RealZxx,1.,false);
  check(!sample(original,Quantity::XY,0,30.).usable(),"Rotation uses masked contributors");
  check(sample(original,Quantity::XY,0,90.).usable(),"Quarter-turn requires an unrelated masked diagonal");
  check(sample(original,Quantity::Berdichevsky,0,30.).usable(),"Invariant needs unnecessary rotation contributors");
  check(!sample(original,Quantity::XY,0,NAN).usable(),"Unavailable rotation silently uses original axes");
}
void averaging() {
  auto a=halfspace(1.,10.), b=halfspace(1.,1000.);
  auto mean=regionalMean({{a},{b}}); near(mean[0].resistivity,100.); near(mean[0].phase,45.); check(mean[0].stations==2,"Wrong count");
  check(mean[0].logRhoSpread>0.,"Spread missing");
  b.period=1.0005;
  check(regionalMean({{a},{b}}).size()==2,"Exact periods merged by default");
  check(regionalMean({{a},{b}},.001)[0].stations==2,"Tolerance did not group nearby periods");
  auto c=b; c.period=1.0011; check(regionalMean({{a},{b},{c}},.001).size()==2,"Tolerance chained beyond limit");
  auto duplicate=a; duplicate.period=1.0001;
  check(regionalMean({{a,duplicate},{b}},.001)[0].stations==1,"A station was counted twice");
  check(std::isnan(regionalMean({{a},{b}},0.,2)[0].resistivity),"Minimum count ignored");
  auto missing=a; missing.impedance={NAN,NAN};
  check(regionalMean({{missing}})[0].stations==0,"Masked period lost");
  a.impedance=std::polar(std::abs(a.impedance),179*pi/180); b=a; b.impedance=std::conj(a.impedance);
  near(std::abs(regionalMean({{a},{b}})[0].phase),180.);
}
void rho() {
  std::vector<Sample> samples;
  for(unsigned i=0;i<24;++i) samples.push_back(halfspace(std::pow(10.,-3.+i*7./23),100.));
  auto result=rhoPlus(samples);
  std::cout << "Halfspace RMS " << result.normalizedRms << ", chi " << result.chiSquare << ", delta " << result.refinementChange << ", converged " << result.converged << std::endl;
  check(result.converged && result.gridStable,"Halfspace solver failed");
  check(result.normalizedRms<2e-4,"Halfspace rho+ fit inaccurate");
  for(const auto &s:result.response) { near(s.resistivity(),100.,1e-4); near(s.phase(),45.,1e-4); }
  // Exact two-layer response from the independent impedance recursion.
  for(auto &s:samples) {
    const double w=2*pi/s.period;
    const auto top=std::sqrt(Complex(0.,w*mu*10.)), bottom=std::sqrt(Complex(0.,w*mu*1000.));
    const auto factor=std::tanh(std::sqrt(Complex(0.,w*mu/10.))*1500.);
    s.impedance=top*(bottom+top*factor)/(top+bottom*factor); s.error=.03*std::abs(s.impedance); propagate(s);
  }
  result=rhoPlus(samples);
  std::cout << "Layered fit: RMS " << result.normalizedRms << ", converged " << result.converged << ", refinement " << result.refinementChange << '\n';
  check(result.converged && result.gridStable && result.normalizedRms<.05,"Layered rho+ fit failed");
  const auto original=samples;
  samples[5].impedance={NAN,NAN}; samples[7].error=0.; propagate(samples[7]);
  result=rhoPlus(samples); check(result.excluded==1 && result.invalidErrors==1,"Missing values/errors not reported");
  samples=original;
  // Positive-quadrant phases alone are insufficient: amplitude/phase must agree.
  samples[10].impedance *= 100.; propagate(samples[10]);
  result=rhoPlus(samples); check(result.normalizedRms>1.,"Amplitude-phase inconsistency accepted");
  // Use the supplied rho and phase errors independently, without a hidden floor.
  for(auto &s:samples) { s.rhoError*=2.; s.phaseError*=3.; }
  result=rhoPlus(samples);
  double chi=0.;
  for(unsigned i=0;i<result.observed.size();++i) {
    const auto &s=result.observed[i], &p=result.fitted[i];
    near(result.rhoResidual[i],std::log(p.resistivity()/s.resistivity())/(s.rhoError/s.resistivity()));
    near(result.phaseResidual[i],(p.phase()-s.phase())/s.phaseError);
    chi+=std::pow(result.rhoResidual[i],2)+std::pow(result.phaseResidual[i],2);
  }
  near(result.chiSquare,chi);
  samples=original;
  for(unsigned i=0;i<samples.size();++i)
    samples[i].impedance *= std::polar(std::exp(.02*std::sin(i*1.7)), .02*std::cos(i*2.3));
  result=rhoPlus(samples); check(result.converged && result.normalizedRms<1.,"Noisy valid response solver failed");
  samples=original;
  for(auto &s:samples) s.impedance=std::conj(s.impedance);
  result=rhoPlus(samples); check(result.normalizedRms>5.,"Non-passive phase accepted");
  bool rejected=false; try { rhoPlus({samples[0],samples[1]}); } catch(const std::invalid_argument &) { rejected=true; }
  check(rejected,"Underdetermined data accepted");
}
void outliers() {
  RhoPlusResult r; r.converged=true; r.gridStable=true;
  r.observed.resize(6);
  r.rhoResidual={3.,-3.01,0.,5.,NAN,0.};
  r.phaseResidual={-3.,0.,3.01,-7.,0.,NAN};
  check(rhoPlusOutliers(r,3.)==std::vector<unsigned>({1,2,3}),"Residual threshold, sign or union is incorrect");
  check(rhoPlusOutliers(r,2.)==std::vector<unsigned>({0,1,2,3}),"Adjustable threshold ignored");
  r.gridStable=false; check(rhoPlusOutliers(r,2.).empty(),"Unstable fit supplied candidates");
  r.gridStable=true; r.converged=false; check(rhoPlusOutliers(r,2.).empty(),"Unconverged fit supplied candidates");
  bool rejected=false; try { rhoPlusOutliers(r,0.); } catch(const std::invalid_argument &) { rejected=true; }
  check(rejected,"Invalid threshold accepted");
}
}
int main() {
  try { invariants(); rotations(); averaging(); rho(); outliers(); std::cout<<"Impedance analysis checks passed.\n"; }
  catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
