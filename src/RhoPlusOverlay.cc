#include "include/RhoPlusOverlay.h"
#include "include/MTDataPlot.h"
#include <cmath>

namespace {
bool same(double a,double b) { return a==b || (std::isnan(a) && std::isnan(b)); }
bool same(const std::vector<ImpedanceAnalysis::Sample> &a,const std::vector<ImpedanceAnalysis::Sample> &b) {
  if(a.size()!=b.size()) return false;
  for(unsigned i=0;i<a.size();++i)
    if(!same(a[i].period,b[i].period) || !same(a[i].error,b[i].error) ||
       !same(a[i].rhoError,b[i].rhoError) || !same(a[i].phaseError,b[i].phaseError) ||
       !same(a[i].impedance.real(),b[i].impedance.real()) || !same(a[i].impedance.imag(),b[i].impedance.imag())) return false;
  return true;
}
const std::array<QString,4> names{{"XX","XY","YX","YY"}};
bool supported(unsigned component) { return component==1 || component==2; }
}
RhoPlusOverlay::RhoPlusOverlay(QCustomPlot *rho, QCustomPlot *phase, std::function<void(const QString &)> status)
  :rho(rho),phase(phase),report(std::move(status)) {
  for(unsigned c=0;c<4;++c) for(auto *plot:{rho,phase}) {
    auto *line=new QCPCurve(plot->xAxis,plot->yAxis);
    line->setObjectName("rhoPlus"+names[c]); line->setName(names[c]+" ρ+");
    line->setPen(QPen(PlotColors::componentColors()[c],.9,Qt::DashLine));
    line->setSelectable(QCP::stNone); line->setVisible(false); line->removeFromLegend();
    (plot==rho?rhoLines:phaseLines)[c]=line;
  }
  timer.setSingleShot(true); timer.setInterval(120);
  QObject::connect(&timer,&QTimer::timeout,&timer,[this]{calculate();});
}
void RhoPlusOverlay::hide() {
  for(unsigned c=0;c<4;++c) { rhoLines[c]->setVisible(false); phaseLines[c]->setVisible(false); }
}
void RhoPlusOverlay::update(const MTStationData *station, bool on, bool wrap,
                           const std::array<bool,4> &rVisible,const std::array<bool,4> &pVisible,
                           const ImpedanceAnalysis::RhoPlusOptions &settings) {
  enabled=on; phaseWrap=wrap; rhoVisible=rVisible; phaseVisible=pVisible;
  const bool changed=settings.intervalsPerDecade!=options.intervalsPerDecade;
  options=settings;
  if(changed) for(auto &fit:components) { fit.response.clear(); fit.pending=true; }
  if(!station || !enabled) {
    timer.stop(); hide(); rho->replot(); phase->replot();
    report(enabled?QObject::tr("Select a station."):QString()); return;
  }
  bool pending=false;
  for(unsigned c=0;c<4;++c) if(supported(c)) {
    auto input=ImpedanceAnalysis::curve(*station,ImpedanceAnalysis::Quantity(c));
    auto &fit=components[c];
    if(changed || !same(input,fit.input)) { fit.input=std::move(input); fit.response.clear(); fit.pending=true; fit.status=QObject::tr("updating"); }
    pending |= fit.pending && (rhoVisible[c] || phaseVisible[c]);
  }
  draw();
  if(pending) timer.start();
}
void RhoPlusOverlay::calculate() {
  if(!enabled) return;
  for(unsigned c=0;c<4;++c) {
    auto &fit=components[c];
    if(!supported(c) || !fit.pending || !(rhoVisible[c] || phaseVisible[c])) continue;
    fit.pending=false;
    try {
      const auto result=ImpedanceAnalysis::rhoPlus(fit.input,options);
      if(result.converged && result.gridStable) {
        fit.response=result.response;
        fit.status=QObject::tr("nRMS %1").arg(result.normalizedRms,0,'g',3);
      } else fit.status=QObject::tr("fit omitted; check convergence in the validity window");
    } catch(const std::exception &e) { fit.status=QString::fromUtf8(e.what()); }
  }
  draw();
}
void RhoPlusOverlay::draw() {
  QStringList status;
  for(unsigned c=0;c<4;++c) {
    const auto &fit=components[c]; QVector<double> periods,rhoValues,phaseValues;
    for(const auto &sample:fit.response) {
      periods.push_back(sample.period); rhoValues.push_back(sample.resistivity());
      // Fits use -Zyx; restore its native phase when phase wrapping is off.
      phaseValues.push_back(sample.phase() - (c==2 && !phaseWrap ? 180. : 0.));
    }
    rhoLines[c]->setData(periods,rhoValues); phaseLines[c]->setData(periods,phaseValues);
    const bool available=enabled && supported(c) && !fit.response.empty();
    rhoLines[c]->setVisible(available && rhoVisible[c]); phaseLines[c]->setVisible(available && phaseVisible[c]);
    if(enabled && supported(c)) status.push_back(names[c]+": "+fit.status);
  }
  // Updating the guides never rescales the observations or the user's view.
  rho->replot(); phase->replot(); report(status.join('\n'));
}
