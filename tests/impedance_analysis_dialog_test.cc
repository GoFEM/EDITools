#include "ImpedanceAnalysisWindows.h"
#include "include/MTSurveyData.h"
#include "include/MTDataPlot.h"
#include "qcustomplot.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QFileDialog>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QProgressDialog>
#include <QSpinBox>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTextStream>
#include <QThread>
#include <QTimer>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool value,const char *text) { if(!value) throw std::runtime_error(text); }
template<class T> T *control(QWidget &w,const char *name) { auto *v=w.findChild<T *>(name); check(v,"Missing analysis control"); return v; }
void settle() {
  QElapsedTimer timer; timer.start();
  while(timer.elapsed()<250) { QApplication::processEvents(QEventLoop::AllEvents,25); QThread::msleep(5); }
}
void masking(const std::shared_ptr<MTSurveyData> &original, const QString &previewDirectory) {
  auto survey=std::make_shared<MTSurveyData>(*original);
  auto &s=survey->get_station_data("A");
  const auto frequency=s.frequencies().at(5);
  double real,error; check(s.scalar_value(RealZxy,5,real,error),"Missing synthetic impedance");
  // A single non-passive phase within an otherwise exact half-space response.
  s.set_data(frequency,{RealZxy,ImagZxy},{real,-real},{error,error});
  for(auto f:s.frequencies()) s.set_data(f,{RealTzx,ImagTzx,RealTzy,ImagTzy},{.1,.2,.3,.4},{.01,.01,.01,.01});
  s.set_data_mask(PTxx,frequency,false);
  s.set_data_mask(RealZxx,s.frequencies().front(),false);
  s.set_data_mask(RealZyy,s.frequencies().front(),false);
  const auto z=s.impedance_mask(),pt=s.phase_tensor_mask(),tipper=s.tipper_mask();
  unsigned notifications=0;
  ImpedanceValidityWindow *window=nullptr;
  ImpedanceValidityWindow validity(nullptr,{},[&]{ ++notifications; window->setData(survey); }); window=&validity;
  validity.setData(survey,"A"); validity.show(); settle();
  auto *preview=control<QCheckBox>(validity,"validityPreviewMasks");
  auto *apply=control<QPushButton>(validity,"validityApplyMasks");
  auto *undo=control<QPushButton>(validity,"validityUndoMasks");
  auto *threshold=control<QDoubleSpinBox>(validity,"validityMaskThreshold");
  auto *quantity=control<QComboBox>(validity,"validityQuantity");
  check(threshold->value()==3.,"Default threshold should be 3 sigma");
  preview->setChecked(true); settle();
  check(apply->isEnabled(),"Outlier preview has no candidates");
  auto *graph=control<QCustomPlot>(validity,"analysisRhoPlot")->findChild<QCPGraph *>("rhoPlusMaskCandidates");
  check(graph && !graph->data()->isEmpty(),"Preview highlights missing");
  if(!previewDirectory.isEmpty()) {
    QApplication::processEvents();
    check(validity.grab().save(previewDirectory+"/mask-preview.png"),"Could not save masking preview");
  }
  check(s.impedance_mask()==z && s.phase_tensor_mask()==pt && s.tipper_mask()==tipper,"Preview changed masks");
  threshold->setValue(2.); check(!apply->isEnabled(),"Threshold left a stale preview enabled");
  threshold->setValue(3.);
  quantity->setCurrentIndex(quantity->findData(int(ImpedanceAnalysis::Quantity::SSQ)));
  check(!preview->isEnabled() && !apply->isEnabled(),"Invariant masking allowed");
  quantity->setCurrentIndex(quantity->findData(-1));
  auto *xy=control<QCheckBox>(validity,"validityXY");
  xy->setChecked(false); preview->setChecked(true); settle();
  check(!apply->isEnabled(),"Unchecked XY component supplied candidates");
  xy->setChecked(true);
  preview->setChecked(true); settle(); apply->click();
  check(!s.impedance_mask()[1][5],"Outlier was not masked");
  check(!s.impedance_mask()[0][5],"Zxy masking did not mask Zxx");
  for(unsigned c:{2u,3u}) check(s.impedance_mask()[c]==z[c],"Masking changed the opposite polarization");
  for(unsigned f=0;f<s.frequencies().size();++f)
    check(s.impedance_mask()[0][f]==(z[0][f] && s.impedance_mask()[1][f]),"Zxx masking affected a different period");
  for(unsigned c=0;c<4;++c) check(!s.phase_tensor_mask()[c][5],"Linked PT mask missing");
  check(s.tipper_mask()==tipper,"Masking changed tippers");
  check(notifications==1 && undo->isEnabled() && !apply->isEnabled(),"Apply refresh or one-pass behavior broken");
  const auto after=s.impedance_mask(); settle(); check(s.impedance_mask()==after,"Refit caused another masking pass");
  undo->click(); check(s.impedance_mask()==z && s.phase_tensor_mask()==pt && s.tipper_mask()==tipper,"Undo did not restore exact previous masks");
  check(notifications==2 && !undo->isEnabled(),"Undo notification missing");
  validity.setLinkTensorMasks(false); preview->setChecked(true); settle(); apply->click();
  check(!s.impedance_mask()[1][5] && !s.impedance_mask()[0][5] && s.phase_tensor_mask()==pt,"Polarization masking incorrectly depends on PT linking");
  undo->click();
  // Reject stale previews even if a data change was not notified to the dialog.
  preview->setChecked(true); settle(); s.set_data_mask(RealZyy,s.frequencies().back(),false);
  const auto external=s.impedance_mask(); apply->click();
  check(s.impedance_mask()==external && !apply->isEnabled(),"Stale preview applied masks");
  control<QPushButton>(validity,"validityRun")->click(); apply->click(); check(undo->isEnabled(),"Second apply failed");
  validity.setData(survey); check(!undo->isEnabled() && !apply->isEnabled(),"External refresh retained undo or preview");
  control<QComboBox>(validity,"validityStation")->setCurrentText("B"); preview->setChecked(true); settle();
  check(!apply->isEnabled(),"Clean half-space supplied candidates");
  control<QComboBox>(validity,"validityStation")->setCurrentText("DISABLED"); preview->setChecked(true); settle();
  check(!apply->isEnabled(),"Disabled station supplied candidates");
  // Shared-period XY/YX outliers must deduplicate linked PT edits and undo both.
  s=original->get_station_data("A");
  s.set_data(frequency,{RealZxy,ImagZxy,RealZyx,ImagZyx},{real,-real,-real,real},{error,error,error,error});
  const auto bothZ=s.impedance_mask(),bothPT=s.phase_tensor_mask();
  validity.setLinkTensorMasks(true); validity.setData(survey,"A"); preview->setChecked(true); settle();
  const auto graphs=control<QCustomPlot>(validity,"analysisRhoPlot")->findChildren<QCPGraph *>("rhoPlusMaskCandidates");
  check(graphs.size()==2 && graphs[0]->data()->size()==graphs[1]->data()->size(),"Combined preview lost a component");
  const int linkedCount=4*graphs[0]->data()->size();
  check(control<QLabel>(validity,"validityMaskStatus")->text().contains(QString("%1 linked PT masks").arg(linkedCount)),"Shared PT masks counted twice");
  check(control<QLabel>(validity,"validityMaskStatus")->text().contains(QString("%1 diagonal masks").arg(linkedCount/2)),"Preview omitted paired diagonal count");
  apply->click();
  for(unsigned c=0;c<4;++c) check(!s.impedance_mask()[c][5],"Combined apply missed a polarization component");
  undo->click(); check(s.impedance_mask()==bothZ && s.phase_tensor_mask()==bothPT,"Combined undo did not restore masks");
  // YX-only masking must pair with YY without touching XX/XY.
  xy->setChecked(false); settle(); apply->click();
  check(!s.impedance_mask()[2][5] && !s.impedance_mask()[3][5],"Zyx masking did not mask Zyy");
  check(s.impedance_mask()[0]==bothZ[0] && s.impedance_mask()[1]==bothZ[1],"YX-only masking changed XX/XY");
  undo->click(); check(s.impedance_mask()==bothZ && s.phase_tensor_mask()==bothPT,"YX undo failed");
  xy->setChecked(true);
  s.mask_type(RealZyx,false); validity.setData(survey); preview->setChecked(true); settle();
  check(apply->isEnabled() && control<QLabel>(validity,"validityMaskStatus")->text().contains("skipped"),"Unavailable YX prevented valid XY preview");
  validity.setData({}); check(!apply->isEnabled() && !undo->isEnabled(),"Unload retained masking actions");
}
void batchMasking(const std::shared_ptr<MTSurveyData> &original) {
  auto survey=std::make_shared<MTSurveyData>(*original);
  for(const auto &name:survey->get_stations_names()) {
    auto &s=survey->get_station_data(name); const bool active=s.active(); s.set_active(true);
    double real,error; check(s.scalar_value(RealZxy,5,real,error),"Missing batch fixture");
    s.set_data(s.frequencies()[5],{RealZxy,ImagZxy,RealZyx,ImagZyx},{real,-real,-real,real},{error,error,error,error});
    s.set_data_mask(PTxx,s.frequencies()[5],false);
    s.set_data_mask(RealZxx,s.frequencies().front(),false);
    s.set_data_mask(RealZyy,s.frequencies().back(),false); s.set_active(active);
  }
  survey->get_station_data("B").mask_type(RealZyx,false);
  const auto before=*survey;
  unsigned notifications=0;
  ImpedanceValidityWindow *window=nullptr;
  ImpedanceValidityWindow validity(nullptr,{},[&]{++notifications; window->setData(survey);}); window=&validity;
  validity.setData(survey,"A"); validity.show(); settle();
  auto *preview=control<QCheckBox>(validity,"validityPreviewMasks");
  auto *station=control<QComboBox>(validity,"validityStation");
  auto *all=control<QPushButton>(validity,"validityApplyAllMasks");
  auto *undo=control<QPushButton>(validity,"validityUndoMasks");
  auto *apply=control<QPushButton>(validity,"validityApplyMasks");
  check(preview->isChecked() && apply->isEnabled(),"Automatic preview missing on open");
  station->setCurrentText("B"); settle();
  check(preview->isChecked() && apply->isEnabled(),"Station navigation did not refresh preview");
  preview->setChecked(false);
  check(!apply->isEnabled() && control<QCustomPlot>(validity,"analysisRhoPlot")->findChildren<QCPGraph *>("rhoPlusMaskCandidates").empty(),"Unchecked preview retained candidates");
  station->setCurrentText("A"); settle(); check(!apply->isEnabled(),"Station navigation re-enabled unchecked preview");
  preview->setChecked(true); check(apply->isEnabled(),"Re-enabling preview did not show candidates");
  bool cancelled=false; QTimer cancel;
  QObject::connect(&cancel,&QTimer::timeout,&validity,[&]{
    if(auto *p=validity.findChild<QProgressDialog *>("validityMaskProgress")) if(p->value()>=1) { cancelled=true; p->cancel(); cancel.stop(); }
  });
  cancel.start(0); all->click(); cancel.stop();
  check(cancelled && notifications==0,"Batch cancellation did not stop before mutation");
  auto unchanged=[&]{
    for(const auto &name:survey->get_stations_names()) {
      const auto &s=survey->get_station_data(name), &old=before.get_station_data(name);
      check(s.impedance_mask()==old.impedance_mask() && s.phase_tensor_mask()==old.phase_tensor_mask() && s.tipper_mask()==old.tipper_mask(),"Batch masks not restored exactly");
    }
  };
  unchanged();
  all->click(); check(notifications==1 && undo->isEnabled(),"Batch apply or undo unavailable");
  for(const auto &name:{"A","B"}) {
    const auto &s=survey->get_station_data(name), &old=before.get_station_data(name);
    check(!s.impedance_mask()[1][5] && !s.phase_tensor_mask()[1][5],"Batch skipped enabled station");
    for(unsigned f=0;f<s.frequencies().size();++f) for(const auto &pair:{std::pair<unsigned,unsigned>{0,1},{3,2}}) {
      const auto diagonal=pair.first, off=pair.second;
      const bool newlyMasked=old.impedance_mask()[off][f] && !s.impedance_mask()[off][f];
      check(s.impedance_mask()[diagonal][f]==(old.impedance_mask()[diagonal][f] && !newlyMasked),"Batch paired the wrong diagonal or period");
    }
    check(s.tipper_mask()==old.tipper_mask(),"Batch changed tippers");
  }
  auto summary=control<QLabel>(validity,"validityMaskStatus")->text();
  check(summary.contains("2 stations") && summary.contains("1 disabled stations") && summary.contains("1 skipped fits"),"Batch summary counts incorrect");
  survey->get_station_data("DISABLED").set_active(true);
  check(survey->get_station_data("DISABLED").impedance_mask()[1][5],"Batch masked a disabled station");
  survey->get_station_data("DISABLED").set_active(false);
  station->setCurrentText("B"); settle(); check(undo->isEnabled(),"Navigation lost batch undo");
  undo->click(); check(notifications==2,"Batch undo did not notify plots"); unchanged();
  validity.setLinkTensorMasks(false);
  control<QCheckBox>(validity,"validityYX")->setChecked(false); settle(); all->click();
  for(const auto &name:{"A","B"}) {
    const auto &s=survey->get_station_data(name), &old=before.get_station_data(name);
    check(!s.impedance_mask()[0][5],"Unlinked batch did not mask the paired diagonal");
    check(s.phase_tensor_mask()==old.phase_tensor_mask() && s.impedance_mask()[2]==old.impedance_mask()[2] && s.impedance_mask()[3]==old.impedance_mask()[3],"Batch ignored link or component controls");
  }
  undo->click(); unchanged();
  control<QCheckBox>(validity,"validityXY")->setChecked(false); check(!all->isEnabled(),"Batch enabled without components");
  validity.setData({}); check(!all->isEnabled(),"Batch enabled without survey");
}
void dataErrors(const std::shared_ptr<MTSurveyData> &original) {
  auto survey=std::make_shared<MTSurveyData>(*original);
  ImpedanceValidityWindow validity;
  for(double floor:{.02,.08}) {
    survey->set_error_floor(floor);
    const auto &station=survey->get_station_data("A");
    std::vector<dvector> rho,rhoErrors,phase,phaseErrors;
    station.get_apparent_resistivity(rho,rhoErrors); station.get_phase(phase,phaseErrors);
    validity.setData(survey,"A"); validity.show(); settle();
    check(!validity.findChild<QDoubleSpinBox *>("validityErrorFloor"),"Independent rho+ floor still exposed");
    for(unsigned c:{1u,2u}) {
      const auto samples=ImpedanceAnalysis::curve(station,ImpedanceAnalysis::Quantity(c));
      for(unsigned i=0;i<samples.size();++i) {
        const unsigned f=samples.size()-1-i;
        check(samples[i].rhoError==rhoErrors[c][f] && samples[i].phaseError==phaseErrors[c][f],"Analysis changed main-data uncertainties");
      }
    }
    for(bool isRho:{true,false}) {
      auto *plot=control<QCustomPlot>(validity,isRho?"analysisRhoPlot":"analysisPhasePlot");
      const auto bars=plot->findChildren<QCPErrorBars *>(); check(bars.size()==2,"Missing component error bars");
      for(unsigned c=0;c<2;++c) for(int i=0;i<bars[c]->data()->size();++i) {
        const auto &errors=isRho?rhoErrors:phaseErrors;
        check(bars[c]->data()->at(i).errorPlus==errors[c+1][errors[c+1].size()-1-i],"Displayed error bars differ from main data");
      }
    }
    control<QPushButton>(validity,"validityRun")->click();
    std::vector<dvector> values,after;
    station.get_apparent_resistivity(values,after); check(after==rhoErrors,"Rho+ overwrote resistivity errors");
    station.get_phase(values,after); check(after==phaseErrors,"Rho+ overwrote phase errors");
  }
}
}
int main(int argc,char **argv) {
  QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc,argv);
  try {
    QTemporaryDir directory; check(directory.isValid(),"No temporary directory");
    const auto path=directory.path()+"/synthetic.gofem";
    QFile file(path); check(file.open(QIODevice::WriteOnly|QIODevice::Text),"Cannot create synthetic survey");
    QTextStream out(&file); out.setRealNumberPrecision(17);
    const double pi=std::acos(-1.);
    for(const auto &station: {QString("A"),QString("B"),QString("DISABLED")}) for(int i=0;i<12;++i) {
      const double f=std::pow(10.,-2.+i*4./11), rho=station=="B"?1000.:100., z=std::sqrt(pi*f*4e-7*pi*rho);
      for(const auto &part:{QString("Real"),QString("Imag")}) for(const auto &component:{QString("xx"),QString("xy"),QString("yx"),QString("yy")}) {
        const double value=component=="xy"?z:component=="yx"?-z:0.;
        out<<part<<'Z'<<component<<' '<<f<<" Plane_wave "<<station<<' '<<value<<' '<<.02*z<<'\n';
      }
    }
    out.flush(); file.close();
    auto survey=std::make_shared<MTSurveyData>(); survey->load_from_gofem(path.toStdString()); survey->set_active_flag("DISABLED",false);
    const auto mask=survey->get_station_data("A").impedance_mask();
    ImpedanceValidityWindow validity; validity.setData(survey,"A"); validity.show(); settle();
    auto *values=control<QTableWidget>(validity,"analysisValues");
    check(values->rowCount()==24,"Both off-diagonal fits should be listed");
    auto *quantity=control<QComboBox>(validity,"validityQuantity");
    check(quantity->currentData().toInt()==-1 && quantity->findData(0)<0 && quantity->findData(3)<0,"Diagonal rho+ choices remain");
    check(values->item(0,0)->text()=="Zxy" && values->item(12,0)->text()=="−Zyx","Table lost component labels");
    check(control<QLabel>(validity,"validityStatus")->text().contains("Fits within"),"Halfspace should fit");
    auto *plot=control<QCustomPlot>(validity,"analysisRhoPlot"); check(plot->graphCount()==4,"Combined rho+ overlays missing");
    for(unsigned c=0;c<2;++c) {
      check(plot->graph(2*c)->pen().color()==PlotColors::componentColors()[c+1] &&
        plot->graph(2*c+1)->pen().color()==PlotColors::componentColors()[c+1],"Observed and fitted component colors differ");
      check(plot->graph(2*c+1)->pen().style()==Qt::DashLine,"Fit line is not dashed");
      auto *phase=control<QCustomPlot>(validity,"analysisPhasePlot");
      check(std::abs(phase->graph(2*c)->data()->constBegin()->value-45.)<1e-8,"Off-diagonal phase conventions differ");
    }
    control<QCheckBox>(validity,"validityRhoPlus")->setChecked(false); check(plot->graphCount()==2,"Rho+ visibility ignored");
    control<QCheckBox>(validity,"validityRhoPlus")->setChecked(true);
    control<QCheckBox>(validity,"validityXY")->setChecked(false); settle();
    check(plot->graphCount()==2 && plot->graph(0)->name()=="−Zyx" && values->rowCount()==12,"Component checkbox ignored");
    control<QCheckBox>(validity,"validityXY")->setChecked(true); settle();
    const auto csvPath=directory.path()+"/validity.csv";
    QTimer::singleShot(0,[csvPath]{
      if(auto *dialog=qobject_cast<QFileDialog *>(QApplication::activeModalWidget())) {
        dialog->selectFile(csvPath); QMetaObject::invokeMethod(dialog,"accept",Qt::QueuedConnection);
      }
    });
    control<QPushButton>(validity,"validityCsv")->click();
    QFile csv(csvPath); check(csv.open(QIODevice::ReadOnly),"Could not read exported fit table");
    const auto exported=QString::fromUtf8(csv.readAll());
    check(exported.contains("\"Zxy\",") && exported.contains("\"−Zyx\",") && exported.count('\n')==25,"CSV lost component labels or rows");
    control<QComboBox>(validity,"validityStation")->setCurrentText("DISABLED"); settle();
    check(values->rowCount()==0,"Disabled station retained old fit");
    control<QComboBox>(validity,"validityStation")->setCurrentText("A"); settle();
    SurveyAverageWindow averages; averages.setData(survey); averages.show(); settle();
    auto *table=control<QTableWidget>(averages,"analysisValues");
    check(table->rowCount()==12,"Regional periods missing");
    check(table->item(0,1)->data(Qt::UserRole).toInt()==2,"Disabled station counted");
    check(std::abs(table->item(0,2)->data(Qt::UserRole).toDouble()-std::sqrt(100.*1000.))<1e-6,"Regional mean wrong");
    auto *stations=control<QListWidget>(averages,"averageStations"); stations->item(1)->setCheckState(Qt::Unchecked); settle();
    check(std::abs(table->item(0,2)->data(Qt::UserRole).toDouble()-100.)<1e-6,"Unchecked station counted");
    control<QSpinBox>(averages,"averageMinimumStations")->setValue(2); settle();
    check(std::isnan(table->item(0,2)->data(Qt::UserRole).toDouble()),"Low coverage mean visible");
    check(survey->get_station_data("A").impedance_mask()==mask && survey->is_active("B"),"Analysis modified the survey");
    control<QSpinBox>(averages,"averageMinimumStations")->setValue(1);
    stations->item(1)->setCheckState(Qt::Checked); settle();
    if(argc>1) { validity.grab().save(QString::fromLocal8Bit(argv[1])+"/validity.png"); averages.grab().save(QString::fromLocal8Bit(argv[1])+"/averages.png"); }
    validity.setData({}); averages.setData({}); settle();
    check(values->rowCount()==0 && table->rowCount()==0,"Unloading survey retained stale results");
    masking(survey,argc>1 ? QString::fromLocal8Bit(argv[1]) : QString());
    batchMasking(survey);
    dataErrors(survey);
    std::cout<<"Impedance analysis dialog checks passed.\n";
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
