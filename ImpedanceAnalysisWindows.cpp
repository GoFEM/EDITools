#include "ImpedanceAnalysisWindows.h"
#include "include/ImpedanceAnalysis.h"
#include "include/MTSurveyData.h"
#include "include/SurveyHistory.h"
#include "include/MTDataPlot.h"
#include "include/ToolTips.h"
#include "qcustomplot.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QHeaderView>
#include <QLabel>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QProgressDialog>
#include <QSaveFile>
#include <QShowEvent>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QTextStream>
#include <QTimer>
#include <QVBoxLayout>
#include <cmath>

namespace {
using namespace ImpedanceAnalysis;
const double pi = std::acos(-1.);
QString number(double x) { return std::isfinite(x) ? QString::number(x, 'g', 8) : QStringLiteral("—"); }
QComboBox *quantities(QWidget *parent, const char *objectName, bool averages) {
  auto *box = new QComboBox(parent); box->setObjectName(objectName);
  if(!averages) box->addItem(QObject::tr("XY / YX"), -1);
  for(unsigned q = 4; q <= 6; ++q) box->addItem(QString::fromUtf8(name(Quantity(q))), int(q));
  ToolTips::set(box, QObject::tr("SSQ = sqrt((Zxx² + Zxy² + Zyx² + Zyy²)/2), using complex squares. Berdichevsky = (Zxy − Zyx)/2. Determinant = sqrt(Zxx Zyy − Zxy Zyx). Square roots use the branch with nonnegative imaginary part (0–180°). Required components must be finite and unmasked; missing diagonals are never replaced by zero. YX is sign-reversed for comparison with XY."));
  return box;
}
struct Plots {
  QCustomPlot *rho, *phase;
  explicit Plots(QWidget *parent) {
    rho = new QCustomPlot(parent); phase = new QCustomPlot(parent);
    rho->setObjectName("analysisRhoPlot"); phase->setObjectName("analysisPhasePlot");
    for(auto *p: {rho, phase}) {
      p->setMinimumSize(300, 240); p->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
      p->xAxis->setScaleType(QCPAxis::stLogarithmic);
      p->xAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
      p->xAxis->setLabel(QObject::tr("Period (s)"));
    }
    rho->yAxis->setScaleType(QCPAxis::stLogarithmic);
    rho->yAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    rho->yAxis->setLabel(QObject::tr("Apparent resistivity (Ω m)")); phase->yAxis->setLabel(QObject::tr("Phase (°)"));
  }
  void clear() { rho->clearPlottables(); phase->clearPlottables(); }
  void finish() {
    for(auto *p: {rho, phase}) {
      if(p->graphCount()) { p->rescaleAxes(); p->yAxis->scaleRange(1.15); }
      else { p->xAxis->setRange(.01,1000.); p->yAxis->setRange(p == rho ? 1. : 0., p == rho ? 1000. : 90.); }
      p->replot();
    }
  }
  void add(const QVector<double> &periods, const QVector<double> &resistivity, const QVector<double> &degrees,
           const QColor &color, const QString &label, bool points = false, double width = 1.) {
    for(auto *p: {rho, phase}) {
      auto *g = p->addGraph(); g->setName(label); g->setPen(QPen(color,width));
      g->setLineStyle(points ? QCPGraph::lsNone : QCPGraph::lsLine);
      if(points) g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, color, Qt::white, 5));
      QVector<double> values = p == rho ? resistivity : degrees;
      if(!points && std::count_if(values.begin(),values.end(),[](double v){return std::isfinite(v);}) == 1)
        g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, color, color, 5));
      if(p == phase && !points) for(int j = values.size()-1; j > 0; --j)
        if(std::abs(values[j] - values[j-1]) > 180.) values[j] = NAN;
      g->setData(periods,values);
    }
  }
  void add(const std::vector<Sample> &samples, const QColor &color, const QString &label, bool points = false, double width = 1.) {
    QVector<double> x, y, p;
    for(const auto &s: samples) { x.push_back(s.period); y.push_back(s.resistivity()); p.push_back(s.phase()); }
    add(x,y,p,color,label,points,width);
  }
};
QTableWidget *table(QWidget *parent, const QStringList &headers) {
  auto *t = new QTableWidget(parent); t->setObjectName("analysisValues"); t->setColumnCount(headers.size());
  t->setHorizontalHeaderLabels(headers); t->setEditTriggers(QAbstractItemView::NoEditTriggers);
  t->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch); t->verticalHeader()->hide(); return t;
}
void row(QTableWidget *t, int index, const std::vector<double> &values) {
  for(unsigned j = 0; j < values.size(); ++j) {
    auto *item = new QTableWidgetItem(number(values[j])); item->setData(Qt::UserRole, values[j]);
    item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter); t->setItem(index,j,item);
  }
}
void exportTable(QWidget *parent, QTableWidget *t, const QString &suggestion) {
  const auto path = QFileDialog::getSaveFileName(parent, QObject::tr("Export analysis"), suggestion, QObject::tr("CSV (*.csv)"));
  if(path.isEmpty()) return;
  QSaveFile file(path);
  if(!file.open(QIODevice::WriteOnly | QIODevice::Text)) { QMessageBox::warning(parent,QObject::tr("Export analysis"),file.errorString()); return; }
  QTextStream out(&file); out.setEncoding(QStringConverter::Utf8); out.setLocale(QLocale::c());
  for(int c=0;c<t->columnCount();++c) { if(c) out << ','; out << '"' << t->horizontalHeaderItem(c)->text() << '"'; } out << '\n';
  for(int r=0;r<t->rowCount();++r) {
    for(int c=0;c<t->columnCount();++c) {
      if(c) out << ',';
      const auto *item=t->item(r,c); if(!item) continue;
      const auto value=item->data(Qt::UserRole);
      if(value.metaType().id()==QMetaType::QString) out << '"' << QString(value.toString()).replace('"',"\"\"") << '"';
      else if(std::isfinite(value.toDouble())) out << QString::number(value.toDouble(),'g',17);
    }
    out << '\n';
  }
  out.flush(); if(!file.commit()) QMessageBox::warning(parent,QObject::tr("Export analysis"),file.errorString());
}
QCheckBox *visibility(QWidget *parent, QHBoxLayout *row, const QString &text, const char *name, const QColor &color) {
  auto *box = new QCheckBox(text,parent); box->setObjectName(name); box->setChecked(true);
  auto palette=box->palette(); palette.setColor(QPalette::WindowText,color); box->setPalette(palette); row->addWidget(box); return box;
}
}
struct ImpedanceValidityWindow::Impl {
  ImpedanceValidityWindow *owner;
  std::function<void(const RhoPlusOptions &)> optionsChanged;
  std::function<void()> masksChanged;
  std::shared_ptr<MTSurveyData> survey;
  QComboBox *station, *quantity;
  QSpinBox *density;
  QCheckBox *observed, *response, *errorBars, *showXY, *showYX;
  QLabel *status;
  QTableWidget *values;
  QTimer timer;
  Plots plots;
  struct Series {
    Quantity quantity;
    std::vector<Sample> data;
    RhoPlusResult result;
    std::vector<unsigned> candidates;
    int firstRow=0;
    QColor color() const {
      return quantity==Quantity::XY || quantity==Quantity::YX ? PlotColors::componentColors()[unsigned(quantity)] : QColor("#2563eb");
    }
    QString label() const { return QString::fromUtf8(name(quantity)); }
  };
  std::vector<Series> series;
  QDoubleSpinBox *threshold;
  QCheckBox *preview;
  QPushButton *apply, *applyAll, *undo;
  unsigned generation=0;
  QLabel *maskStatus;
  bool linkTensorMasks=true, editingMasks=false;
  struct MaskEdit { RealDataType type; double frequency; };
  std::vector<MaskEdit> previewEdits;
  struct MaskState {
    std::vector<double> frequencies;
    std::vector<std::vector<bool>> z, pt;
    explicit MaskState(const MTStationData &s) : frequencies(s.frequencies()),z(s.impedance_mask()),pt(s.phase_tensor_mask()) {}
    bool matches(const MTStationData &s) const {
      return s.active() && frequencies==s.frequencies() && z==s.impedance_mask() && pt==s.phase_tensor_mask();
    }
  };
  std::unique_ptr<MaskState> previewState;
  struct StationEdits {
    std::string name;
    std::vector<MaskEdit> edits;
    MaskState state;
    StationEdits(const MTStationData &s, std::vector<MaskEdit> changes)
      : name(s.name()), edits(std::move(changes)), state(s) {}
  };
  SurveyHistory localHistory{[this] { notifyMasksChanged(); }};
  SurveyHistory *history = &localHistory;
  bool maskable() const { return quantity->currentData().toInt()==-1; }
  void clearPreview() {
    for(auto &s:series) s.candidates.clear();
    previewEdits.clear(); previewState.reset(); apply->setEnabled(false);
    maskStatus->setText(QObject::tr("Preview checked components."));
  }
  void invalidate() {
    ++generation;
    clearPreview(); preview->setEnabled(maskable() && (showXY->isChecked() || showYX->isChecked()));
    applyAll->setEnabled(survey && survey->n_stations()>0 && preview->isEnabled());
    showXY->setVisible(maskable()); showYX->setVisible(maskable());
    series.clear(); values->setRowCount(0); draw();
  }
  explicit Impl(ImpedanceValidityWindow *w) : owner(w), plots(w) {
    auto *layout = new QVBoxLayout(w); auto *controls = new QHBoxLayout;
    station = new QComboBox(w); station->setObjectName("validityStation");
    quantity = quantities(w,"validityQuantity",false);
    controls->addWidget(new QLabel(QObject::tr("Station:"),w)); controls->addWidget(station,1);
    controls->addWidget(new QLabel(QObject::tr("Response:"),w)); controls->addWidget(quantity);
    density = new QSpinBox(w); density->setObjectName("validityGridDensity"); density->setRange(4,48); density->setValue(12); density->setKeyboardTracking(false);
    controls->addWidget(ToolTips::label(density,QObject::tr("Grid/decade:"),QObject::tr("Spectral intervals per decade. Rho+ uses bounded least squares with exact interval integrals and both spectral tails, then repeats with twice this density. Increase if the fit is grid-sensitive. Finite-grid results are approximations."))); controls->addWidget(density);
    layout->addLayout(controls);
    auto *legend = new QHBoxLayout;
    showXY=visibility(w,legend,QObject::tr("XY"),"validityXY",PlotColors::componentColors()[1]);
    showYX=visibility(w,legend,QObject::tr("YX"),"validityYX",PlotColors::componentColors()[2]);
    for(auto *box:{showXY,showYX}) ToolTips::set(box,QObject::tr("Show this component's observations, fit and table rows, and include it in masking previews. YX uses −Zyx so its phase can be compared directly with XY."));
    observed=visibility(w,legend,QObject::tr("Observed"),"validityObserved",QColor("#333333"));
    response=visibility(w,legend,QObject::tr("Rho+"),"validityRhoPlus",QColor("#333333"));
    errorBars = new QCheckBox(QObject::tr("Data errors"),w); errorBars->setObjectName("validityErrors"); errorBars->setChecked(true); legend->addWidget(errorBars);
    ToolTips::set(errorBars,QObject::tr("Propagated resistivity and phase errors from the main data, including its error floor. Rho+ does not add or change error floors. Missing or non-positive errors are excluded from fitting and masking."));
    ToolTips::set(response,QObject::tr("Show the fitted rho+ resistivity and phase. This tests consistency with a passive 1-D response, not validity of an arbitrary 3-D impedance. Diagonal components and tensor averages are not generally 1-D. Poor fits can reflect structure, conventions or underestimated errors, not only bad data. Only Apply masks changes observations' masks."));
    legend->addStretch(); auto *run=new QPushButton(QObject::tr("Recheck"),w); run->setObjectName("validityRun"); legend->addWidget(run);
    auto *csv=new QPushButton(QObject::tr("CSV…"),w); csv->setObjectName("validityCsv"); legend->addWidget(csv); layout->addLayout(legend);
    auto *maskControls=new QHBoxLayout;
    threshold=new QDoubleSpinBox(w); threshold->setObjectName("validityMaskThreshold"); threshold->setRange(.5,20.); threshold->setDecimals(1); threshold->setSingleStep(.5); threshold->setValue(3.); threshold->setKeyboardTracking(false);
    maskControls->addWidget(ToolTips::label(threshold,QObject::tr("Residual threshold (σ):"),QObject::tr("Flag a period if either |r lnρ| or |r φ| exceeds this threshold, using the main data’s propagated resistivity and phase errors. This is relative to measurement errors, not a confidence band around the fitted curve. 2σ is useful for review; 3σ is the default for masking. Each checked XY/YX component of the selected station is fitted and assessed independently."))); maskControls->addWidget(threshold);
    preview=new QCheckBox(QObject::tr("Preview masks"),w); preview->setObjectName("validityPreviewMasks"); preview->setChecked(true); maskControls->addWidget(preview);
    apply=new QPushButton(QObject::tr("Apply masks"),w); apply->setObjectName("validityApplyMasks"); apply->setEnabled(false); maskControls->addWidget(apply);
    applyAll=new QPushButton(QObject::tr("Apply all stations"),w); applyAll->setObjectName("validityApplyAllMasks"); applyAll->setEnabled(false); maskControls->addWidget(applyAll);
    undo=new QPushButton(QObject::tr("Undo masks"),w); undo->setObjectName("validityUndoMasks"); undo->setEnabled(false); maskControls->addWidget(undo);
    ToolTips::set(preview,QObject::tr("Automatically highlight candidates when stations, components or fit settings change. Requires a converged, grid-stable XY/YX fit. Preview never changes masks."));
    ToolTips::set(apply,QObject::tr("Mask the previewed component/periods in one pass, including both resistivity and phase. Each Zxy candidate also masks Zxx; each Zyx candidate also masks Zyy at the same period. Respects Data → Link Z / PT masks; when linked, all PT entries at each candidate period are masked. Tippers are unchanged. No repeated automatic rejection."));
    ToolTips::set(applyAll,QObject::tr("Fit the checked XY/YX components independently at every enabled station using the current threshold, grid density and main-data errors, then apply one masking pass. Zxy also masks Zxx; Zyx also masks Zyy. Disabled stations and unavailable or unstable fits are skipped. Cancel stops the batch without changing any masks. Undo restores the whole batch."));
    ToolTips::set(undo,QObject::tr("Undo the most recent survey edit. When opened from the main window, this uses the shared Edit history, including edits made in other windows."));
    maskControls->addStretch(); layout->addLayout(maskControls);
    maskStatus=new QLabel(w); maskStatus->setObjectName("validityMaskStatus"); maskStatus->setWordWrap(true); layout->addWidget(maskStatus);
    clearPreview();
    auto *split=new QSplitter(Qt::Vertical,w); auto *panels=new QWidget(split); auto *pair=new QHBoxLayout(panels); pair->setContentsMargins(0,0,0,0); pair->addWidget(plots.rho); pair->addWidget(plots.phase);
    values=table(split,{QObject::tr("Component"),QObject::tr("Period (s)"),QObject::tr("ρ obs. (Ω m)"),QObject::tr("φ obs. (°)"),QObject::tr("ρ+ (Ω m)"),QObject::tr("φ+ (°)"),QObject::tr("r lnρ"),QObject::tr("r φ")});
    values->setMinimumHeight(130); split->addWidget(panels); split->addWidget(values); split->setStretchFactor(0,3); split->setStretchFactor(1,1); layout->addWidget(split,1);
    for(int c: {6,7}) values->horizontalHeaderItem(c)->setToolTip(QObject::tr("Predicted minus observed, divided by the effective standard error. Values beyond ±2 are highlighted; this is a diagnostic, not automatic rejection."));
    status=new QLabel(w); status->setObjectName("validityStatus"); status->setWordWrap(true); layout->addWidget(status);
    ToolTips::set(status,QObject::tr("χ² is the sum of squared normalized ln(rho) and phase residuals; nRMS = sqrt(χ² / (2N)). This is not reduced chi-square or a formal p-value. A finite numerical grid may overestimate the best possible fit. Convergence and grid stability are reported separately. Rho+ curves are drawn only within the usable period range."));
    timer.setSingleShot(true); timer.setInterval(150);
    QObject::connect(&timer,&QTimer::timeout,w,[this]{ calculate(); });
    auto schedule=[this]{ invalidate(); if(owner->isVisible()) timer.start(); };
    for(auto *box:{showXY,showYX}) QObject::connect(box,&QCheckBox::toggled,w,[schedule](bool){ schedule(); });
    for(auto *box:{station,quantity}) QObject::connect(box,QOverload<int>::of(&QComboBox::currentIndexChanged),w,[schedule](int){ schedule(); });
    auto changed=[this,schedule]{
      RhoPlusOptions options; options.intervalsPerDecade=density->value();
      if(optionsChanged) optionsChanged(options);
      schedule();
    };
    QObject::connect(density,QOverload<int>::of(&QSpinBox::valueChanged),w,[changed](int){ changed(); });
    QObject::connect(run,&QPushButton::clicked,w,[this]{ timer.stop(); calculate(); });
    for(auto *box:{observed,response,errorBars}) QObject::connect(box,&QCheckBox::toggled,w,[this]{ draw(); });
    QObject::connect(csv,&QPushButton::clicked,w,[this]{ exportTable(owner,values,"impedance-validity.csv"); });
    QObject::connect(threshold,QOverload<double>::of(&QDoubleSpinBox::valueChanged),w,[schedule](double){ schedule(); });
    QObject::connect(preview,&QCheckBox::toggled,w,[this](bool){ refreshPreview(); draw(); });
    QObject::connect(applyAll,&QPushButton::clicked,w,[this]{ applyAllMasks(); });
    QObject::connect(apply,&QPushButton::clicked,w,[this]{ applyMasks(); });
    QObject::connect(undo,&QPushButton::clicked,w,[this]{ undoMasks(); });
    QObject::connect(&history->stack(), &QUndoStack::canUndoChanged, w, [this](bool on) { undo->setEnabled(on); });
  }
  void calculate() {
    clearPreview(); series.clear(); values->setRowCount(0);
    if(!survey || station->currentText().isEmpty() || !survey->is_station_present(station->currentText().toStdString())) { status->setText(QObject::tr("Load a survey to check impedance responses.")); draw(); return; }
    std::vector<Quantity> selected;
    if(maskable()) {
      if(showXY->isChecked()) selected.push_back(Quantity::XY);
      if(showYX->isChecked()) selected.push_back(Quantity::YX);
    } else selected.push_back(Quantity(quantity->currentData().toInt()));
    QStringList summary;
    for(auto q:selected) {
      Series entry; entry.quantity=q;
      entry.data=curve(survey->get_station_data(station->currentText().toStdString()),q);
      entry.firstRow=values->rowCount();
      auto &result=entry.result;
      try {
        RhoPlusOptions options; options.intervalsPerDecade=density->value();
        result=rhoPlus(entry.data,options);
        const QString outcome = !result.converged ? QObject::tr("Solver limit — provisional fit") : !result.gridStable ? QObject::tr("Grid-sensitive fit — increase resolution") :
          result.normalizedRms<=1. ? QObject::tr("Fits within the selected error level") : QObject::tr("Residuals exceed the selected error level");
        auto metric=[](double value){ return QString::number(value,'g',3); };
        summary.push_back(QObject::tr("%1: %2 · N=%3 · nRMS=%4 · χ²=%5 · Δχ²=%6 · %7 omitted · %8 invalid-error samples")
          .arg(entry.label(),outcome).arg(result.observed.size()).arg(metric(result.normalizedRms)).arg(metric(result.chiSquare)).arg(metric(result.refinementChange)).arg(result.excluded).arg(result.invalidErrors));
        values->setRowCount(entry.firstRow+result.observed.size());
        for(unsigned i=0;i<result.observed.size();++i) {
          const int r=entry.firstRow+i;
          row(values,r,{NAN,result.observed[i].period,result.observed[i].resistivity(),result.observed[i].phase(),result.fitted[i].resistivity(),result.fitted[i].phase(),result.rhoResidual[i],result.phaseResidual[i]});
          values->item(r,0)->setText(entry.label()); values->item(r,0)->setData(Qt::UserRole,entry.label());
          values->item(r,0)->setForeground(entry.color());
        }
      } catch(const std::exception &e) { summary.push_back(entry.label()+": "+QString::fromUtf8(e.what())); }
      series.push_back(std::move(entry));
    }
    status->setText(summary.isEmpty() ? QObject::tr("Check XY or YX to show a response.") : summary.join('\n'));
    refreshPreview(); draw();
  }
  std::vector<MaskEdit> collectEdits(const MTStationData &s, std::vector<Series> &entries, unsigned &linkedCount, unsigned &diagonalCount) {
    const MaskState state(s);
    std::vector<MaskEdit> edits;
    linkedCount=0; diagonalCount=0;
    std::set<std::pair<RealDataType,double>> uniqueEdits;
    for(auto &entry:entries) {
      const auto &result=entry.result;
      if(!result.converged || !result.gridStable) continue;
      entry.candidates=rhoPlusOutliers(result,threshold->value());
      const auto type=entry.quantity==Quantity::XY ? RealZxy : RealZyx;
      const auto diagonal=entry.quantity==Quantity::XY ? RealZxx : RealZyy;
      std::vector<RealDataType> types{type,diagonal};
      if(linkTensorMasks) for(auto related:MTStationData::linked_mask_types(type)) types.push_back(related);
      for(auto index:entry.candidates) {
        const auto &freqs=s.frequencies();
        const auto found=std::find_if(freqs.begin(),freqs.end(),[&](double f){return 1./f==result.observed[index].period;});
        if(found==freqs.end()) throw std::runtime_error("Candidate period is no longer present.");
        const auto f=std::distance(freqs.begin(),found);
        for(auto target:types) {
          const bool tensor=target!=type && target!=diagonal;
          const auto &masks=tensor ? state.pt : state.z;
          if(masks.at(type_to_column_table.at(target)).at(f) && uniqueEdits.insert({target,*found}).second) {
            edits.push_back({target,*found});
            if(tensor) ++linkedCount;
            else if(target==diagonal) ++diagonalCount;
          }
        }
      }
    }
    return edits;
  }
  void refreshPreview() {
    clearPreview();
    if(!preview->isChecked() || !maskable() || series.empty() || !survey) return;
    auto &s=survey->get_station_data(station->currentText().toStdString());
    previewState.reset(new MaskState(s));
    unsigned linkedCount=0, diagonalCount=0;
    previewEdits=collectEdits(s,series,linkedCount,diagonalCount);
    QStringList counts;
    for(const auto &entry:series) {
      if(!entry.result.converged || !entry.result.gridStable)
        counts.push_back(QObject::tr("%1: skipped (fit unavailable or unstable)").arg(entry.label()));
      else counts.push_back(QObject::tr("%1: %2/%3").arg(entry.label()).arg(entry.candidates.size()).arg(entry.result.observed.size()));
    }
    maskStatus->setText(QObject::tr("%1 candidates · %2 diagonal masks · %3 linked PT masks").arg(counts.join(" · ")).arg(diagonalCount).arg(linkedCount));
    apply->setEnabled(!previewEdits.empty());
  }
  bool unchangedCurve(const MTStationData &s) const {
    auto same=[](double a,double b){return a==b || (std::isnan(a) && std::isnan(b));};
    for(const auto &entry:series) {
      const auto current=curve(s,entry.quantity);
      const auto &data=entry.data;
      if(current.size()!=data.size()) return false;
      for(unsigned i=0;i<data.size();++i)
        if(!same(current[i].period,data[i].period) || !same(current[i].error,data[i].error) ||
           !same(current[i].rhoError,data[i].rhoError) || !same(current[i].phaseError,data[i].phaseError) ||
           !same(current[i].impedance.real(),data[i].impedance.real()) || !same(current[i].impedance.imag(),data[i].impedance.imag())) return false;
    }
    return true;
  }
  void notifyMasksChanged() {
    editingMasks=true;
    if(masksChanged) masksChanged();
    editingMasks=false;
    timer.stop(); calculate();
  }
  void applyMasks() {
    if(!previewState || previewEdits.empty() || !survey || !survey->is_station_present(station->currentText().toStdString())) return;
    auto &s=survey->get_station_data(station->currentText().toStdString());
    if(!previewState->matches(s) || !unchangedCurve(s)) {
      clearPreview(); draw(); maskStatus->setText(QObject::tr("Data changed; preview again.")); return;
    }
    unsigned count=0; for(const auto &entry:series) count+=entry.candidates.size();
    const auto changes=previewEdits;
    for(const auto &edit:changes) s.set_data_mask(edit.type,edit.frequency,false,false);
    history->record(QObject::tr("Rho+ masks"));
    notifyMasksChanged();
    maskStatus->setText(QObject::tr("Masked %1 off-diagonal points and their paired diagonals in one pass.").arg(count));
  }
  void applyAllMasks() {
    if(!survey || !maskable() || !(showXY->isChecked() || showYX->isChecked())) return;
    timer.stop();
    const auto source=survey;
    const auto revision=generation;
    const auto names=source->get_stations_names();
    std::vector<Quantity> selected;
    if(showXY->isChecked()) selected.push_back(Quantity::XY);
    if(showYX->isChecked()) selected.push_back(Quantity::YX);
    RhoPlusOptions options; options.intervalsPerDecade=density->value();
    QProgressDialog progress(QObject::tr("Fitting stations…"),QObject::tr("Cancel"),0,int(names.size()),owner);
    progress.setObjectName("validityMaskProgress"); progress.setWindowModality(Qt::ApplicationModal); progress.setMinimumDuration(0);
    progress.setAutoClose(false); progress.setAutoReset(false);
    std::vector<StationEdits> pending;
    unsigned disabled=0, skipped=0, candidates=0, linked=0, diagonals=0, index=0;
    for(const auto &name:names) {
      progress.setLabelText(QObject::tr("Fitting %1").arg(QString::fromStdString(name)));
      progress.setValue(index++);
      if(progress.wasCanceled() || generation!=revision || survey!=source) {
        maskStatus->setText(QObject::tr("Batch cancelled; no masks changed.")); return;
      }
      const auto &s=source->get_station_data(name);
      if(!s.active()) { ++disabled; continue; }
      std::vector<Series> fits;
      for(auto q:selected) {
        Series entry; entry.quantity=q;
        try { entry.result=rhoPlus(curve(s,q),options); } catch(const std::exception &) { ++skipped; continue; }
        if(!entry.result.converged || !entry.result.gridStable) { ++skipped; continue; }
        fits.push_back(std::move(entry));
      }
      unsigned ptCount=0, diagonalCount=0;
      auto edits=collectEdits(s,fits,ptCount,diagonalCount); linked+=ptCount; diagonals+=diagonalCount;
      for(const auto &entry:fits) candidates+=entry.candidates.size();
      if(!edits.empty()) pending.emplace_back(s,std::move(edits));
    }
    progress.setValue(names.size());
    if(progress.wasCanceled() || generation!=revision || survey!=source) {
      maskStatus->setText(QObject::tr("Batch cancelled; no masks changed.")); return;
    }
    // Stage every fit before mutating anything; validate the entire batch first.
    for(const auto &entry:pending) if(!source->is_station_present(entry.name) || !entry.state.matches(source->get_station_data(entry.name))) {
      maskStatus->setText(QObject::tr("Data changed; no masks applied.")); return;
    }
    for(auto &entry:pending) {
      auto &s=source->get_station_data(entry.name);
      for(const auto &edit:entry.edits) s.set_data_mask(edit.type,edit.frequency,false,false);
    }
    const auto changed=pending.size();
    if(changed) { history->record(QObject::tr("Rho+ masks (all stations)")); notifyMasksChanged(); }
    maskStatus->setText(QObject::tr("Masked %1 off-diagonal points across %2 stations · %3 diagonal masks · %4 linked PT masks · %5 disabled stations · %6 skipped fits")
      .arg(candidates).arg(changed).arg(diagonals).arg(linked).arg(disabled).arg(skipped));
  }
  void undoMasks() {
    history->stack().undo();
  }
  void draw() {
    plots.clear();
    // Restore residual highlighting after changing or dismissing a preview.
    for(int r=0;r<values->rowCount();++r) for(int c=0;c<values->columnCount();++c)
      if(auto *item=values->item(r,c)) item->setBackground(c>=6 && std::abs(item->data(Qt::UserRole).toDouble())>2. ? QBrush(QColor("#ffedd5")) : QBrush());
    for(const auto &entry:series) {
      const auto color=entry.color();
      if(observed->isChecked()) {
        plots.add(entry.data,color,entry.label(),true);
        if(errorBars->isChecked()) for(auto *plot:{plots.rho,plots.phase}) {
          QVector<double> errors;
          for(const auto &s:entry.data) {
            const double error=plot==plots.rho ? s.rhoError : s.phaseError;
            errors.push_back(s.usable() && std::isfinite(error) && error>0. ? error : NAN);
          }
          auto *bars=new QCPErrorBars(plot->xAxis,plot->yAxis); bars->setDataPlottable(plot->graph(plot->graphCount()-1)); bars->setData(errors); bars->setPen(QPen(color));
        }
      }
      if(response->isChecked() && !entry.result.response.empty()) {
        plots.add(entry.result.response,color,entry.label()+" ρ+",false,.9);
        for(auto *plot:{plots.rho,plots.phase}) plot->graph(plot->graphCount()-1)->setPen(QPen(color,.9,Qt::DashLine));
      }
      if(!entry.candidates.empty()) {
        std::vector<Sample> flagged;
        for(auto i:entry.candidates) {
          flagged.push_back(entry.result.observed[i]);
          for(int c=0;c<values->columnCount();++c) values->item(entry.firstRow+i,c)->setBackground(QColor("#fee2e2"));
        }
        plots.add(flagged,color,entry.label()+QObject::tr(" mask candidates"),true);
        for(auto *plot:{plots.rho,plots.phase}) {
          auto *g=plot->graph(plot->graphCount()-1); g->setObjectName("rhoPlusMaskCandidates");
          g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCross,QPen(color,2.),Qt::NoBrush,10.));
        }
      }
    }
    plots.finish();
  }

};
ImpedanceValidityWindow::ImpedanceValidityWindow(QWidget *parent, std::function<void(const RhoPlusOptions &)> optionsChanged,
  std::function<void()> masksChanged)
  :QDialog(parent,Qt::Window),impl(new Impl(this)) {
  impl->optionsChanged=std::move(optionsChanged);
  impl->masksChanged=std::move(masksChanged);
  setObjectName("impedanceValidityWindow"); setWindowTitle(tr("Impedance validity — rho+")); resize(1150,780);
}
ImpedanceValidityWindow::~ImpedanceValidityWindow() = default;
void ImpedanceValidityWindow::setHistory(SurveyHistory *history) {
  if(!history || history == impl->history) return;
  disconnect(&impl->history->stack(), nullptr, this, nullptr);
  impl->history=history;
  impl->undo->setText(tr("Undo"));
  impl->undo->setEnabled(history->stack().canUndo());
  connect(&history->stack(), &QUndoStack::canUndoChanged, this, [this](bool on) { impl->undo->setEnabled(on); });
}
void ImpedanceValidityWindow::setFitOptions(const RhoPlusOptions &options) {
  const QSignalBlocker densityBlock(impl->density);
  impl->density->setValue(options.intervalsPerDecade);
  impl->invalidate();
  if(isVisible()) impl->timer.start();
}
void ImpedanceValidityWindow::setLinkTensorMasks(bool linked) {
  if(impl->linkTensorMasks==linked) return;
  impl->linkTensorMasks=linked; impl->invalidate(); if(isVisible()) impl->timer.start();
}
void ImpedanceValidityWindow::setData(const std::shared_ptr<MTSurveyData> &survey, const QString &selected) {
  const QString previous=selected.isEmpty()?impl->station->currentText():selected;
  const QSignalBlocker block(impl->station); impl->survey=survey; impl->station->clear();
  if(!impl->editingMasks && impl->history == &impl->localHistory) impl->localHistory.reset(survey);
  if(survey) for(const auto &name:survey->get_stations_names()) impl->station->addItem(QString::fromStdString(name));
  if(impl->station->count()) impl->station->setCurrentIndex(std::max(0,impl->station->findText(previous)));
  impl->invalidate();
  if(isVisible()) impl->timer.start();
}
void ImpedanceValidityWindow::showEvent(QShowEvent *event) { QDialog::showEvent(event); impl->timer.start(); }

struct SurveyAverageWindow::Impl {
  SurveyAverageWindow *owner;
  std::shared_ptr<MTSurveyData> survey;
  QListWidget *stations;
  QComboBox *quantity;
  QDoubleSpinBox *tolerance;
  QSpinBox *minimum;
  QCheckBox *showStations,*showMean;
  QLabel *status;
  QTableWidget *values;
  QTimer timer;
  Plots plots;
  explicit Impl(SurveyAverageWindow *w):owner(w),plots(w) {
    auto *layout=new QVBoxLayout(w); auto *controls=new QHBoxLayout;
    quantity=quantities(w,"averageQuantity",true); controls->addWidget(new QLabel(QObject::tr("Response:"),w)); controls->addWidget(quantity);
    tolerance=new QDoubleSpinBox(w); tolerance->setObjectName("averageTolerance"); tolerance->setRange(0.,25.); tolerance->setDecimals(4); tolerance->setSingleStep(.01); tolerance->setKeyboardTracking(false);
    controls->addWidget(ToolTips::label(tolerance,QObject::tr("Tolerance (%):"),QObject::tr("Zero uses exact periods. Otherwise the entire group must fit this relative period span, without chaining. Each station contributes once; groups with repeated samples from a station omit that station. No interpolation or extrapolation. The group is shown at its lower median contributing period."))); controls->addWidget(tolerance);
    minimum=new QSpinBox(w); minimum->setObjectName("averageMinimumStations"); minimum->setRange(1,100000); minimum->setValue(1);
    controls->addWidget(ToolTips::label(minimum,QObject::tr("Min. stations:"),QObject::tr("Minimum number of checked, enabled stations with all required components unmasked and finite. Below this count the mean is omitted, leaving a gap."))); controls->addWidget(minimum); controls->addStretch(); layout->addLayout(controls);
    auto *legend=new QHBoxLayout;
    showStations=visibility(w,legend,QObject::tr("Station curves"),"averageShowStations",QColor("#737373"));
    showMean=visibility(w,legend,QObject::tr("Regional geometric mean"),"averageShowMean",QColor("#111827"));
    ToolTips::set(showMean,QObject::tr("Equal weight per checked station. Mean resistivity = exp(mean(ln(rho))); impedance phases are aligned to a common branch before averaging. With exact periods this is the complex geometric mean impedance. With a period tolerance, measured resistivities are averaged at their own periods. Determinant regional averaging is also known as the Berdichevsky regional average. A changing station count can change the curve independently of geology."));
    legend->addStretch(); auto *csv=new QPushButton(QObject::tr("CSV…"),w); legend->addWidget(csv); layout->addLayout(legend);
    auto *split=new QSplitter(w); auto *selection=new QWidget(split); auto *sl=new QVBoxLayout(selection); sl->setContentsMargins(0,0,0,0);
    stations=new QListWidget(selection); stations->setObjectName("averageStations"); stations->setMaximumWidth(270); sl->addWidget(stations);
    ToolTips::set(stations,QObject::tr("Check stations to include in both the displayed curves and regional mean. Disabled survey stations cannot contribute. These checkboxes do not change the survey."));
    auto *buttons=new QHBoxLayout;
    for(bool on:{true,false}) {
      auto *button=new QPushButton(on?QObject::tr("All"):QObject::tr("None"),selection); buttons->addWidget(button);
      QObject::connect(button,&QPushButton::clicked,w,[this,on]{ const QSignalBlocker block(stations); for(int i=0;i<stations->count();++i) if(stations->item(i)->flags().testFlag(Qt::ItemIsEnabled)) stations->item(i)->setCheckState(on?Qt::Checked:Qt::Unchecked); timer.start(); });
    }
    sl->addLayout(buttons);
    auto *right=new QSplitter(Qt::Vertical,split); auto *panels=new QWidget(right); auto *pair=new QHBoxLayout(panels); pair->setContentsMargins(0,0,0,0); pair->addWidget(plots.rho); pair->addWidget(plots.phase);
    values=table(right,{QObject::tr("Period (s)"),QObject::tr("Stations"),QObject::tr("Mean ρ (Ω m)"),QObject::tr("Mean φ (°)"),QObject::tr("SD lnρ")});
    values->horizontalHeaderItem(3)->setToolTip(QObject::tr("Impedance phases aligned to a common branch, then averaged. Undefined when phase directions cancel; the mean resistivity can still be shown."));
    values->horizontalHeaderItem(4)->setToolTip(QObject::tr("Sample standard deviation of log resistivity between stations, not uncertainty of the mean. Undefined for one station."));
    right->addWidget(panels); right->addWidget(values); right->setStretchFactor(0,3); right->setStretchFactor(1,1);
    split->addWidget(selection); split->addWidget(right); split->setStretchFactor(1,1); layout->addWidget(split,1);
    status=new QLabel(w); status->setObjectName("averageStatus"); status->setWordWrap(true); layout->addWidget(status);
    timer.setSingleShot(true); timer.setInterval(100); QObject::connect(&timer,&QTimer::timeout,w,[this]{ draw(); });
    QObject::connect(stations,&QListWidget::itemChanged,w,[this](QListWidgetItem *){timer.start();});
    QObject::connect(quantity,QOverload<int>::of(&QComboBox::currentIndexChanged),w,[this](int){timer.start();});
    QObject::connect(tolerance,QOverload<double>::of(&QDoubleSpinBox::valueChanged),w,[this](double){timer.start();});
    QObject::connect(minimum,QOverload<int>::of(&QSpinBox::valueChanged),w,[this](int){timer.start();});
    for(auto *box:{showStations,showMean}) QObject::connect(box,&QCheckBox::toggled,w,[this]{timer.start();});
    QObject::connect(csv,&QPushButton::clicked,w,[this]{exportTable(owner,values,"regional-average.csv");});
  }
  void draw() {
    plots.clear(); values->setRowCount(0);
    std::vector<std::vector<Sample>> curves;
    if(survey) for(int i=0;i<stations->count();++i) if(stations->item(i)->checkState()==Qt::Checked) {
      const auto name=stations->item(i)->text().toStdString(); if(!survey->is_station_present(name)) continue;
      curves.push_back(curve(survey->get_station_data(name),Quantity(quantity->currentData().toInt())));
      if(showStations->isChecked()) plots.add(curves.back(),QColor(140,140,140,130),QString::fromStdString(name));
    }
    const auto mean=regionalMean(curves,tolerance->value()/100.,minimum->value());
    QVector<double> periods,rho,phase; values->setRowCount(mean.size()); unsigned available=0;
    for(unsigned i=0;i<mean.size();++i) {
      const auto &v=mean[i]; periods.push_back(v.period); rho.push_back(v.resistivity); phase.push_back(v.phase);
      row(values,i,{v.period,double(v.stations),v.resistivity,v.phase,v.logRhoSpread}); if(std::isfinite(v.resistivity)) ++available;
    }
    if(showMean->isChecked() && !mean.empty()) plots.add(periods,rho,phase,QColor("#111827"),QObject::tr("Regional geometric mean"),false,2.5);
    status->setText(QObject::tr("%1 checked stations · %2 usable mean periods / %3 groups · equal station weights · no interpolation")
      .arg(curves.size()).arg(available).arg(mean.size()));
    plots.finish();
  }
};
SurveyAverageWindow::SurveyAverageWindow(QWidget *parent):QDialog(parent,Qt::Window),impl(new Impl(this)) {
  setObjectName("surveyAverageWindow"); setWindowTitle(tr("Survey impedance averages")); resize(1250,780);
}
SurveyAverageWindow::~SurveyAverageWindow() = default;
void SurveyAverageWindow::setData(const std::shared_ptr<MTSurveyData> &survey) {
  std::map<QString,bool> previous;
  if(survey==impl->survey) for(int i=0;i<impl->stations->count();++i) previous[impl->stations->item(i)->text()]=impl->stations->item(i)->checkState()==Qt::Checked;
  const QSignalBlocker block(impl->stations); impl->survey=survey; impl->stations->clear();
  if(survey) for(const auto &name:survey->get_stations_names()) {
    const auto text=QString::fromStdString(name); auto *item=new QListWidgetItem(text,impl->stations);
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    const bool enabled=survey->is_active(name); item->setCheckState(enabled && (!previous.count(text) || previous[text])?Qt::Checked:Qt::Unchecked);
    if(!enabled) {item->setFlags(item->flags() & ~Qt::ItemIsEnabled); item->setToolTip(tr("Disabled in the survey."));}
  }
  if(isVisible()) impl->timer.start();
}
void SurveyAverageWindow::showEvent(QShowEvent *event) { QDialog::showEvent(event); impl->timer.start(); }
