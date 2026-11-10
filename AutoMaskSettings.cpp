#include "AutoMaskSettings.h"
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QSpinBox>
#include <QToolButton>
#include <QVBoxLayout>

QStringList AutoMaskSettings::criterionNames() {
  return {QObject::tr("Rho+ / dispersion"), QObject::tr("Period outliers"), QObject::tr("Neighbours"), QObject::tr("Uncertainty")};
}
QStringList AutoMaskSettings::criterionHelp() {
  return {
    QObject::tr("Joint resistivity/phase consistency using the existing rho+ fit and main-data propagated errors. Only converged, grid-stable fits contribute. Rho+ already constrains dispersion; it is one criterion, not an independent Kramers–Kronig vote. This is a 1-D consistency test, not a general 3-D validity verdict."),
    QObject::tr("Leave-one-out robust linear trend in log(period), using two samples on each side. Tests log resistivity and circular phase; the larger standardized residual is one vote. Tippers use real and imaginary values, with no log or angle transform. Includes reference scatter and measured errors. No extrapolation, no missing/masked reference samples, no large period gaps."),
    QObject::tr("Impedance phase, or complex tipper values, using enabled nearby stations throughout the survey. Requires at least three matching neighbours and 75% agreement: circular phase for impedance, complex distance for tippers. Includes measured errors and robust spatial scatter. Unknown coordinates, ambiguous periods or weak agreement provide no vote. Geological differences can be real; review candidates."),
    QObject::tr("Flag poorly resolved data when the main-data resistivity error, phase error, or absolute tipper error exceeds its limit. Unlike normalized residuals, large errors increase this score. Includes endpoints and isolated periods; requires no fit or neighbours. Limits are in Settings. The score and threshold are multiples of the limits, not sigma. Data and error floors stay unchanged.")};
}
AutoMaskSettings::AutoMaskSettings(const AutoMask::Options &options, QWidget *parent) : QWidget(parent), initial(options) {
  auto *layout = new QVBoxLayout(this); layout->setContentsMargins(0, 0, 0, 0);
  auto *row = new QHBoxLayout; layout->addLayout(row);
  const QStringList labels{"XY", "YX", "Tzx", "Tzy"};
  for(unsigned c = 0; c < components.size(); ++c) {
    auto *box = new QCheckBox(labels[c], this); box->setObjectName("autoMask" + labels[c]);
    box->setChecked(initial.components[c]); components[c] = box; row->addWidget(box);
  }
  auto *form = new QFormLayout; layout->addLayout(form);
  const auto criterionNames = this->criterionNames(), help = criterionHelp();
  for(unsigned k = 0; k < AutoMask::criterionCount; ++k) {
    criteria[k] = new QCheckBox(criterionNames[k], this); criteria[k]->setObjectName("autoMaskCriterion" + QString::number(k));
    criteria[k]->setChecked(initial.criteria[k]); criteria[k]->setToolTip(help[k]);
    thresholds[k] = new QDoubleSpinBox(this); thresholds[k]->setObjectName("autoMaskThreshold" + QString::number(k));
    thresholds[k]->setRange(.1, 100.); thresholds[k]->setValue(initial.threshold[k]); thresholds[k]->setSuffix(k == 3 ? " × limit" : " σ"); thresholds[k]->setToolTip(help[k]);
    form->addRow(criteria[k], thresholds[k]);
  }
  votes = new QSpinBox(this); votes->setObjectName("autoMaskMinimumFlags"); votes->setRange(1, AutoMask::criterionCount); votes->setValue(initial.minimumFlags);
  votes->setToolTip(QObject::tr("A suggested mask needs at least this many enabled criteria to flag the same component and period. Unsupported tests do not count. The default is one: any enabled criterion can nominate points for review. Increase it to require agreement between checks; this can exclude very noisy points that only fail the uncertainty check. You can override individual suggestions in the table."));
  form->addRow(QObject::tr("Minimum flags"), votes);
  auto *advancedButton = new QToolButton(this); advancedButton->setText(QObject::tr("Settings")); advancedButton->setCheckable(true);
  advancedButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon); advancedButton->setArrowType(Qt::RightArrow); layout->addWidget(advancedButton);
  auto *advanced = new QWidget(this); auto *settings = new QFormLayout(advanced); settings->setContentsMargins(0, 0, 0, 0);
  auto spin = [&](const char *name, const QString &label, double low, double high, double value, const QString &tip) {
    auto *box = new QDoubleSpinBox(advanced); box->setObjectName(name); box->setRange(low, high); box->setValue(value);
    box->setKeyboardTracking(false); box->setToolTip(tip); settings->addRow(label, box); return box;
  };
  gap = spin("autoMaskGap", QObject::tr("Max period ratio"), 1.01, 100., initial.maximumGapRatio, QObject::tr("Maximum adjacent-period ratio within the five-sample local window. Each sample must be distinct. Large gaps invalidate the local test."));
  radius = spin("autoMaskRadius", QObject::tr("Radius (km)"), .01, 20000., initial.radiusKm, QObject::tr("Maximum WGS84 distance to reference stations."));
  neighbours = new QSpinBox(advanced); neighbours->setObjectName("autoMaskNeighbours"); neighbours->setRange(3, 50); neighbours->setValue(initial.neighbours); settings->addRow(QObject::tr("Neighbours"), neighbours);
  tolerance = spin("autoMaskTolerance", QObject::tr("Period match (%)"), 0., 10., 100.*initial.periodTolerance, QObject::tr("Relative period tolerance. Matches must be unique on both sides; no frequency interpolation."));
  departure = spin("autoMaskDeparture", QObject::tr("Min phase offset (°)"), .1, 180., initial.minimumPhaseDifference, QObject::tr("Spatial phase departure must meet this minimum as well as exceed the sigma threshold."));
  agreement = spin("autoMaskAgreement", QObject::tr("Phase agreement (°)"), .1, 60., initial.phaseAgreement, QObject::tr("At least three and 75% of matched neighbours must agree within this angle of their circular medoid."));
  rhoError = spin("autoMaskRhoError", QObject::tr("Max ρ error (%)"), .1, 10000., 100.*initial.maximumRelativeRhoError,
    QObject::tr("Base limit for the main-data one-sigma resistivity error divided by resistivity. The uncertainty threshold multiplies this limit. This filters poor precision without changing error floors."));
  phaseError = spin("autoMaskPhaseError", QObject::tr("Max phase error (°)"), .1, 360., initial.maximumPhaseError,
    QObject::tr("Base limit for the main-data one-sigma phase error. Either resistivity or phase error exceeding its limit supplies one uncertainty vote, including at endpoints."));
  density = new QSpinBox(advanced); density->setObjectName("autoMaskDensity"); density->setRange(4, 48); density->setValue(initial.rhoPlus.intervalsPerDecade);
  density->setToolTip(QObject::tr("Rho+ intervals per decade. Error floors always come from the main data and are never changed here.")); settings->addRow(QObject::tr("Rho+ grid"), density);
  tipperError = spin("autoMaskTipperError", QObject::tr("Max tipper error"), .001, 10., initial.maximumTipperError,
    QObject::tr("Absolute one-sigma error limit for either part of Tzx/Tzy. Does not divide by amplitude: zero tippers remain valid."));
  tipperError->setDecimals(3); tipperError->setValue(initial.maximumTipperError);
  tipperDeparture = spin("autoMaskTipperDeparture", QObject::tr("Min tipper offset"), .001, 10., initial.minimumTipperDifference,
    QObject::tr("Minimum complex tipper difference from neighbours, in addition to the standardized residual threshold."));
  tipperDeparture->setDecimals(3); tipperDeparture->setValue(initial.minimumTipperDifference);
  tipperAgreement = spin("autoMaskTipperAgreement", QObject::tr("Tipper agreement"), .001, 10., initial.tipperAgreement,
    QObject::tr("At least three and 75% of matched neighbours must agree within this complex tipper distance."));
  tipperAgreement->setDecimals(3); tipperAgreement->setValue(initial.tipperAgreement);
  advanced->hide(); layout->addWidget(advanced);
  QObject::connect(advancedButton, &QToolButton::toggled, advanced, [advanced, advancedButton](bool on) {
    advanced->setVisible(on); advancedButton->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
  });

  const auto enable = [this] {
    for(unsigned k = 0; k < criteria.size(); ++k) thresholds[k]->setEnabled(criteria[k]->isChecked());
    density->setEnabled(criteria[0]->isChecked()); gap->setEnabled(criteria[1]->isChecked());
    for(QWidget *w : std::vector<QWidget *>{radius, neighbours, tolerance, departure, agreement, tipperDeparture, tipperAgreement}) w->setEnabled(criteria[2]->isChecked());
    for(auto *w : {rhoError, phaseError, tipperError}) w->setEnabled(criteria[3]->isChecked());
  };
  enable();
  for(auto *box : criteria) connect(box, &QCheckBox::toggled, this, [this, enable] { enable(); emit changed(); });
  for(auto *box : components) connect(box, &QCheckBox::toggled, this, &AutoMaskSettings::changed);
  for(auto *box : thresholds) connect(box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &AutoMaskSettings::changed);
  for(auto *box : {gap, radius, tolerance, departure, agreement, rhoError, phaseError, tipperError, tipperDeparture, tipperAgreement})
    connect(box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &AutoMaskSettings::changed);
  for(auto *box : {votes, neighbours, density}) connect(box, QOverload<int>::of(&QSpinBox::valueChanged), this, &AutoMaskSettings::changed);
}
AutoMask::Options AutoMaskSettings::options() const {
  auto result = initial;
  for(unsigned k = 0; k < criteria.size(); ++k) { result.criteria[k] = criteria[k]->isChecked(); result.threshold[k] = thresholds[k]->value(); }
  for(unsigned c = 0; c < components.size(); ++c) result.components[c] = components[c]->isChecked();
  result.minimumFlags = votes->value(); result.neighbours = neighbours->value(); result.rhoPlus.intervalsPerDecade = density->value();
  result.maximumGapRatio = gap->value(); result.radiusKm = radius->value(); result.periodTolerance = tolerance->value()/100.;
  result.minimumPhaseDifference = departure->value(); result.phaseAgreement = agreement->value();
  result.maximumRelativeRhoError = rhoError->value()/100.; result.maximumPhaseError = phaseError->value();
  result.maximumTipperError = tipperError->value(); result.minimumTipperDifference = tipperDeparture->value(); result.tipperAgreement = tipperAgreement->value();
  return result;
}
bool editAutoMaskSettings(QWidget *parent, AutoMask::Options &options) {
  QDialog dialog(parent); dialog.setObjectName("autoMaskSettingsDialog"); dialog.setWindowTitle(QObject::tr("Auto-mask settings"));
  auto *layout = new QVBoxLayout(&dialog); auto *editor = new AutoMaskSettings(options, &dialog); layout->addWidget(editor);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog); layout->addWidget(buttons);
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
  if(dialog.exec() != QDialog::Accepted) return false;
  options = editor->options(); return true;
}
