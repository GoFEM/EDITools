#include "ui/workspace/TFViewControls.h"
#include "ui/plots/PlotColorMaps.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSignalBlocker>
#include <cmath>

namespace {
QGridLayout *grid(QWidget *widget)
{
  auto *layout = new QGridLayout(widget);
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setSpacing(3);
  return layout;
}
void valid(QWidget *widget, bool ok, const QString &help)
{
  widget->setStyleSheet(ok ? QString() : "border: 1px solid #b91c1c;");
  widget->setToolTip(help);
  widget->setProperty("validInput", ok);
}
double bounded(double value, double minimum, double maximum, double fallback)
{
  return std::isfinite(value) ? qBound(minimum, value, maximum) : fallback;
}
}

TFViewControls::TFViewControls(QWidget *owner, std::function<TFViewRange(int)> range,
                             std::function<void(int, TFViewRange)> applyRange,
                             std::function<void()> appearanceChanged)
  : QObject(owner), range(std::move(range)), applyRange(std::move(applyRange)),
    appearanceChanged(std::move(appearanceChanged))
{
  ranges = new QWidget(owner);
  auto *layout = grid(ranges);
  axis = new QComboBox(ranges); axis->setObjectName("viewRangeAxis");
  axis->addItems({tr("Period (s)"), tr("Resistivity XY / YX (Ω m)"), tr("Phase XY / YX (°)"),
                  tr("Resistivity XX / YY (Ω m)"), tr("Phase XX / YY (°)")});
  axis->setAccessibleName(tr("Axis range"));
  layout->addWidget(axis, 0, 0, 1, 2);
  automatic = new QCheckBox(tr("Auto"), ranges); automatic->setObjectName("viewRangeAuto");
  automatic->setToolTip(tr("Fit this axis when changing stations. Turn off to retain a fixed range."));
  layout->addWidget(automatic, 0, 2);
  lower = new QLineEdit(ranges); upper = new QLineEdit(ranges);
  lower->setObjectName("viewRangeMin"); upper->setObjectName("viewRangeMax");
  lower->setPlaceholderText(tr("Minimum")); upper->setPlaceholderText(tr("Maximum"));
  lower->setAccessibleName(tr("Axis minimum")); upper->setAccessibleName(tr("Axis maximum"));
  for(auto *edit: {lower, upper}) edit->setFixedWidth(100);
  layout->addWidget(lower, 1, 0); layout->addWidget(upper, 1, 1);
  auto *apply = new QPushButton(tr("Apply"), ranges); apply->setObjectName("viewRangeApply");
  layout->addWidget(apply, 1, 2);
  const auto commitRange = [this] {
    bool loOk, hiOk;
    TFViewRange value{false, lower->text().toDouble(&loOk), upper->text().toDouble(&hiOk)};
    const bool logarithmic = axis->currentIndex() == 0 || axis->currentIndex() == 1 || axis->currentIndex() == 3;
    const bool ok = loOk && hiOk && std::isfinite(value.lower) && std::isfinite(value.upper) &&
      value.lower < value.upper && QCPRange::validRange(value.lower, value.upper) && (!logarithmic || value.lower > 0.);
    const auto help = tr("Enter minimum < maximum; period and resistivity must be positive. Scientific notation (1e-3) is accepted. Apply sets a fixed range.");
    for(auto *edit: {lower, upper}) valid(edit, ok, help);
    if(!ok) return;
    this->applyRange(axis->currentIndex(), value); refreshRange();
  };
  connect(apply, &QPushButton::clicked, this, commitRange);
  for(auto *edit: {lower, upper}) connect(edit, &QLineEdit::returnPressed, this, commitRange);
  connect(axis, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { refreshRange(); });
  connect(automatic, &QCheckBox::toggled, this, [this](bool on) {
    auto value = this->range(axis->currentIndex()); value.automatic = on;
    this->applyRange(axis->currentIndex(), value); refreshRange();
  });

  colors = new QWidget(owner); layout = grid(colors);
  colorTarget = new QComboBox(colors); colorTarget->setObjectName("viewColorTarget");
  colorTarget->addItems({tr("XX"), tr("XY"), tr("YX"), tr("YY"), tr("Masked Z / PT"),
                         tr("Real tipper"), tr("Imaginary tipper"), tr("Masked tipper")});
  colorTarget->setAccessibleName(tr("Color component"));
  color = new QComboBox(colors); color->setObjectName("viewColor"); color->setEditable(true);
  color->setInsertPolicy(QComboBox::NoInsert);
  for(const auto &hex: {"#579bc9", "#1d5fa2", "#b45309", "#d08738", "#111827", "#6b7280",
                       "#9ca3af", "#bfc4cc", "#dc2626", "#15803d", "#7c3aed", "#0891b2"}) {
    QPixmap swatch(16, 16); swatch.fill(QColor(hex)); color->addItem(QIcon(swatch), hex);
  }
  color->setAccessibleName(tr("Component color"));
  color->setFixedWidth(155); colorTarget->setFixedWidth(155);
  layout->addWidget(colorTarget, 0, 0); layout->addWidget(color, 1, 0);
  connect(colorTarget, qOverload<int>(&QComboBox::currentIndexChanged), this, [this] { refreshColor(); });
  const auto commitColor = [this] {
    const QColor selected(color->currentText().trimmed());
    valid(color, selected.isValid(), tr("Choose a color or type a name / #RRGGBB value and press Enter."));
    if(!selected.isValid()) return;
    appearance.colors[colorTarget->currentIndex()] = selected.name().toStdString();
    this->appearanceChanged();
  };
  connect(color, qOverload<int>(&QComboBox::activated), this, commitColor);
  connect(color->lineEdit(), &QLineEdit::editingFinished, this, commitColor);

  sizes = new QWidget(owner); layout = grid(sizes);
  const QStringList labels{tr("Markers"), tr("Lines"), tr("Arrows"), tr("Ellipses")};
  const std::array<const char *, 4> names{{"viewMarkerSize", "viewLineWidth", "viewArrowScale", "viewEllipseSize"}};
  const std::array<double, 4> minima{{1., .5, 10., 6.}}, maxima{{20., 5., 200., 80.}};
  const std::array<QString, 4> help{{tr("Observation marker diameter in pixels"), tr("Response, error-bar and arrow line width in pixels"),
    tr("Arrow length in pixels for |T| = 1; the reference arrow uses the same scale"), tr("Phase-tensor major diameter in pixels")}};
  for(unsigned i = 0; i < 4; ++i) {
    auto *spin = new QDoubleSpinBox(sizes); sizeControls[i] = spin;
    spin->setObjectName(names[i]); spin->setRange(minima[i], maxima[i]); spin->setDecimals(i < 2 ? 1 : 0);
    spin->setSingleStep(i < 2 ? .5 : 2.); spin->setKeyboardTracking(false); spin->setFixedWidth(62);
    spin->setToolTip(help[i]); spin->setAccessibleName(labels[i]);
    auto *label = new QLabel(labels[i], sizes); label->setBuddy(spin);
    layout->addWidget(label, i % 2, i / 2 * 2); layout->addWidget(spin, i % 2, i / 2 * 2 + 1);
    connect(spin, qOverload<double>(&QDoubleSpinBox::valueChanged), this, [this, i](double value) {
      const std::array<double *, 4> values{{&appearance.markerSize, &appearance.lineWidth, &appearance.arrowScale, &appearance.ellipseSize}};
      *values[i] = value; this->appearanceChanged();
    });
  }

  tensor = new QWidget(owner); layout = grid(tensor);
  colorMap = new QComboBox(tensor); colorMap->setObjectName("viewColorMap"); colorMap->addItems(PlotColorMaps::names());
  colorMap->setAccessibleName(tr("Phase tensor colormap"));
  reverse = new QCheckBox(tr("Reverse"), tensor); reverse->setObjectName("viewColorMapReverse");
  layout->addWidget(colorMap, 0, 0); layout->addWidget(reverse, 0, 1);
  colorMin = new QDoubleSpinBox(tensor); colorMax = new QDoubleSpinBox(tensor);
  colorMin->setObjectName("viewColorMin"); colorMax->setObjectName("viewColorMax");
  colorMin->setPrefix(tr("Min ")); colorMax->setPrefix(tr("Max "));
  for(auto *spin: {colorMin, colorMax}) {
    spin->setRange(-360, 360); spin->setDecimals(1); spin->setSuffix("°"); spin->setFixedWidth(105);
    spin->setKeyboardTracking(false); spin->setToolTip(tr("Minimum principal phase color limits; values outside use the end colors."));
  }
  colorMin->setAccessibleName(tr("Colormap minimum phase")); colorMax->setAccessibleName(tr("Colormap maximum phase"));
  layout->addWidget(colorMin, 1, 0); layout->addWidget(colorMax, 1, 1);
  const auto commitMap = [this] {
    const bool ok = colorMin->value() < colorMax->value();
    for(auto *spin: {colorMin, colorMax}) valid(spin, ok, tr("Minimum principal phase color limits: minimum must be smaller than maximum."));
    if(!ok) return;
    appearance.colorMap = colorMap->currentText().toStdString(); appearance.reverseColorMap = reverse->isChecked();
    appearance.colorMin = colorMin->value(); appearance.colorMax = colorMax->value();
    this->appearanceChanged();
  };
  connect(colorMap, qOverload<int>(&QComboBox::currentIndexChanged), this, commitMap);
  connect(reverse, &QCheckBox::toggled, this, commitMap);
  for(auto *spin: {colorMin, colorMax}) connect(spin, qOverload<double>(&QDoubleSpinBox::valueChanged), this, commitMap);
  setStyle(appearance); refreshRange();
}

std::vector<CommandRibbon::Group> TFViewControls::groups() const
{
  return {{tr("Axis ranges"), {}, ranges}, {tr("Colors"), {}, colors},
          {tr("Sizes (pixels)"), {}, sizes}, {tr("Phase tensor · φmin"), {}, tensor}};
}

void TFViewControls::refreshRange()
{
  const auto value = range(axis->currentIndex());
  const QSignalBlocker block(automatic); automatic->setChecked(value.automatic);
  // Do not replace a partially typed range when another plot is repainted.
  for(const auto &entry: {std::make_pair(lower, value.lower), std::make_pair(upper, value.upper)})
    if(!entry.first->hasFocus()) {
      entry.first->setText(QString::number(entry.second, 'g', 10));
      valid(entry.first, true, tr("Enter a range, then press Enter or Apply. Scientific notation is accepted."));
    }
}

void TFViewControls::refreshColor()
{
  const QSignalBlocker block(color);
  color->setCurrentText(QString::fromStdString(appearance.colors[colorTarget->currentIndex()]));
  valid(color, true, tr("Choose a color or type a name / #RRGGBB value and press Enter."));
}

void TFViewControls::setStyle(TFViewStyle style)
{
  const TFViewStyle defaults;
  for(unsigned i = 0; i < style.colors.size(); ++i)
    if(!QColor(QString::fromStdString(style.colors[i])).isValid()) style.colors[i] = defaults.colors[i];
  if(!PlotColorMaps::names().contains(QString::fromStdString(style.colorMap))) style.colorMap = defaults.colorMap;
  style.markerSize = bounded(style.markerSize, 1, 20, defaults.markerSize);
  style.lineWidth = bounded(style.lineWidth, .5, 5, defaults.lineWidth);
  style.arrowScale = bounded(style.arrowScale, 10, 200, defaults.arrowScale);
  style.ellipseSize = bounded(style.ellipseSize, 6, 80, defaults.ellipseSize);
  style.colorMin = bounded(style.colorMin, -360, 360, defaults.colorMin);
  style.colorMax = bounded(style.colorMax, -360, 360, defaults.colorMax);
  if(style.colorMin >= style.colorMax) { style.colorMin = defaults.colorMin; style.colorMax = defaults.colorMax; }
  appearance = style;
  const std::array<double, 4> values{{style.markerSize, style.lineWidth, style.arrowScale, style.ellipseSize}};
  for(unsigned i = 0; i < 4; ++i) { const QSignalBlocker block(sizeControls[i]); sizeControls[i]->setValue(values[i]); }
  const QSignalBlocker mapBlock(colorMap), reverseBlock(reverse), minBlock(colorMin), maxBlock(colorMax);
  colorMap->setCurrentText(QString::fromStdString(style.colorMap)); reverse->setChecked(style.reverseColorMap);
  colorMin->setValue(style.colorMin); colorMax->setValue(style.colorMax);
  for(auto *spin: {colorMin, colorMax}) valid(spin, true, tr("Minimum principal phase color limits."));
  refreshColor();
}
