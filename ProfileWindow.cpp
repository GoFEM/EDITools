#include "ProfileWindow.h"
#include "include/SharedMapInteraction.h"
#include <QJsonArray>
#include "ProfilePlots.h"
#include "include/FileLabels.h"
#include "include/ImpedanceAnalysis.h"
#include "include/MapBackground.h"
#include "include/PlotColorMaps.h"
#include "include/ToolTips.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMenu>
#include <QMessageBox>
#include <QPdfWriter>
#include <QPushButton>
#include <QSaveFile>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <set>

namespace {
using Point = ProfileAnalysis::Point;
using Quantity = ImpedanceAnalysis::Quantity;
class Corridor : public QCPLayerable {
public:
  explicit Corridor(QCustomPlot *p) : QCPLayerable(p) { setLayer("main"); }
  std::vector<Point> line;
  double distance = 0.;

protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(true); }
  void draw(QCPPainter *p) override {
    if (line.size() < 2)
      return;
    QPainterPath path;
    for (unsigned i = 0; i < line.size(); ++i) {
      const QPointF point(mParentPlot->xAxis->coordToPixel(line[i][0]),
                          mParentPlot->yAxis->coordToPixel(line[i][1]));
      if (i)
        path.lineTo(point);
      else
        path.moveTo(point);
    }
    const double pixels =
        std::abs(mParentPlot->xAxis->coordToPixel(distance) - mParentPlot->xAxis->coordToPixel(0.));
    p->setBrush(Qt::NoBrush);
    p->setPen(QPen(QColor(37, 99, 235, 35), 2 * pixels, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    p->drawPath(path);
  }
};
// Overlay only: it never owns or clears survey-map graphs, labels or backgrounds.
class SharedProfileOverlay : public QCPLayerable {
public:
  explicit SharedProfileOverlay(QCustomPlot *map) : QCPLayerable(map) {
    setObjectName("profileMapOverlay"); setLayer("overlay"); setVisible(false);
  }
  QPolygonF line, stations, highlighted;
  QList<QPolygonF> corridor;
protected:
  QRect clipRect() const override { return mParentPlot->axisRect()->rect(); }
  void applyDefaultAntialiasingHint(QCPPainter *p) const override { p->setAntialiasing(true); }
  QPolygonF pixels(const QPolygonF &points) const {
    QPolygonF result;
    for(const auto &point: points) result << QPointF(mParentPlot->xAxis->coordToPixel(point.x()), mParentPlot->yAxis->coordToPixel(point.y()));
    return result;
  }
  void draw(QCPPainter *p) override {
    QPainterPath area; area.setFillRule(Qt::WindingFill);
    for(const auto &polygon: corridor) area.addPolygon(pixels(polygon));
    p->setPen(Qt::NoPen); p->setBrush(QColor(37, 99, 235, 35)); p->drawPath(area);
    p->setPen(QPen(QColor("#c2410c"), 1.)); p->setBrush(QColor("#c2410c"));
    for(const auto &point: pixels(stations)) p->drawEllipse(point, 3.5, 3.5);
    p->setPen(QPen(QColor("#2563eb"), 1.5)); p->setBrush(Qt::NoBrush);
    const auto points = pixels(line); p->drawPolyline(points);
    for(const auto &point: points) p->drawEllipse(point, 4., 4.);
    p->setPen(QPen(QColor("#d81b60"), 2.)); p->setBrush(Qt::white);
    for(const auto &point: pixels(highlighted)) p->drawEllipse(point, 6., 6.);
  }
};
QDoubleSpinBox *spin(QWidget *owner, const char *name, double low, double high, double value,
                     int decimals = 2) {
  auto *s = new QDoubleSpinBox(owner);
  s->setObjectName(name);
  s->setDecimals(decimals);
  s->setRange(low, high);
  s->setValue(value);
  s->setKeyboardTracking(false);
  return s;
}
QComboBox *combo(QWidget *owner, const char *name, const QStringList &items) {
  auto *c = new QComboBox(owner);
  c->setObjectName(name);
  c->addItems(items);
  return c;
}
void label(QHBoxLayout *row, QWidget *control, const QString &text, const QString &help = {}) {
  row->addWidget(ToolTips::label(control, text, help));
  row->addWidget(control);
}
QString quoted(QString value) { return '"' + value.replace('"', "\"\"") + '"'; }
} // namespace

struct ProfileWindow::Impl {
  ProfileWindow *owner;
  std::function<void(std::shared_ptr<MTSurveyData>)> extracted;
  std::shared_ptr<MTSurveyData> survey;
  std::map<QString, const MTSurveyData *> responses;
  SurveyCoordinates coordinates;
  std::vector<ProfileAnalysis::Station> stations;
  std::vector<ProfileAnalysis::Offset> offsets;
  std::vector<Point> vertices;
  bool automaticDistance = true, updatingPlots = false;
  double stationSpacing = 1.;
  QCustomPlot *map;
  QPointer<QCustomPlot> sharedMap;
  QPointer<SharedProfileOverlay> sharedOverlay;
  SurveyCoordinates sharedCoordinates;
  SharedMapInteraction sharedInteraction;
  MapBackground *background;
  Corridor *corridor;
  QCPCurve *line;
  QComboBox *dataset, *vertical, *periodMin, *periodMax, *palette, *quantities, *tensorColor, *roseBand,
      *bins;
  QCheckBox *drawing, *names, *cells, *equalStations, *filterTensor;
  std::array<QCheckBox *, 3> directions;
  QDoubleSpinBox *distance, *gap, *rhoMin, *rhoMax, *phaseMin, *phaseMax, *tensorMin, *tensorMax, *glyphSize,
      *skew, *ellipticity;
  QTableWidget *table;
  QLabel *status, *coordinateLabel;
  QPushButton *extractButton, *csvButton, *pdfButton;
  std::array<ProfileSectionPlot *, 4> sections;
  ProfileSectionPlot *tensors;
  ProfileRosePlot *roses;
  QPoint press;
  QString selectedStation;
  QCheckBox *rotate;
  QComboBox *strikeMode;
  QDoubleSpinBox *strikeAngle;
  QPushButton *swapStrike;
  QLabel *strikeInfo;
  ProfileAnalysis::Strike inferredStrike;
  bool swappedStrike = false, manualStrikeSet = false;
  double manualStrike = 0., rotation = 0.;

  Impl(ProfileWindow *w, std::function<void(std::shared_ptr<MTSurveyData>)> callback, QCustomPlot *shared)
      : owner(w), extracted(std::move(callback)), sharedMap(shared), sharedInteraction(shared) {
    auto *layout = new QVBoxLayout(w);
    auto *top = new QHBoxLayout;
    dataset = combo(w, "profileDataset", {});
    label(top, dataset, ProfileWindow::tr("Data:"));
    vertical =
        combo(w, "profileVertical", {ProfileWindow::tr("Period (s)"), ProfileWindow::tr("Frequency (Hz)")});
    label(top, vertical, ProfileWindow::tr("Vertical:"));
    periodMin = combo(w, "profilePeriodMin", {});
    periodMax = combo(w, "profilePeriodMax", {});
    label(top, periodMin, ProfileWindow::tr("Periods (s):"),
          ProfileWindow::tr("Measured period range used by sections and rose diagrams. Frequency mode "
                            "changes the vertical axis, not the underlying selection."));
    top->addWidget(periodMax);
    palette = combo(w, "profilePalette", PlotColorMaps::names());

    pdfButton = new QPushButton(ProfileWindow::tr("Save PDF…"), w);
    pdfButton->setObjectName("profileSavePdf");
    top->addWidget(pdfButton);
    layout->addLayout(top);
    auto *rotationRow = new QHBoxLayout;
    rotate = new QCheckBox(ProfileWindow::tr("Rotate to strike"), w);
    rotate->setObjectName("profileRotateStrike");
    rotate->setToolTip(ProfileWindow::tr(
        "Rotate the displayed data axes clockwise from geographic north, assuming stored x/y channels use "
        "north/east axes. The drawn line, station "
        "coordinates, original data, errors and masks stay unchanged. XY uses electric field parallel to the "
        "chosen strike axis. Phase-tensor strike has a 90° ambiguity; +90° swaps the axes."));
    rotationRow->addWidget(rotate);
    strikeMode =
        combo(w, "profileStrikeMode", {ProfileWindow::tr("Auto (observed)"), ProfileWindow::tr("Manual")});
    rotationRow->addWidget(strikeMode);
    strikeAngle = spin(w, "profileStrikeAngle", -181., 180., 0., 2);
    strikeAngle->setSuffix("°");
    label(rotationRow, strikeAngle, ProfileWindow::tr("Strike:"),
          ProfileWindow::tr("Clockwise from north. Auto uses the selected corridor's observed phase tensors, "
                            "period range and rose band, with the skew/ellipticity thresholds and station "
                            "weighting below. Manual accepts a chosen angle."));
    swapStrike = new QPushButton("+90°", w);
    swapStrike->setObjectName("profileSwapStrike");
    rotationRow->addWidget(swapStrike);
    strikeInfo = new QLabel(w);
    strikeInfo->setObjectName("profileStrikeInfo");
    rotationRow->addWidget(strikeInfo);
    rotationRow->addStretch();
    layout->addLayout(rotationRow);
    auto *split = new QSplitter(w);
    layout->addWidget(split, 1);
    auto *left = new QWidget(split);
    auto *ml = new QVBoxLayout(left);
    ml->setContentsMargins(0, 0, 0, 0);
    auto *tools = new QHBoxLayout;
    drawing = new QCheckBox(ProfileWindow::tr("Draw line"), w);
    drawing->setObjectName("profileDraw");
    drawing->setChecked(true);
    drawing->setToolTip(
        ProfileWindow::tr("Click successive vertices on the map. Right-click to finish. Turn off Draw line "
                          "to pan. Undo removes the last vertex; Reverse changes the distance origin."));
    tools->addWidget(drawing);
    auto button = [&](const QString &text, const char *name, std::function<void()> action) {
      auto *b = new QPushButton(text, w);
      b->setObjectName(name);
      tools->addWidget(b);
      QObject::connect(b, &QPushButton::clicked, w, std::move(action));
      return b;
    };
    button(ProfileWindow::tr("Undo vertex"), "profileUndo", [this] {
      if (!vertices.empty())
        vertices.pop_back();
      select();
    });
    button(ProfileWindow::tr("Reverse"), "profileReverse", [this] {
      std::reverse(vertices.begin(), vertices.end());
      select();
    });
    button(ProfileWindow::tr("Clear"), "profileClear", [this] {
      vertices.clear();
      select();
    });
    button(ProfileWindow::tr("Fit"), "profileFitMap", [this] { fitMap(); });
    ml->addLayout(tools);
    auto *widthRow = new QHBoxLayout;
    distance = spin(w, "profileDistance", 0., 100000., 1., 3);
    label(widthRow, distance, ProfileWindow::tr("Distance ± (km):"),
          ProfileWindow::tr("Maximum distance from the finite polyline, on each side and around its "
                            "endpoints. Default: median nonzero nearest-neighbour spacing of enabled survey "
                            "stations, in the displayed UTM projection. Wide surveys inherit that "
                            "projection's distortion. No station is included twice at bends or crossings."));
    auto *automatic = new QPushButton(ProfileWindow::tr("Auto"), w);
    automatic->setObjectName("profileAutoDistance");
    widthRow->addWidget(automatic);
    ml->addLayout(widthRow);
    auto *mapTools = new QHBoxLayout;
    names = new QCheckBox(ProfileWindow::tr("Station names"), w);
    names->setObjectName("profileNames");
    mapTools->addWidget(names);
    map = new QCustomPlot(left);
    map->setObjectName(sharedMap ? "profileExportMap" : "profileMap");
    map->setMinimumSize(320, 260);
    background = new MapBackground(map);
    corridor = new Corridor(map);
    mapTools->addWidget(MapBackground::button(left, background->layers(), [this](int layers) {
      background->setLayers(layers);
      map->replot();
    }));
    mapTools->addStretch();
    if(!sharedMap) ml->addLayout(mapTools);
    for (int i = 0; i < 3; ++i) {
      auto *g = map->addGraph();
      g->setLineStyle(QCPGraph::lsNone);
      g->setSelectable(QCP::stNone);
      g->setScatterStyle(QCPScatterStyle(i == 2 ? QCPScatterStyle::ssCircle : QCPScatterStyle::ssDisc,
                                         QColor(i == 1   ? "#c2410c"
                                                : i == 2 ? "#2563eb"
                                                         : "#808080"),
                                         i == 0 ? 4 : 7));
    }
    line = new QCPCurve(map->xAxis, map->yAxis);
    line->setPen(QPen(QColor("#2563eb"), 1.4));
    line->setSelectable(QCP::stNone);
    line->removeFromLegend();
    auto *highlighted = map->addGraph();
    highlighted->setObjectName("profileHighlightedStation");
    highlighted->setLineStyle(QCPGraph::lsNone);
    highlighted->setSelectable(QCP::stNone);
    highlighted->setScatterStyle(
        QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor("#d81b60"), 2.), QBrush(Qt::white), 12));
    highlighted->setLayer("overlay");
    map->xAxis->setLabel(ProfileWindow::tr("Easting (km)"));
    map->yAxis->setLabel(ProfileWindow::tr("Northing (km)"));
    map->setInteractions(QCP::iRangeZoom);
    map->setCursor(Qt::CrossCursor);
    if(sharedMap) map->hide(); else ml->addWidget(map, 3);
    coordinateLabel = new QLabel(w);
    coordinateLabel->setWordWrap(true);
    ml->addWidget(coordinateLabel);
    table = new QTableWidget(w);
    table->setObjectName("profileStations");
    table->setColumnCount(3);
    table->setHorizontalHeaderLabels(
        {ProfileWindow::tr("Station"), ProfileWindow::tr("Along (km)"), ProfileWindow::tr("Offset (km)")});
    table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::SingleSelection);
    table->verticalHeader()->hide();
    ml->addWidget(table, 1);
    auto *outputs = new QHBoxLayout;
    extractButton = new QPushButton(ProfileWindow::tr("Open as survey"), w);
    extractButton->setObjectName("profileExtract");
    extractButton->setToolTip(
        ProfileWindow::tr("Open a copy of the corridor stations with their original data and masks. Display "
                          "rotation is not applied to the copied survey."));
    outputs->addWidget(extractButton);
    csvButton = new QPushButton(ProfileWindow::tr("Station CSV…"), w);
    csvButton->setObjectName("profileCsv");
    outputs->addWidget(csvButton);
    ml->addLayout(outputs);
    auto *tabs = new QTabWidget(split);
    tabs->setObjectName("profileTabs");
    auto *scalarPage = new QWidget(tabs);
    auto *sl = new QVBoxLayout(scalarPage);
    auto *scalarTools = new QHBoxLayout;
    quantities = combo(w, "profileQuantities",
                       {ProfileWindow::tr("XY / YX"), ProfileWindow::tr("Berdichevsky / SSQ"),
                        ProfileWindow::tr("Determinant")});
    label(scalarTools, quantities, ProfileWindow::tr("Response:"));
    quantities->setToolTip(ProfileWindow::tr(
        "Uses the same complex impedance quantities and masks as impedance analysis. YX is sign-reversed. "
        "The "
        "Rotate to strike option changes the data axes; otherwise XY/YX remain the stored components."));
    cells = new QCheckBox(ProfileWindow::tr("Cells"), w);
    cells->setObjectName("profileCells");
    cells->setChecked(true);
    scalarTools->addWidget(cells);
    cells->setToolTip(
        ProfileWindow::tr("Color each measured sample's bounded cell. No interpolation. Turn off for sample "
                          "markers. Masked values and large station/period gaps stay blank. Stations with "
                          "identical along-profile offsets overlap; narrow the corridor to separate them."));
    gap = spin(w, "profileGap", 1.01, 100., 2.);
    label(scalarTools, gap, ProfileWindow::tr("Gap ratio:"),
          ProfileWindow::tr("Maximum period ratio spanned by one sample cell. Larger gaps remain blank; this "
                            "never fills missing or masked data."));
    label(scalarTools, palette, ProfileWindow::tr("Colormap:"));
    sl->addLayout(scalarTools);
    auto *ranges = new QHBoxLayout;
    rhoMin = spin(w, "profileRhoMin", .0001, 1e9 - 1, 1., 4);
    rhoMax = spin(w, "profileRhoMax", .0002, 1e9, 10000., 4);
    phaseMin = spin(w, "profilePhaseMin", -180., 179.99, 0.);
    phaseMax = spin(w, "profilePhaseMax", -179.99, 180., 90.);
    label(ranges, rhoMin, ProfileWindow::tr("ρa:"), ProfileWindow::tr("Logarithmic color limits in Ω m."));
    ranges->addWidget(rhoMax);
    label(ranges, phaseMin, ProfileWindow::tr("Phase:"),
          ProfileWindow::tr("Linear color limits in degrees; out-of-range values use endpoint colors."));
    ranges->addWidget(phaseMax);
    sl->addLayout(ranges);
    auto *grid = new QGridLayout;
    sl->addLayout(grid, 1);
    for (unsigned i = 0; i < 4; ++i) {
      sections[i] = new ProfileSectionPlot(w);
      sections[i]->setObjectName(QString("profileSection%1").arg(i));
      grid->addWidget(sections[i], i / 2, i % 2);
    }
    tabs->addTab(scalarPage, ProfileWindow::tr("Resistivity / phase"));
    auto *ptPage = new QWidget(tabs);
    auto *pl = new QVBoxLayout(ptPage);
    auto *ptTools = new QHBoxLayout;
    tensorColor = combo(
        w, "profileTensorColor",
        {ProfileWindow::tr("φmax (°)"), ProfileWindow::tr("φmin (°)"), ProfileWindow::tr("Skew β (°)")});
    label(ptTools, tensorColor, ProfileWindow::tr("Color:"));
    tensorMin = spin(w, "profileTensorMin", -180., 179.99, 0.);
    tensorMax = spin(w, "profileTensorMax", -179.99, 180., 90.);
    ptTools->addWidget(tensorMin);
    ptTools->addWidget(tensorMax);
    glyphSize = spin(w, "profileGlyphSize", 3., 80., 14., 0);
    label(ptTools, glyphSize, ProfileWindow::tr("Size:"),
          ProfileWindow::tr("Ellipse major diameter in pixels. Geometry is independent of axis aspect ratio. "
                            "Up within each ellipse is the selected data x axis (north unless rotated); "
                            "right is +90°. Vertical "
                            "placement is period/frequency, not depth."));
    auto *tensorPalette = combo(w, "profileTensorPalette", PlotColorMaps::names());
    label(ptTools, tensorPalette, ProfileWindow::tr("Colormap:"));
    QObject::connect(tensorPalette, QOverload<int>::of(&QComboBox::currentIndexChanged), palette,
                     &QComboBox::setCurrentIndex);
    QObject::connect(palette, QOverload<int>::of(&QComboBox::currentIndexChanged), tensorPalette,
                     &QComboBox::setCurrentIndex);
    pl->addLayout(ptTools);
    auto *tensorSplit = new QSplitter(Qt::Vertical, ptPage);
    tensorSplit->setObjectName("profileTensorRoseSplit");
    pl->addWidget(tensorSplit, 1);
    tensors = new ProfileSectionPlot(w);
    tensors->setObjectName("profileTensorPlot");
    tensorSplit->addWidget(tensors);
    tabs->addTab(ptPage, ProfileWindow::tr("Phase tensors / strike"));
    auto *rosePage = new QWidget(tensorSplit);
    auto *rl = new QVBoxLayout(rosePage);
    auto *roseTools = new QHBoxLayout;
    const QStringList titles{ProfileWindow::tr("PT azimuth"), ProfileWindow::tr("Real tipper normal"),
                             ProfileWindow::tr("Imaginary tipper normal")};
    for (unsigned i = 0; i < 3; ++i) {
      directions[i] = new QCheckBox(titles[i], w);
      directions[i]->setObjectName(QString("profileDirection%1").arg(i));
      directions[i]->setChecked(i < 2);
      roseTools->addWidget(directions[i]);
    }
    equalStations = new QCheckBox(ProfileWindow::tr("Equal stations"), w);
    equalStations->setObjectName("profileEqualStations");
    equalStations->setChecked(true);
    roseTools->addWidget(equalStations);
    equalStations->setToolTip(ProfileWindow::tr(
        "Give each contributing station total weight one in roses and automatic strike, so denser period "
        "sampling does not dominate. "
        "Otherwise every valid station-period sample has equal weight. Roses are axial (0–180°), mirrored on "
        "the opposite side; PT strike has an additional 90° ambiguity. Tipper normals are perpendicular to "
        "the induction direction, not an independently determined geological strike."));
    rl->addLayout(roseTools);
    auto *filterTools = new QHBoxLayout;
    roseBand = combo(w, "profileRoseBand", {});
    label(filterTools, roseBand, ProfileWindow::tr("Band:"));
    bins = combo(w, "profileRoseBins", {"5°", "10°", "15°", "30°"});
    bins->setCurrentIndex(1);
    label(filterTools, bins, ProfileWindow::tr("Bins:"));
    filterTensor = new QCheckBox(ProfileWindow::tr("Filter PT"), w);
    filterTensor->setObjectName("profileFilterTensor");
    filterTensor->setChecked(true);
    filterTools->addWidget(filterTensor);
    skew = spin(w, "profileMaxSkew", 0., 90., 3.);
    ellipticity = spin(w, "profileMinEllipticity", 0., 1., .1);
    label(filterTools, skew, ProfileWindow::tr("|β| ≤"),
          ProfileWindow::tr("Auto strike and filtered PT rose: maximum absolute skew, in degrees. Does not "
                            "change masks or section visibility."));
    label(filterTools, ellipticity, ProfileWindow::tr("Ellipticity ≥"),
          ProfileWindow::tr(
              "Auto strike and filtered PT rose: (1 − axis ratio)/(1 + axis ratio). Nearly circular tensors "
              "have poorly defined "
              "azimuths. Exactly circular tensors and zero tipper vectors are always excluded from roses."));
    rl->addLayout(filterTools);
    roses = new ProfileRosePlot(w);
    roses->setObjectName("profileRosePlot");
    rl->addWidget(roses, 1);
    tensorSplit->addWidget(rosePage);
    tensorSplit->setSizes({450, 350});
    std::vector<ProfileSectionPlot *> linked(sections.begin(), sections.end());
    linked.push_back(tensors);
    for (auto *plot : linked)
      QObject::connect(plot, &ProfileSectionPlot::stationClicked, w,
                       [this](const QString &name) { highlight(name); });
    QObject::connect(table, &QTableWidget::currentCellChanged, w, [this](int row, int, int, int) {
      highlight(row >= 0 && table->item(row, 0) ? table->item(row, 0)->text() : QString());
    });
    for (auto *plot : linked)
      for (bool horizontal : {true, false}) {
        auto *axis = horizontal ? plot->xAxis : plot->yAxis;
        QObject::connect(axis, QOverload<const QCPRange &>::of(&QCPAxis::rangeChanged), w,
                         [this, linked, plot, horizontal](const QCPRange &range) {
                           if (updatingPlots)
                             return;
                           updatingPlots = true;
                           for (auto *other : linked)
                             if (other != plot) {
                               (horizontal ? other->xAxis : other->yAxis)->setRange(range);
                               other->replot(QCustomPlot::rpQueuedReplot);
                             }
                           updatingPlots = false;
                         });
      }
    split->setSizes(sharedMap ? QList<int>{260, 1160} : QList<int>{420, 1000});
    status = new QLabel(w);
    status->setObjectName("profileStatus");
    status->setWordWrap(true);
    layout->addWidget(status);

    QObject::connect(rotate, &QCheckBox::toggled, w, [this] { plots(false); });
    QObject::connect(strikeMode, QOverload<int>::of(&QComboBox::currentIndexChanged), w, [this] {
      if (strikeMode->currentIndex() == 1 && !manualStrikeSet) {
        manualStrike =
            std::isfinite(inferredStrike.degrees) ? inferredStrike.degrees + (swappedStrike ? 90. : 0.) : 0.;
        manualStrikeSet = true;
      }
      plots(false);
    });
    QObject::connect(strikeAngle, QOverload<double>::of(&QDoubleSpinBox::valueChanged), w,
                     [this](double angle) {
                       if (strikeMode->currentIndex() == 1) {
                         manualStrike = angle;
                         manualStrikeSet = true;
                         plots(false);
                       }
                     });
    QObject::connect(swapStrike, &QPushButton::clicked, w, [this] {
      if (strikeMode->currentIndex() == 0)
        swappedStrike = !swappedStrike;
      else
        manualStrike = std::remainder(manualStrike + 90., 180.);
      plots(false);
    });
    QObject::connect(drawing, &QCheckBox::toggled, w, [this](bool on) {
      map->setInteraction(QCP::iRangeDrag, !on);
      map->setCursor(on ? Qt::CrossCursor : Qt::ArrowCursor);
      sharedInteractions();
    });
    QObject::connect(map, &QCustomPlot::mousePress, w,
                     [this](QMouseEvent *e) { press = e->position().toPoint(); });
    QObject::connect(map, &QCustomPlot::mouseRelease, w, [this](QMouseEvent *e) {
      if (!drawing->isChecked() || !map->axisRect()->rect().contains(e->position().toPoint()))
        return;
      if (e->button() == Qt::RightButton) {
        drawing->setChecked(false);
        return;
      }
      if (e->button() != Qt::LeftButton || (e->position().toPoint() - press).manhattanLength() > 4)
        return;
      const Point point{
          {map->xAxis->pixelToCoord(e->position().x()), map->yAxis->pixelToCoord(e->position().y())}};
      if (vertices.empty() || std::hypot(point[0] - vertices.back()[0], point[1] - vertices.back()[1]) > 1e-6)
        vertices.push_back(point);
      select();
    });
    QObject::connect(map, &QCustomPlot::afterLayout, w, [this] {
      const double width = map->axisRect()->width(), height = map->axisRect()->height();
      if (width <= 0 || height <= 0)
        return;
      if (map->xAxis->range().size() / width > map->yAxis->range().size() / height)
        map->yAxis->setScaleRatio(map->xAxis);
      else
        map->xAxis->setScaleRatio(map->yAxis);
    });
    QObject::connect(distance, QOverload<double>::of(&QDoubleSpinBox::valueChanged), w, [this] {
      automaticDistance = false;
      select();
    });
    QObject::connect(automatic, &QPushButton::clicked, w, [this] {
      automaticDistance = true;
      defaultDistance();
      select();
    });
    QObject::connect(names, &QCheckBox::toggled, w, [this] { updateMap(); });
    QObject::connect(dataset, QOverload<int>::of(&QComboBox::currentIndexChanged), w, [this] {
      periods();
      plots(true);
    });
    QObject::connect(vertical, QOverload<int>::of(&QComboBox::currentIndexChanged), w,
                     [this] { plots(true); });
    for (auto *c : {palette, quantities, tensorColor, roseBand, bins})
      QObject::connect(c, QOverload<int>::of(&QComboBox::currentIndexChanged), w, [this, c] {
        if (c == tensorColor) {
          QSignalBlocker a(tensorMin), b(tensorMax);
          tensorMin->setValue(c->currentIndex() == 2 ? -10. : c->currentIndex() == 1 ? -90. : 0.);
          tensorMax->setValue(c->currentIndex() == 2 ? 10. : 90.);
        }
        plots(false);
      });
    for (auto *c : {periodMin, periodMax})
      QObject::connect(c, QOverload<int>::of(&QComboBox::currentIndexChanged), w, [this, c] {
        if (periodMin->currentIndex() > periodMax->currentIndex()) {
          auto *other = c == periodMin ? periodMax : periodMin;
          QSignalBlocker b(other);
          other->setCurrentIndex(c->currentIndex());
        }
        bands();
        plots(true);
      });
    for (auto *c : {cells, equalStations, filterTensor, directions[0], directions[1], directions[2]})
      QObject::connect(c, &QCheckBox::toggled, w, [this] { plots(false); });
    for (auto *s :
         {gap, rhoMin, rhoMax, phaseMin, phaseMax, tensorMin, tensorMax, glyphSize, skew, ellipticity})
      QObject::connect(s, QOverload<double>::of(&QDoubleSpinBox::valueChanged), w, [this, s] {
        for (auto range : {std::make_pair(rhoMin, rhoMax), std::make_pair(phaseMin, phaseMax),
                           std::make_pair(tensorMin, tensorMax)})
          if (range.first->value() >= range.second->value()) {
            auto *other = s == range.first ? range.second : range.first;
            QSignalBlocker block(other);
            other->setValue(s->value() + (s == range.first ? 1. : -1.) * std::pow(10., -other->decimals()));
          }
        plots(false);
      });
    QObject::connect(extractButton, &QPushButton::clicked, w, [this] {
      if (!survey || offsets.empty())
        return;
      auto copy = std::make_shared<MTSurveyData>(*survey);
      std::set<std::string> keep;
      for (const auto &s : offsets)
        keep.insert(s.name);
      for (const auto &name : survey->get_stations_names())
        if (!keep.count(name))
          copy->remove_station(name);
      if (extracted)
        extracted(std::move(copy));
    });
    QObject::connect(csvButton, &QPushButton::clicked, w, [this] { saveCsv(); });
    QObject::connect(pdfButton, &QPushButton::clicked, w, [this] { savePdf(); });
    split->setObjectName("profileMainSplit");
    w->addToolGroup(ProfileWindow::tr("Profile tools"), ProfileWindow::tr("Data and periods"), top);
    w->addToolGroup(ProfileWindow::tr("Profile tools"), ProfileWindow::tr("Line"), tools);
    w->addToolGroup(ProfileWindow::tr("Profile tools"), ProfileWindow::tr("Corridor"), widthRow);
    if(!sharedMap) w->addToolGroup(ProfileWindow::tr("Profile tools"), ProfileWindow::tr("Map"), mapTools);
    else {
      while(auto *item = mapTools->takeAt(0)) { if(item->widget()) item->widget()->hide(); delete item; }
      delete mapTools;
    }
    w->addToolGroup(ProfileWindow::tr("Profile tools"), ProfileWindow::tr("Output"), outputs);
    w->addToolGroup(ProfileWindow::tr("Profile appearance"), ProfileWindow::tr("Scalar sections"), scalarTools);
    w->addToolGroup(ProfileWindow::tr("Profile appearance"), ProfileWindow::tr("Color ranges"), ranges);
    w->addToolGroup(ProfileWindow::tr("Profile appearance"), ProfileWindow::tr("Phase tensors"), ptTools);
    w->addToolGroup(ProfileWindow::tr("Strike tools"), ProfileWindow::tr("Rotation"), rotationRow);
    w->addToolGroup(ProfileWindow::tr("Strike tools"), ProfileWindow::tr("Rose directions"), roseTools);
    w->addToolGroup(ProfileWindow::tr("Strike tools"), ProfileWindow::tr("Quality and weighting"), filterTools);
    if(sharedMap) {
      sharedOverlay = new SharedProfileOverlay(sharedMap);
      QObject::connect(sharedMap, &QCustomPlot::mousePress, w, [this](QMouseEvent *e) {
        if(sharedInteraction.active()) press = e->position().toPoint();
      });
      QObject::connect(sharedMap, &QCustomPlot::mouseRelease, w, [this](QMouseEvent *e) {
        if(!sharedInteraction.active() || !drawing->isChecked() || !sharedMap->axisRect()->rect().contains(e->position().toPoint())) return;
        if(e->button() == Qt::RightButton) { drawing->setChecked(false); return; }
        if(e->button() != Qt::LeftButton || (e->position().toPoint() - press).manhattanLength() > 4) return;
        try {
          const std::array<double, 3> point{{sharedMap->yAxis->pixelToCoord(e->position().y()), sharedMap->xAxis->pixelToCoord(e->position().x()), 0.}};
          const auto projected = coordinates.transform(sharedCoordinates.inverse({point})).front();
          const Point vertex{{projected[1] / 1000., projected[0] / 1000.}};
          if(vertices.empty() || std::hypot(vertex[0] - vertices.back()[0], vertex[1] - vertices.back()[1]) > 1e-6) vertices.push_back(vertex);
          select();
        } catch(const std::exception &e) { status->setText(ProfileWindow::tr("Cannot project profile vertex: %1").arg(e.what())); }
      });
    }
    w->setExportButton("pdf", pdfButton); w->setExportButton("csv", csvButton);
    w->setStateRefresh([this] {
      manualStrike = strikeAngle->value();
      map->setInteraction(QCP::iRangeDrag, !drawing->isChecked());
      map->setCursor(drawing->isChecked() ? Qt::CrossCursor : Qt::ArrowCursor);
      periods(); select(); sharedInteractions();
    });

  }

  ~Impl() {
    deactivateSharedMap();
    delete sharedOverlay.data();
  }
  void deactivateSharedMap() {
    if(sharedOverlay) sharedOverlay->setVisible(false);
    sharedInteraction.deactivate();
  }
  void sharedInteractions() {
    if(drawing->isChecked()) {
      const QCursor cursor(Qt::CrossCursor);
      sharedInteraction.configure(QCP::iRangeDrag | QCP::iSelectPlottables, &cursor);
    } else sharedInteraction.configure();
  }
  QPolygonF toSharedMap(const std::vector<Point> &points) const {
    std::vector<std::array<double, 3>> projected;
    for(const auto &p: points) projected.push_back({{p[1] * 1000., p[0] * 1000., 0.}});
    QPolygonF result;
    for(const auto &p: sharedCoordinates.transform(coordinates.inverse(projected))) result << QPointF(p[1], p[0]);
    return result;
  }
  void updateSharedMap() {
    if(!sharedInteraction.active() || !sharedMap || !sharedOverlay) return;
    try {
      sharedOverlay->line = toSharedMap(vertices);
      sharedOverlay->corridor.clear();
      const double radius = distance->value(), pi = std::acos(-1.);
      for(unsigned i = 1; i < vertices.size() && radius > 0.; ++i) {
        const auto a = vertices[i - 1], b = vertices[i];
        const double angle = std::atan2(b[1] - a[1], b[0] - a[0]);
        std::vector<Point> polygon;
        for(int end = 0; end < 2; ++end) for(int step = 0; step <= 24; ++step) {
          const auto center = end ? a : b;
          const double t = angle - pi / 2. + end * pi + step * pi / 24.;
          polygon.push_back({{center[0] + radius * std::cos(t), center[1] + radius * std::sin(t)}});
        }
        sharedOverlay->corridor << toSharedMap(polygon);
      }
      std::vector<Point> selected, highlighted;
      for(const auto &station: stations) {
        if(std::any_of(offsets.begin(), offsets.end(), [&](const ProfileAnalysis::Offset &p) { return p.name == station.name; })) selected.push_back(station.position);
        if(QString::fromStdString(station.name) == selectedStation) highlighted.push_back(station.position);
      }
      sharedOverlay->stations = toSharedMap(selected); sharedOverlay->highlighted = toSharedMap(highlighted);
      sharedOverlay->setProperty("vertexCount", int(vertices.size()));
      sharedOverlay->setVisible(true);
    } catch(const std::exception &e) {
      sharedOverlay->setVisible(false);
      status->setText(ProfileWindow::tr("Cannot display profile on the survey map: %1").arg(e.what()));
    }
    sharedMap->replot(QCustomPlot::rpQueuedReplot);
  }
  const MTSurveyData *source() const {
    const auto key = dataset->currentData().toString();
    const auto found = responses.find(key);
    return key.isEmpty() ? survey.get() : (found == responses.end() ? nullptr : found->second);
  }
  void defaultDistance() {
    QSignalBlocker block(distance);
    distance->setValue(stationSpacing);
  }
  void periods() {
    const double low = periodMin->currentData().toDouble(), high = periodMax->currentData().toDouble();
    QSignalBlocker a(periodMin), b(periodMax);
    periodMin->clear();
    periodMax->clear();
    if (source())
      for (double p : source()->get_unique_periods())
        if (std::isfinite(p) && p > 0.)
          for (auto *c : {periodMin, periodMax})
            c->addItem(QString::number(p, 'g', 8), p);
    periodMin->setCurrentIndex(std::max(0, periodMin->findData(low)));
    const int index = periodMax->findData(high);
    periodMax->setCurrentIndex(index < 0 ? periodMax->count() - 1 : index);
    bands();
  }
  void bands() {
    const auto old = roseBand->currentText();
    QSignalBlocker block(roseBand);
    roseBand->clear();
    roseBand->addItem(ProfileWindow::tr("All selected periods"));
    const double low = periodMin->currentData().toDouble(), high = periodMax->currentData().toDouble();
    if (low > 0. && high >= low)
      for (int decade = std::max(-20, int(std::floor(std::log10(low))));
           decade <= std::min(20, int(std::floor(std::log10(high)))); ++decade)
        roseBand->addItem(QString("%1–%2 s")
                              .arg(std::pow(10., decade), 0, 'g', 4)
                              .arg(std::pow(10., decade + 1), 0, 'g', 4),
                          decade);
    roseBand->setCurrentIndex(std::max(0, roseBand->findText(old)));
  }
  void highlight(const QString &name, bool notify = true) {
    int row = -1;
    for (unsigned i = 0; i < offsets.size(); ++i)
      if (QString::fromStdString(offsets[i].name) == name) {
        row = int(i);
        break;
      }
    selectedStation = row < 0 ? QString() : name;
    if(notify && !selectedStation.isEmpty()) emit owner->stationActivated(selectedStation);
    const QSignalBlocker block(table);
    table->clearSelection();
    table->setCurrentCell(row, row < 0 ? -1 : 0);
    if (row >= 0) {
      table->selectRow(row);
      table->scrollToItem(table->item(row, 0));
    }
    for (auto *plot : sections)
      plot->setSelectedStation(selectedStation);
    tensors->setSelectedStation(selectedStation);
    updateMap();
  }
  void select() {
    offsets = ProfileAnalysis::select(stations, vertices, distance->value());
    corridor->line = vertices;
    corridor->distance = distance->value();
    const QSignalBlocker block(table);
    table->setRowCount(int(offsets.size()));
    for (unsigned i = 0; i < offsets.size(); ++i) {
      const auto &s = offsets[i];
      table->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(s.name)));
      table->setItem(i, 1, new QTableWidgetItem(QString::number(s.along, 'f', 3)));
      table->setItem(i, 2, new QTableWidgetItem(QString::number(s.across, 'f', 3)));
    }
    extractButton->setEnabled(!offsets.empty());
    csvButton->setEnabled(!offsets.empty());
    pdfButton->setEnabled(!offsets.empty());
    plots(true);
    highlight(selectedStation);
  }
  void updateMap() {
    QVector<double> x, y, sx, sy, vx, vy, order, hx, hy;
    std::set<std::string> selected;
    for (const auto &s : offsets)
      selected.insert(s.name);
    map->clearItems();
    for (const auto &s : stations) {
      x.push_back(s.position[0]);
      y.push_back(s.position[1]);
      if (selected.count(s.name)) {
        sx.push_back(s.position[0]);
        sy.push_back(s.position[1]);
      }
      const bool highlighted = QString::fromStdString(s.name) == selectedStation;
      if (highlighted) {
        hx.push_back(s.position[0]);
        hy.push_back(s.position[1]);
      }
      if (names->isChecked() || highlighted) {
        auto *label = new QCPItemText(map);
        label->setText(QString::fromStdString(s.name));
        label->position->setCoords(s.position[0], s.position[1]);
        label->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom);
        label->setFont(QFont(owner->font().family(), 8));
        label->setSelectable(false);
        if (highlighted)
          label->setColor(QColor("#d81b60"));
      }
    }
    for (unsigned i = 0; i < vertices.size(); ++i) {
      vx.push_back(vertices[i][0]);
      vy.push_back(vertices[i][1]);
      order.push_back(i);
    }
    map->graph(0)->setData(x, y);
    map->graph(1)->setData(sx, sy);
    map->graph(2)->setData(vx, vy);
    map->graph(3)->setData(hx, hy);
    line->setData(order, vx, vy);
    map->replot();
    updateSharedMap();
  }
  void fitMap() {
    if (stations.empty()) {
      map->xAxis->setRange(-1, 1);
      map->yAxis->setRange(-1, 1);
    } else {
      double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
      for (const auto &s : stations) {
        x0 = std::min(x0, s.position[0]);
        x1 = std::max(x1, s.position[0]);
        y0 = std::min(y0, s.position[1]);
        y1 = std::max(y1, s.position[1]);
      }
      const double margin = std::max(.1, .08 * std::max(x1 - x0, y1 - y0));
      map->xAxis->setRange(x0 - margin, x1 + margin);
      map->yAxis->setRange(y0 - margin, y1 + margin);
    }
    map->replot();
    if(sharedInteraction.active() && sharedMap) {
      std::vector<Point> points = vertices;
      for(const auto &station: stations) points.push_back(station.position);
      const auto bounds = toSharedMap(points).boundingRect();
      if(!points.empty()) {
        const double margin = std::max(sharedCoordinates.utm ? 100. : .001, .08 * std::max(bounds.width(), bounds.height()));
        sharedMap->xAxis->setRange(bounds.left() - margin, bounds.right() + margin);
        sharedMap->yAxis->setRange(bounds.top() - margin, bounds.bottom() + margin);
        sharedMap->replot();
      }
    }
  }
  bool inRoseBand(double period) const {
    if (roseBand->currentIndex() <= 0)
      return true;
    const int decade = roseBand->currentData().toInt();
    return period >= std::pow(10., decade) * (1. - 1e-10) &&
           period < std::pow(10., decade + 1) * (1. - 1e-10);
  }
  bool tensorForStrike(const MTMapData::PhaseTensor &tensor) const {
    const double shape = (1. - tensor.axisRatio) / (1. + tensor.axisRatio);
    return shape > 1e-8 && shape >= ellipticity->value() && std::abs(tensor.skew) <= skew->value();
  }
  QString frame() const {
    if (!rotate->isChecked())
      return ProfileWindow::tr("geographic north");
    return std::isfinite(rotation)
               ? ProfileWindow::tr("strike %1° clockwise from north").arg(rotation, 0, 'f', 2)
               : ProfileWindow::tr("strike unavailable");
  }
  void updateStrike(double low, double high) {
    const bool rotating = rotate->isChecked(), automaticMode = strikeMode->currentIndex() == 0;
    for(auto *control : std::vector<QWidget *>{strikeMode, swapStrike, strikeInfo}) owner->setToolVisible(control, rotating);
    owner->setToolVisible(strikeAngle, rotating && !automaticMode);
    const bool quality = filterTensor->isChecked() || (rotating && automaticMode);
    owner->setToolVisible(skew, quality); owner->setToolVisible(ellipticity, quality);
    owner->setToolVisible(gap, cells->isChecked());

    std::vector<ProfileAnalysis::Direction> directions;
    if (survey)
      for (const auto &site : offsets) {
        const auto &station = survey->get_station_data(site.name);
        for (unsigned f = 0; f < station.frequencies().size(); ++f) {
          const double period = 1. / station.frequencies()[f];
          MTMapData::PhaseTensor tensor;
          if (!std::isfinite(period) || period < low * (1. - 1e-10) || period > high * (1. + 1e-10) ||
              !inRoseBand(period))
            continue;
          if (MTMapData::phase_tensor(station, f, tensor) && tensorForStrike(tensor))
            directions.push_back({site.name, tensor.azimuth});
        }
      }
    inferredStrike = ProfileAnalysis::estimateStrike(directions, equalStations->isChecked());
    const bool automatic = strikeMode->currentIndex() == 0;
    const double chosen = automatic ? inferredStrike.degrees + (swappedStrike ? 90. : 0.) : manualStrike;
    {
      QSignalBlocker block(strikeAngle);
      strikeAngle->setMinimum(automatic ? -181. : -180.);
      strikeAngle->setSpecialValueText(automatic ? QStringLiteral("—") : QString());
      strikeAngle->setReadOnly(automatic);
      strikeAngle->setValue(std::isfinite(chosen) ? chosen : -181.);
    }
    swapStrike->setEnabled(std::isfinite(chosen));
    rotation = rotate->isChecked() ? chosen : 0.;
    strikeInfo->setText(automatic
                            ? (!std::isfinite(chosen) ? ProfileWindow::tr("No estimate — choose Manual")
                                                      : ProfileWindow::tr("%1 sites · R₄ %2%3")
                                                            .arg(inferredStrike.stations)
                                                            .arg(inferredStrike.concentration, 0, 'f', 2)
                                                            .arg(inferredStrike.concentration < .25
                                                                     ? ProfileWindow::tr(" · weak agreement")
                                                                     : QString()))
                            : ProfileWindow::tr("Manual axes"));
    strikeInfo->setToolTip(
        ProfileWindow::tr(
            "Observed phase tensors: %1 samples from %2 stations. The fourth circular moment averages "
            "principal axes modulo 90°. R₄ measures directional agreement, not statistical confidence. +90° "
            "chooses the other ambiguous axis. Skew and ellipticity thresholds always apply to automatic "
            "strike, even when the rose filter is off.")
            .arg(inferredStrike.samples)
            .arg(inferredStrike.stations));
    status->setText(
        vertices.size() < 2
            ? ProfileWindow::tr("Click at least two map vertices to define a profile. Right-click to finish.")
            : ProfileWindow::tr(
                  "%1 stations · %2 km profile · ±%3 km corridor · UTM distances · data x axis: %4")
                  .arg(offsets.size())
                  .arg(ProfileAnalysis::length(vertices), 0, 'f', 2)
                  .arg(distance->value(), 0, 'f', 3)
                  .arg(frame()));
    pdfButton->setEnabled(!offsets.empty() && std::isfinite(rotation));
  }
  void plots(bool fit) {
    updatingPlots = true;
    std::array<std::vector<ProfileGlyph>, 5> data;
    std::array<std::vector<ProfileAnalysis::Direction>, 3> angles;
    const auto *active = source();
    const double low = periodMin->currentData().toDouble(), high = periodMax->currentData().toDouble();
    updateStrike(low, high);
    const bool frequency = vertical->currentIndex() == 1;
    const std::array<Quantity, 2> components =
        quantities->currentIndex() == 0 ? std::array<Quantity, 2>{{Quantity::XY, Quantity::YX}}
        : quantities->currentIndex() == 1
            ? std::array<Quantity, 2>{{Quantity::Berdichevsky, Quantity::SSQ}}
            : std::array<Quantity, 2>{{Quantity::Determinant, Quantity::Determinant}};
    std::vector<double> along;
    for (const auto &s : offsets)
      along.push_back(s.along);
    const auto horizontal = ProfileAnalysis::cells(along, stationSpacing);
    for (auto *plot : sections)
      plot->setStations(offsets, stationSpacing);
    tensors->setStations(offsets, stationSpacing);
    for (unsigned i = 0; i < offsets.size(); ++i) {
      const auto &offset = offsets[i];
      if (!std::isfinite(rotation) || !active || !active->is_station_present(offset.name))
        continue;
      const auto &station = active->get_station_data(offset.name);
      if (!station.active())
        continue;
      std::vector<std::pair<double, unsigned>> periods;
      for (unsigned f = 0; f < station.frequencies().size(); ++f) {
        double p = 1. / station.frequencies()[f];
        if (std::isfinite(p) && p > 0.)
          periods.emplace_back(std::log10(p), f);
      }
      std::sort(periods.begin(), periods.end());
      std::vector<double> logs;
      for (const auto &p : periods)
        logs.push_back(p.first);
      const auto bounds = ProfileAnalysis::cells(logs, std::log10(gap->value()));
      for (unsigned j = 0; j < periods.size(); ++j) {
        const double p = std::pow(10., periods[j].first);
        if (p < low * (1. - 1e-10) || p > high * (1. + 1e-10))
          continue;
        const unsigned f = periods[j].second;
        ProfileGlyph glyph{QString::fromStdString(offset.name),
                           offset.along,
                           p,
                           horizontal[i][0],
                           horizontal[i][1],
                           bounds[j][0],
                           bounds[j][1],
                           NAN,
                           {}};
        for (unsigned c = 0; c < 2; ++c) {
          const auto sample = ImpedanceAnalysis::sample(station, components[c], f, rotation);
          glyph.value = sample.resistivity();
          if (std::isfinite(glyph.value) && glyph.value > 0.)
            data[c].push_back(glyph);
          glyph.value = sample.phase();
          if (std::isfinite(glyph.value))
            data[c + 2].push_back(glyph);
        }
        const bool tensor = MTMapData::phase_tensor(station, f, glyph.tensor);
        if (tensor) {
          glyph.tensor.azimuth = std::fmod(glyph.tensor.azimuth - rotation + 540., 180.);
          glyph.value = tensorColor->currentIndex() == 0   ? glyph.tensor.phiMax
                        : tensorColor->currentIndex() == 1 ? glyph.tensor.phiMin
                                                           : glyph.tensor.skew;
          data[4].push_back(glyph);
        }
        if (!inRoseBand(p))
          continue;
        if (tensor && glyph.tensor.axisRatio < 1. - 1e-8 &&
            (!filterTensor->isChecked() || tensorForStrike(glyph.tensor)))
          angles[0].push_back({offset.name, glyph.tensor.azimuth});
        for (unsigned component = 0; component < 2; ++component) {
          std::array<double, 2> v;
          if (MTMapData::induction_vector(station, f, component == 1, v) && std::hypot(v[0], v[1]) > 1e-12)
            angles[component + 1].push_back(
                {offset.name, std::atan2(v[1], v[0]) * 180. / std::acos(-1.) + 90. - rotation});
        }
      }
    }
    auto colors = PlotColorMaps::gradient(palette->currentText());
    for (unsigned i = 0; i < 4; ++i) {
      const bool rho = i < 2;
      sections[i]->setVisible(i % 2 == 0 || quantities->currentIndex() != 2);
      const auto title =
          QString(rho ? "ρa %1 (Ω m)" : "Phase %1 (°)")
              .arg(
                  QString::fromUtf8(ImpedanceAnalysis::name(components[i % 2])) +
                  (rotate->isChecked() && quantities->currentIndex() == 0 ? QStringLiteral("′") : QString()));
      sections[i]->setSamples(
          std::move(data[i]), frequency, false, cells->isChecked(), glyphSize->value(), title, colors,
          rho ? QCPRange(rhoMin->value(), rhoMax->value()) : QCPRange(phaseMin->value(), phaseMax->value()),
          rho);
    }
    tensors->setReference(frame());
    tensors->setSamples(std::move(data[4]), frequency, true, false, glyphSize->value(),
                        tensorColor->currentText(), colors, QCPRange(tensorMin->value(), tensorMax->value()),
                        false);
    if (fit) {
      const double left = along.empty() ? 0. : along.front() - .5 * stationSpacing;
      const double right = along.empty() ? 1. : along.back() + .5 * stationSpacing;
      const double bottom = low > 0. ? low / std::sqrt(gap->value()) : .1;
      const double top = high > 0. ? high * std::sqrt(gap->value()) : 100.;
      std::vector<ProfileSectionPlot *> linked(sections.begin(), sections.end());
      linked.push_back(tensors);
      for (auto *plot : linked) {
        plot->xAxis->setRange(left, right);
        plot->yAxis->setRange(frequency ? 1. / top : bottom, frequency ? 1. / bottom : top);
        plot->replot();
      }
    }
    updatingPlots = false;
    std::vector<ProfileRosePlot::Panel> panels;
    const unsigned binValues[]{5, 10, 15, 30};

    for (unsigned i = 0; i < 3; ++i)
      if (directions[i]->isChecked())
        panels.push_back(
            {directions[i]->text(), QColor::fromRgb(colors.color(.2 + .3 * i, QCPRange(0., 1.))),
             ProfileAnalysis::rose(angles[i], binValues[bins->currentIndex()], equalStations->isChecked())});
    roses->setPanels(std::move(panels),
                     ProfileWindow::tr("%1 · %2–%3 s · 0° = %4")
                         .arg(roseBand->currentText())
                         .arg(low, 0, 'g', 5)
                         .arg(high, 0, 'g', 5)
                         .arg(frame()),
                     rotate->isChecked());
  }
  void saveCsv() {
    const auto path = QFileDialog::getSaveFileName(owner, ProfileWindow::tr("Export profile stations"),
                                                   "profile-stations.csv", ProfileWindow::tr("CSV (*.csv)"));
    if (path.isEmpty())
      return;
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QMessageBox::warning(owner, ProfileWindow::tr("Export"), file.errorString());
      return;
    }
    QTextStream out(&file);
    out.setEncoding(QStringConverter::Utf8);
    out.setLocale(QLocale::c());
    out.setRealNumberPrecision(17);
    out << "# WGS84 / UTM " << coordinates.zone << (coordinates.north ? "N" : "S")
        << "; origin_easting_m=" << coordinates.origin_easting
        << "; origin_northing_m=" << coordinates.origin_northing << '\n';
    out << "# corridor_half_width_km=" << distance->value() << '\n';
    out << "# displayed_data_x_axis=" << frame()
        << "; strike_mode=" << (rotate->isChecked() ? strikeMode->currentText() : ProfileWindow::tr("Off"))
        << '\n';
    for (const auto &p : vertices)
      out << "# vertex_east_km=" << p[0] << ",north_km=" << p[1] << '\n';
    out << "station,along_km,offset_km,latitude,longitude\n";
    for (const auto &s : offsets) {
      const auto &p = survey->get_station_data(s.name).position();
      out << quoted(QString::fromStdString(s.name)) << ',' << s.along << ',' << s.across << ',' << p[0] << ','
          << p[1] << '\n';
    }
    out.flush();
    if (!file.commit())
      QMessageBox::warning(owner, ProfileWindow::tr("Export"), file.errorString());
  }
  void savePdf() {
    const auto path = QFileDialog::getSaveFileName(owner, ProfileWindow::tr("Export profile analysis"),
                                                   owner->exportName("pdf"), ProfileWindow::tr("PDF (*.pdf)"));
    if (path.isEmpty())
      return;
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
      QMessageBox::warning(owner, ProfileWindow::tr("Export"), file.errorString());
      return;
    }
    bool ok = true;
    const auto mapX = map->xAxis->range(), mapY = map->yAxis->range();
    {
      QPdfWriter pdf(&file);
      pdf.setCreator("EDITools");
      pdf.setTitle(owner->displayTitle());
      pdf.setResolution(120);
      pdf.setPageSize(QPageSize(QPageSize::A4));
      pdf.setPageOrientation(QPageLayout::Landscape);
      QCPPainter painter;
      if (!painter.begin(&pdf))
        ok = false;
      else {
        std::vector<std::pair<QString, QCustomPlot *>> pages{{ProfileWindow::tr("Profile corridor"), map}};
        for (unsigned i = 0; i < 4; ++i)
          if (i % 2 == 0 || quantities->currentIndex() != 2)
            pages.push_back({i < 2 ? ProfileWindow::tr("Apparent resistivity") : ProfileWindow::tr("Phase"),
                             sections[i]});
        pages.push_back({ProfileWindow::tr("Phase tensor ellipses — glyph up = %1").arg(frame()), tensors});
        pages.push_back({ProfileWindow::tr("Strike direction roses"), roses});
        for (unsigned i = 0; i < pages.size(); ++i) {
          if (i && !pdf.newPage()) {
            ok = false;
            break;
          }
          painter.setPen(Qt::black);
          painter.setFont(QFont("Sans Serif", 10));
          painter.drawText(QRect(0, 0, pdf.width(), 50), Qt::AlignCenter,
                           pages[i].first + "\n" + (i == 0 ? coordinateLabel->text() + " · " : QString()) +
                               status->text());
          painter.save();
          painter.translate(0, 60);
          pages[i].second->toPainter(&painter, pdf.width(), pdf.height() - 65);
          painter.restore();
        }
        painter.end();
      }
    }
    if (!ok || !file.commit())
      QMessageBox::warning(owner, ProfileWindow::tr("Export"),
                           ProfileWindow::tr("Could not write the profile PDF."));
    map->xAxis->setRange(mapX);
    map->yAxis->setRange(mapY);
    map->replot();
    for (auto *plot : sections)
      plot->replot();
    tensors->replot();
    roses->replot();
  }
};

ProfileWindow::ProfileWindow(QWidget *parent, std::function<void(std::shared_ptr<MTSurveyData>)> extracted, QCustomPlot *sharedMap)
    : WorkspaceView(parent), impl(new Impl(this, std::move(extracted), sharedMap)) {
  setObjectName("profileWindow");
  setWindowTitle(tr("Profile analysis"));
  resize(1500, 950);
}
ProfileWindow::~ProfileWindow() = default;
void ProfileWindow::setSharedMapActive(bool active, const SurveyCoordinates &coordinates) {
  auto &s = *impl;
  if(!active) { s.deactivateSharedMap(); return; }
  if(!s.sharedInteraction.activate()) return;
  s.sharedCoordinates = coordinates;
  s.sharedInteractions(); s.updateSharedMap();
}
void ProfileWindow::selectStation(const QString &name) { impl->highlight(name, false); }
void ProfileWindow::setMapLayers(int layers) {
  impl->background->setLayers(layers);
  for (auto *action : control<QToolButton>("mapLayers")->menu()->actions())
    if (action->isCheckable()) {
      QSignalBlocker block(action);
      action->setChecked(layers & action->data().toInt());
    }
  impl->map->replot();
}
void ProfileWindow::setData(const std::shared_ptr<MTSurveyData> &survey,
                            const std::map<std::string, MTSurveyData> &responses) {
  auto &s = *impl;
  const bool changed = s.survey != survey;
  const auto previous = s.coordinates.description();
  std::vector<std::array<double, 3>> oldLine;
  if(!changed && !s.vertices.empty()) {
    for(const auto &p: s.vertices) oldLine.push_back({{p[1] * 1000., p[0] * 1000., 0.}});
    oldLine = s.coordinates.inverse(oldLine);
  }
  if (changed) {
    s.vertices.clear();
    s.selectedStation.clear();
    s.automaticDistance = true;
  }
  s.survey = survey;
  s.responses.clear();
  s.stations.clear();
  {
    const auto key = s.dataset->currentData().toString();
    QSignalBlocker block(s.dataset);
    s.dataset->clear();
    s.dataset->addItem(tr("Observed data"), QString());
    for (const auto &entry : responses) {
      const auto name = QString::fromStdString(entry.first);
      s.responses[name] = &entry.second;
      s.dataset->addItem(FileLabels::fileName(name), name);
    }
    s.dataset->setCurrentIndex(std::max(0, s.dataset->findData(key)));
  }
  s.background->setTopography(survey ? survey->topography() : Topography::State{});
  QString error;
  try {
    if (survey) {
      std::vector<std::array<double, 3>> points;
      std::vector<std::string> names;
      for (const auto &name : survey->get_stations_names())
        if (survey->is_active(name)) {
          const auto &p = survey->get_station_data(name).position();
          if (std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) && p[0] >= -80. &&
              p[0] <= 84. && std::abs(p[1]) <= 180.) {
            points.push_back(p);
            names.push_back(name);
          }
        }
      if (!points.empty()) {
        s.coordinates = survey->coordinates();
        if (!s.coordinates.utm)
          s.coordinates = SurveyCoordinates::suggested(points);
        const auto projected = s.coordinates.transform(points);
        for (unsigned i = 0; i < projected.size(); ++i)
          s.stations.push_back({names[i], {{projected[i][1] / 1000., projected[i][0] / 1000.}}});
        s.background->setCoordinates(s.coordinates, .001);
      }
    }
  } catch (const std::exception &e) {
    error = QString::fromUtf8(e.what());
    s.stations.clear();
  }
  if(previous != s.coordinates.description() && !oldLine.empty()) {
    s.vertices.clear();
    for(const auto &p: s.coordinates.transform(oldLine)) s.vertices.push_back({{p[1] / 1000., p[0] / 1000.}});
  }
  s.coordinateLabel->setText(tr("WGS84 / UTM %1%2 · origin E %3, N %4 m")
                                 .arg(s.coordinates.zone)
                                 .arg(s.coordinates.north ? 'N' : 'S')
                                 .arg(s.coordinates.origin_easting, 0, 'f', 1)
                                 .arg(s.coordinates.origin_northing, 0, 'f', 1));
  s.stationSpacing = ProfileAnalysis::medianSpacing(s.stations);
  if (s.stationSpacing <= 0.)
    s.stationSpacing = 1.;
  if (s.automaticDistance)
    s.defaultDistance();
  s.periods();
  s.select();
  if (changed || previous != s.coordinates.description())
    s.fitMap();
  if (!error.isEmpty())
    s.status->setText(tr("Cannot project profile stations: %1").arg(error));
}

QJsonObject ProfileWindow::saveViewState() const {
  auto state = WorkspaceView::saveViewState();
  std::vector<std::array<double, 3>> projected;
  for(const auto &point: impl->vertices) projected.push_back({{point[1] * 1000., point[0] * 1000., 0.}});
  QJsonArray vertices;
  for(const auto &p: impl->coordinates.inverse(projected)) vertices.append(QJsonArray{p[0], p[1]});
  state["vertices"] = vertices; state["autoDistance"] = impl->automaticDistance;
  state["manualStrike"] = impl->manualStrike; state["manualStrikeSet"] = impl->manualStrikeSet;
  state["swappedStrike"] = impl->swappedStrike; state["station"] = impl->selectedStation;
  state["mapLayers"] = impl->background->layers();
  return state;
}
void ProfileWindow::restoreViewState(const QJsonObject &state) {
  WorkspaceView::restoreViewState(state);
  auto &s = *impl;
  s.automaticDistance = state["autoDistance"].toBool(true);
  s.manualStrike = state["manualStrike"].toDouble(); s.manualStrikeSet = state["manualStrikeSet"].toBool();
  s.swappedStrike = state["swappedStrike"].toBool();
  std::vector<std::array<double, 3>> geographic;
  for(const auto &value: state["vertices"].toArray()) {
    const auto point = value.toArray(); if(point.size() != 2) continue;
    geographic.push_back({{point[0].toDouble(), point[1].toDouble(), 0.}});
  }
  s.vertices.clear();
  try { for(const auto &p: s.coordinates.transform(geographic)) s.vertices.push_back({{p[1] / 1000., p[0] / 1000.}}); }
  catch(const std::exception &e) { s.status->setText(tr("Could not restore profile: %1").arg(e.what())); }
  s.select(); s.highlight(state["station"].toString());
  if(state.contains("mapLayers")) setMapLayers(state["mapLayers"].toInt());
  restorePlotRanges(this, state["plots"].toObject());
}
