#include "SurveyCoordinatesDialog.h"
#include "include/MTSurveyData.h"
#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QVBoxLayout>

SurveyCoordinatesDialog::SurveyCoordinatesDialog(MTSurveyData &data, QWidget *parent)
  : QDialog(parent), survey(data)
{
  setWindowTitle(tr("Convert coordinates to UTM"));
  resize(760, 500);
  auto *layout = new QVBoxLayout(this);
  auto *form = new QFormLayout;
  layout->addLayout(form);
  const auto settings = survey.coordinates().utm ? survey.coordinates() : SurveyCoordinates::suggested(survey.geographic_locations());
  zone = new QSpinBox(this); zone->setRange(1, 60); zone->setValue(settings.zone);
  zone->setObjectName("utmZone");
  hemisphere = new QComboBox(this); hemisphere->addItems({tr("North"), tr("South")});
  hemisphere->setCurrentIndex(settings.north ? 0 : 1);
  hemisphere->setObjectName("utmHemisphere");
  form->addRow(tr("WGS84 / UTM zone"), zone);
  form->addRow(tr("Hemisphere"), hemisphere);
  center = new QCheckBox(tr("Center coordinates on the survey"), this);
  center->setObjectName("utmCenter");
  center->setChecked(settings.centered);
  form->addRow(center);
  easting = new QLineEdit(this); northing = new QLineEdit(this);
  easting->setObjectName("utmOriginEasting"); northing->setObjectName("utmOriginNorthing");
  easting->setReadOnly(true); northing->setReadOnly(true);
  form->addRow(tr("Origin UTM easting (m)"), easting);
  form->addRow(tr("Origin UTM northing (m)"), northing);
  copy = new QPushButton(tr("Copy origin"), this);
  copy->setObjectName("copyUtmOrigin");
  form->addRow(copy);
  auto *help = new QLabel(tr("X = northing − origin northing; Y = easting − origin easting.\nThe center is the midpoint of all stations' UTM bounds. Elevations and response tensors are unchanged."), this);
  help->setWordWrap(true);
  layout->addWidget(help);
  table = new QTableWidget(this);
  table->setColumnCount(5);
  table->setHorizontalHeaderLabels({tr("Station"), tr("UTM easting (m)"), tr("UTM northing (m)"), tr("X / North (m)"), tr("Y / East (m)")});
  table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  table->setEditTriggers(QAbstractItemView::NoEditTriggers);
  layout->addWidget(table);
  status = new QLabel(this); status->setWordWrap(true); layout->addWidget(status);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
  apply = buttons->button(QDialogButtonBox::Ok); apply->setText(tr("Use these coordinates"));
  layout->addWidget(buttons);
  connect(zone, QOverload<int>::of(&QSpinBox::valueChanged), this, [this] {updatePreview();});
  connect(hemisphere, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] {updatePreview();});
  connect(center, &QCheckBox::toggled, this, [this] {updatePreview();});
  connect(copy, &QPushButton::clicked, this, [this] {QApplication::clipboard()->setText(QString::fromStdString(preview.description()));});
  connect(buttons, &QDialogButtonBox::accepted, this, [this] {survey.set_coordinates(preview); accept();});
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  updatePreview();
}

void SurveyCoordinatesDialog::updatePreview()
{
  try {
    const auto geographic = survey.geographic_locations();
    preview = SurveyCoordinates::calculate(geographic, zone->value(), hemisphere->currentIndex() == 0, center->isChecked());
    easting->setText(QString::number(preview.origin_easting, 'g', 17));
    northing->setText(QString::number(preview.origin_northing, 'g', 17));
    const auto positions = preview.transform(geographic);
    const auto names = survey.get_stations_names();
    table->setRowCount(int(names.size()));
    for(unsigned i = 0; i < names.size(); ++i) {
      table->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(names[i])));
      const double values[] = {positions[i][1]+preview.origin_easting, positions[i][0]+preview.origin_northing, positions[i][0], positions[i][1]};
      for(unsigned c = 0; c < 4; ++c) table->setItem(i, c+1, new QTableWidgetItem(QString::number(values[c], 'f', 3)));
    }
    status->setText(tr("Applies to maps and data exports. Original latitude/longitude are retained."));
    apply->setEnabled(true); copy->setEnabled(true);
  } catch(const std::exception &e) {
    status->setText(QString::fromUtf8(e.what()));
    apply->setEnabled(false); copy->setEnabled(false);
    table->setRowCount(0); easting->clear(); northing->clear();
  }
}
