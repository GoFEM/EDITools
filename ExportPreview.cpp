#include "ExportPreview.h"
#include "include/MTSurveyData.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QLabel>
#include <QMessageBox>
#include <QSaveFile>
#include <QTableWidget>
#include <QTabWidget>
#include <QVBoxLayout>
#include <map>
#include <cmath>
#include <set>

ExportPreview::ExportPreview(QWidget *parent) : QWidget(parent) {
  setObjectName("exportPreview");
  auto *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);
  summary = new QLabel(this);
  summary->setObjectName("exportSummary");
  summary->setWordWrap(true);
  layout->addWidget(summary);
  auto *tabs = new QTabWidget(this);
  layout->addWidget(tabs);
  auto table = [tabs](const QString &title, const QStringList &headers) {
    auto *t = new QTableWidget(tabs);
    t->setColumnCount(headers.size());
    t->setHorizontalHeaderLabels(headers);
    t->setEditTriggers(QAbstractItemView::NoEditTriggers);
    t->verticalHeader()->hide();
    t->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    tabs->addTab(t, title);
    return t;
  };
  stations = table(tr("Stations"), {tr("Station"), tr("Scalars"), tr("Periods")});
  stations->setObjectName("exportStations");
  periods = table(tr("Periods"), {tr("Period (s)"), tr("Scalars"), tr("Stations")});
  periods->setObjectName("exportPeriods");
  components = table(tr("Components"), {tr("Scalar"), tr("Count"), tr("Min σ"), tr("Max σ")});
  components->setObjectName("exportComponents");
  details = new QLabel(this);
  details->setObjectName("exportDetails");
  details->setWordWrap(true);
  details->setTextInteractionFlags(Qt::TextSelectableByMouse);
  layout->addWidget(details);
  tabs->setMinimumHeight(120);
  tabs->setMaximumHeight(190);
  setToolTip(tr("Counts describe exactly the prepared output, after station, period and component selection "
                "and masking. Complex values count as two scalars. Uncertainties use the writer's effective "
                "errors, in the same units as each scalar. Invalid values or errors prevent export."));
}
bool ExportPreview::update(const MTSurveyData &survey, const NativeMT::Options &options,
                           SurveyExport::Format format) {
  for (auto *t : {stations, periods, components})
    t->setRowCount(0);
  prepared = {};
  details->clear();
  try {
    prepared = SurveyExport::prepare(survey, options, format);
    struct Count {
      size_t scalars = 0;
      std::set<double> frequencies;
      std::set<std::string> names;
      double low = INFINITY, high = 0.;
    };
    std::map<std::string, Count> byStation;
    std::map<double, Count> byPeriod;
    std::map<RealDataType, Count> byComponent;
    for (const auto &row : prepared.observations) {
      for (auto *count : {&byStation[row.receiver], &byPeriod[1. / row.frequency], &byComponent[row.type]}) {
        ++count->scalars;
        count->frequencies.insert(row.frequency);
        count->names.insert(row.receiver);
        count->low = std::min(count->low, row.error);
        count->high = std::max(count->high, row.error);
      }
    }
    auto append = [](QTableWidget *table, const QStringList &cells) {
      int r = table->rowCount();
      table->insertRow(r);
      for (int c = 0; c < cells.size(); ++c)
        table->setItem(r, c, new QTableWidgetItem(cells[c]));
    };
    for (const auto &entry : byStation)
      append(stations, {QString::fromStdString(entry.first), QString::number(entry.second.scalars),
                        QString::number(entry.second.frequencies.size())});
    for (const auto &entry : byPeriod)
      append(periods, {QString::number(entry.first, 'g', 12), QString::number(entry.second.scalars),
                       QString::number(entry.second.names.size())});
    for (const auto &entry : byComponent)
      append(components, {QString::fromStdString(Datum::convert_type_to_string(entry.first)),
                          QString::number(entry.second.scalars), QString::number(entry.second.low, 'g', 6),
                          QString::number(entry.second.high, 'g', 6)});
    unsigned disabled = 0;
    for (const auto &name : survey.get_stations_names())
      disabled += !survey.is_active(name);
    summary->setStyleSheet({});
    summary->setText(tr("%1 scalars · %2 stations · %3 periods\nOmitted: %4 disabled stations · %5 masked "
                        "scalars · %6 missing · %7 undefined phases")
                         .arg(prepared.observations.size())
                         .arg(prepared.receivers.size())
                         .arg(prepared.frequencies.size())
                         .arg(disabled)
                         .arg(prepared.masked)
                         .arg(prepared.missing)
                         .arg(prepared.undefined_phase));
    QString description = QString::fromStdString(prepared.coordinate_description);
    details->setText(tr("Z: Ω · ρ: Ω m · phase: degrees · PT/tipper: dimensionless\nNorth / East (m) · z %1 "
                        "· conjugate %2 · reverse tipper %3")
                         .arg(format == SurveyExport::Format::Native ? tr("down") : tr("up"),
                              options.negative_time_convention ? tr("yes") : tr("no"),
                              options.vertical_magnetic_up ? tr("yes") : tr("no")));
    details->setToolTip(description);
    // Keep the origin visible; full coordinate provenance is available on hover.
    const auto lines = description.split('\n');
    for (const auto &line : lines)
      if (line.startsWith("WGS84") || line.startsWith("UTM origin"))
        details->setText(details->text() + "\n" + line);
    for (const auto &warning : prepared.warnings)
      details->setText(details->text() + "\n" + QString::fromStdString(warning));
    return true;
  } catch (const std::exception &e) {
    summary->setStyleSheet("color: #a45216");
    summary->setText(QString::fromUtf8(e.what()));
    return false;
  }
}
QString saveSurveyExport(QWidget *parent, const NativeMT::Export &result, SurveyExport::Format format,
                         const QString &directory) {
  QString path = QFileDialog::getSaveFileName(
      parent, QObject::tr("Export observations"), directory,
      format == SurveyExport::Format::Native ? QObject::tr("MT observations (*.data);;All files (*)")
                                             : QObject::tr("GoFEM observations (*.gofem);;All files (*)"));
  if (path.isEmpty())
    return {};
  if (QFileInfo(path).suffix().isEmpty())
    path += format == SurveyExport::Format::Native ? ".data" : ".gofem";
  const QFileInfo info(path);
  const QString base =
      format == SurveyExport::Format::Native ? info.absolutePath() + "/" + info.completeBaseName() : path;
  const std::array<QString, 3> paths{{path, base + ".recvs", base + ".freqs"}};
  if (paths[0] == paths[1] || paths[0] == paths[2])
    throw std::runtime_error("Observation filename conflicts with its companion files.");
  QStringList existing;
  for (unsigned i = 1; i < paths.size(); ++i)
    if (QFileInfo::exists(paths[i]))
      existing << paths[i];
  if (!existing.empty() && QMessageBox::question(parent, QObject::tr("Replace companion files?"),
                                                 existing.join("\n")) != QMessageBox::Yes)
    return {};
  const std::array<std::string, 3> texts{{result.data_text, result.receiver_text, result.frequency_text}};
  std::array<std::unique_ptr<QSaveFile>, 3> files;
  for (unsigned i = 0; i < paths.size(); ++i) {
    files[i].reset(new QSaveFile(paths[i]));
    if (!files[i]->open(QIODevice::WriteOnly) ||
        files[i]->write(texts[i].data(), qint64(texts[i].size())) != qint64(texts[i].size()))
      throw std::runtime_error("Cannot stage export: " + files[i]->errorString().toStdString());
  }
  for (auto &file : files)
    if (!file->commit())
      throw std::runtime_error(
          "Cannot finalize export; some files may have been replaced. Retry the complete export. " +
          file->errorString().toStdString());
  return path;
}
