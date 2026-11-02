#ifndef EDITOOLS_EXPORT_PREVIEW_H
#define EDITOOLS_EXPORT_PREVIEW_H
#include "include/SurveyExport.h"
#include <QWidget>
class QLabel;
class QCheckBox;
class QTabWidget;
class QTableWidget;
class ExportPreview : public QWidget {
  Q_OBJECT
public:
  explicit ExportPreview(QWidget *parent = nullptr);
  bool update(const MTSurveyData &survey, const NativeMT::Options &options, SurveyExport::Format format);
  const NativeMT::Export &result() const { return prepared; }
signals:
  void optionsChanged();

private:
  NativeMT::Export prepared;
  QLabel *summary, *details;
  QCheckBox *compatibleNames;
  QTableWidget *stations, *periods, *components;
};
// Returns the chosen observation path, or empty on cancellation. Stages all companions before committing.
QString saveSurveyExport(QWidget *parent, const NativeMT::Export &result, SurveyExport::Format format,
                         const QString &directory);
#endif
