#ifndef EDITOOLS_EXPORT_PREVIEW_H
#define EDITOOLS_EXPORT_PREVIEW_H
#include "include/SurveyExport.h"
#include <QWidget>
class QLabel;
class QTabWidget;
class QTableWidget;
class ExportPreview : public QWidget {
public:
  explicit ExportPreview(QWidget *parent = nullptr);
  bool update(const MTSurveyData &survey, const NativeMT::Options &options, SurveyExport::Format format);
  const NativeMT::Export &result() const { return prepared; }

private:
  NativeMT::Export prepared;
  QLabel *summary, *details;
  QTableWidget *stations, *periods, *components;
};
// Returns the chosen observation path, or empty on cancellation. Stages all companions before committing.
QString saveSurveyExport(QWidget *parent, const NativeMT::Export &result, SurveyExport::Format format,
                         const QString &directory);
#endif
