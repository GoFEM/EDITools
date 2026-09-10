#ifndef EXPORT_NATIVE_MT_DIALOG_H
#define EXPORT_NATIVE_MT_DIALOG_H

#include <QDialog>
#include "include/NativeMT.h"

class MTSurveyData;
class QListWidget;
class QCheckBox;

class ExportNativeMTDialog : public QDialog
{
public:
  ExportNativeMTDialog(MTSurveyData &survey, const QString &directory, QWidget *parent = nullptr);
  QString exportedFile() const { return exported_file; }

private:
  void exportFiles();
  MTSurveyData &survey;
  QString directory, exported_file;
  QListWidget *periods, *types;
  QCheckBox *conjugate, *reverse_vertical, *distortion;
};
#endif
