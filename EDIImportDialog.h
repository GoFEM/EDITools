#ifndef EDI_IMPORT_DIALOG_H
#define EDI_IMPORT_DIALOG_H
#include "include/MTSurveyData.h"
class QWidget;
// Review a staged import. Cancellation leaves the candidate and live survey untouched.
bool reviewEDIImport(QWidget *parent, MTSurveyData &candidate,
                     const std::set<std::string> &incoming, const std::vector<std::string> &duplicates);
#endif
