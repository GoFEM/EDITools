#ifndef EDITOOLS_SURVEY_EXPORT_H
#define EDITOOLS_SURVEY_EXPORT_H
#include "NativeMT.h"
class MTSurveyData;
namespace SurveyExport {
enum class Format { Native, GoFEM };
// Prepare once for preview and writing. Does not change survey coordinates or data.
NativeMT::Export prepare(const MTSurveyData &survey, NativeMT::Options options, Format format);
} // namespace SurveyExport
#endif
