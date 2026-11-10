#pragma once
#include "include/AutoMask.h"
class QWidget;

// Review a single analysis pass, then apply checked candidates. Caller records history.
bool reviewAutoMasks(QWidget *parent, MTSurveyData &survey, const std::set<std::string> &selected,
                     bool linkTensorMasks, const ImpedanceAnalysis::RhoPlusOptions &rhoPlus = {},
                     AutoMask::Options *preferences = nullptr, bool *pairPreference = nullptr);
