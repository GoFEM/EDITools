#pragma once
#include <set>
#include <string>
class QWidget;
class MTSurveyData;

// Applies only after acceptance; returns whether any mask changed.
bool editGroupMasks(QWidget *parent, MTSurveyData &survey,
                    const std::set<std::string> &stations, bool linkTensorMasks);
