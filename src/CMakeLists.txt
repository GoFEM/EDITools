# Compile the numerical/data code once, without a Qt dependency.
add_library(editools_data STATIC
  core/AutoMask.cpp
  core/AutoMask.h
  core/Datum.cpp
  core/Datum.h
  core/EDIFileReader.cpp
  core/EDIFileReader.h
  core/EDIPeriodMerge.cpp
  core/EDIPeriodMerge.h
  core/FitStatistics.cpp
  core/FitStatistics.h
  core/ImpedanceAnalysis.cpp
  core/ImpedanceAnalysis.h
  core/MTMapData.cpp
  core/MTMapData.h
  core/MTResponseData.cpp
  core/MTResponseData.h
  core/MTSpectra.cpp
  core/MTSpectra.h
  core/MTStationData.cpp
  core/MTStationData.h
  core/MTSurveyData.cpp
  core/MTSurveyData.h
  core/NativeMT.cpp
  core/NativeMT.h
  core/PeriodResampling.cpp
  core/PeriodResampling.h
  core/PeriodResamplingInfo.h
  core/PhaseStatistics.h
  core/PolarityAnalysis.cpp
  core/PolarityAnalysis.h
  core/ProfileAnalysis.cpp
  core/ProfileAnalysis.h
  core/SurveyCoordinates.cpp
  core/SurveyCoordinates.h
  core/SurveyExport.cpp
  core/SurveyExport.h
  core/SurveyGeometry.h
  core/Topography.cpp
  core/Topography.h
  core/UTMAnalysis.cpp
  core/UTMAnalysis.h
)
target_compile_features(editools_data PUBLIC cxx_std_14)
if(GDAL_VERSION VERSION_GREATER_EQUAL 3.9)
  target_compile_features(editools_data PUBLIC cxx_std_17)
endif()
target_include_directories(editools_data PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(editools_data SYSTEM PRIVATE ${EIGEN_PATH})
target_link_libraries(editools_data PUBLIC Boost::serialization PROJ::proj PRIVATE Boost::filesystem Boost::system GDAL::GDAL)

# QCustomPlot is shared by the application and every graphical test.
add_library(editools_qcustomplot STATIC "${QCUSTOMPLOT_PATH}/qcustomplot.cpp" "${QCUSTOMPLOT_PATH}/qcustomplot.h")
set_target_properties(editools_qcustomplot PROPERTIES AUTOMOC ON)
target_include_directories(editools_qcustomplot SYSTEM PUBLIC ${QCUSTOMPLOT_PATH})
target_link_libraries(editools_qcustomplot PUBLIC Qt6::Widgets Qt6::PrintSupport)

add_library(editools_ui STATIC
  "${PROJECT_SOURCE_DIR}/resources/icons/icons.qrc"
  app/MainWindow.cpp
  app/MainWindow.h
  app/MainWindow.ui
  app/MainWindowView.cpp
  app/MainWindowWorkspace.cpp
  ui/FileLabels.h
  ui/SurveyReport.cpp
  ui/SurveyReport.h
  ui/ToolTips.h
  ui/analysis/FitStatisticsWindow.cpp
  ui/analysis/FitStatisticsWindow.h
  ui/analysis/ImpedanceAnalysisWindows.cpp
  ui/analysis/ImpedanceAnalysisWindows.h
  ui/analysis/PeriodLayoutWindow.cpp
  ui/analysis/PeriodLayoutWindow.h
  ui/analysis/PeriodMapWindow.cpp
  ui/analysis/PeriodMapWindow.h
  ui/analysis/PolarityWindow.cpp
  ui/analysis/PolarityWindow.h
  ui/analysis/ProfileWindow.cpp
  ui/analysis/ProfileWindow.h
  ui/analysis/StationComparisonWindow.cpp
  ui/analysis/StationComparisonWindow.h
  ui/dialogs/EDIImportDialog.cpp
  ui/dialogs/EDIImportDialog.h
  ui/dialogs/ExportGOFEMDialog.cpp
  ui/dialogs/ExportGOFEMDialog.h
  ui/dialogs/ExportGOFEMDialog.ui
  ui/dialogs/ExportNativeMTDialog.cpp
  ui/dialogs/ExportNativeMTDialog.h
  ui/dialogs/ExportPreview.cpp
  ui/dialogs/ExportPreview.h
  ui/dialogs/SurveyCoordinatesDialog.cpp
  ui/dialogs/SurveyCoordinatesDialog.h
  ui/dialogs/TopographyDialog.cpp
  ui/dialogs/TopographyDialog.h
  ui/maps/MapBackground.cpp
  ui/maps/MapBackground.h
  ui/maps/MapDisplaySettings.h
  ui/maps/MapPlot.cpp
  ui/maps/MapPlot.h
  ui/maps/OsmBasemap.cpp
  ui/maps/OsmBasemap.h
  ui/maps/OsmTiles.cpp
  ui/maps/OsmTiles.h
  ui/maps/SharedMapInteraction.h
  ui/maps/TopographyLayer.cpp
  ui/maps/TopographyLayer.h
  ui/masking/AutoMaskDialog.cpp
  ui/masking/AutoMaskDialog.h
  ui/masking/AutoMaskSettings.cpp
  ui/masking/AutoMaskSettings.h
  ui/masking/AutoMaskWorkflow.cpp
  ui/masking/AutoMaskWorkflow.h
  ui/masking/GroupMaskDialog.cpp
  ui/masking/GroupMaskDialog.h
  ui/masking/GroupMaskPlots.cpp
  ui/masking/GroupMaskPlots.h
  ui/plots/ApparentResistivityPlot.cpp
  ui/plots/ApparentResistivityPlot.h
  ui/plots/MTDataPlot.cpp
  ui/plots/MTDataPlot.h
  ui/plots/PhasePlot.cpp
  ui/plots/PhasePlot.h
  ui/plots/PhaseTensorPlot.cpp
  ui/plots/PhaseTensorPlot.h
  ui/plots/PlotColorMaps.cpp
  ui/plots/PlotColorMaps.h
  ui/plots/ProfilePlots.cpp
  ui/plots/ProfilePlots.h
  ui/plots/RhoPlusOverlay.cpp
  ui/plots/RhoPlusOverlay.h
  ui/plots/TFViewStyle.h
  ui/plots/TipperPlot.cpp
  ui/plots/TipperPlot.h
  ui/workspace/CommandIcons.cpp
  ui/workspace/CommandIcons.h
  ui/workspace/CommandRibbon.cpp
  ui/workspace/CommandRibbon.h
  ui/workspace/StationListDelegate.cpp
  ui/workspace/StationListDelegate.h
  ui/workspace/SurveyHistory.cpp
  ui/workspace/SurveyHistory.h
  ui/workspace/TFViewControls.cpp
  ui/workspace/TFViewControls.h
  ui/workspace/WorkspaceState.cpp
  ui/workspace/WorkspaceTabs.cpp
  ui/workspace/WorkspaceTabs.h
  ui/workspace/WorkspaceView.cpp
  ui/workspace/WorkspaceView.h
)
set_target_properties(editools_ui PROPERTIES AUTOMOC ON AUTOUIC ON AUTORCC ON)
target_link_libraries(editools_ui PUBLIC editools_data editools_qcustomplot Qt6::Network)

add_executable(main app/main.cpp)
if(WIN32)
  set_target_properties(main PROPERTIES WIN32_EXECUTABLE ON)
endif()
target_link_libraries(main PRIVATE editools_ui)
