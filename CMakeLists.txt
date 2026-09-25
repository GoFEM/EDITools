#
# The EDI Tools application.
#
# Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(EDITools LANGUAGES CXX)

find_package(Qt5 5.12 REQUIRED COMPONENTS Widgets PrintSupport Network)
find_package(PROJ REQUIRED CONFIG)
find_package(Boost REQUIRED CONFIG COMPONENTS serialization filesystem system)

set(QCUSTOMPLOT_PATH "${QCUSTOMPLOT_PATH}" CACHE PATH "Directory containing qcustomplot.cpp and qcustomplot.h")
set(EIGEN_PATH "${EIGEN_PATH}" CACHE PATH "Directory containing Eigen headers")
if(NOT EXISTS "${QCUSTOMPLOT_PATH}/qcustomplot.cpp" OR NOT EXISTS "${QCUSTOMPLOT_PATH}/qcustomplot.h")
  message(FATAL_ERROR "Set QCUSTOMPLOT_PATH to the QCustomPlot source directory.")
endif()
if(NOT EXISTS "${EIGEN_PATH}/Eigen/Core")
  message(FATAL_ERROR "Set EIGEN_PATH to the Eigen include directory.")
endif()

# Compile the numerical/data code once, without a Qt dependency.
add_library(editools_data STATIC
  src/Datum.cc
  src/MTSpectra.cc
  src/MTStationData.cc
  src/NativeMT.cc
  src/MTResponseData.cc
  src/FitStatistics.cc
  src/MTMapData.cc
  src/PeriodResampling.cc
  src/EDIPeriodMerge.cc
  src/SurveyCoordinates.cc
  src/MTSurveyData.cc
  src/EDIFileReader.cc
)
target_compile_features(editools_data PUBLIC cxx_std_14)
target_include_directories(editools_data PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(editools_data SYSTEM PRIVATE ${EIGEN_PATH})
target_link_libraries(editools_data PUBLIC Boost::serialization PROJ::proj PRIVATE Boost::filesystem Boost::system)

# QCustomPlot is shared by the application and every graphical test.
add_library(editools_qcustomplot STATIC "${QCUSTOMPLOT_PATH}/qcustomplot.cpp" "${QCUSTOMPLOT_PATH}/qcustomplot.h")
set_target_properties(editools_qcustomplot PROPERTIES AUTOMOC ON)
target_include_directories(editools_qcustomplot SYSTEM PUBLIC ${QCUSTOMPLOT_PATH})
target_link_libraries(editools_qcustomplot PUBLIC Qt5::Widgets Qt5::PrintSupport)

add_library(editools_ui STATIC
  src/MTDataPlot.cc
  src/ApparentResistivityPlot.cc
  src/PhasePlot.cc
  src/TipperPlot.cc
  src/PhaseTensorPlot.cc
  src/PlotColorMaps.cc
  src/OsmTiles.cc
  src/OsmBasemap.cc
  src/MapBackground.cc
  src/MapPlot.cc
  ExportGOFEMDialog.cpp
  ExportNativeMTDialog.cpp
  SurveyCoordinatesDialog.cpp
  FitStatisticsWindow.cpp
  PeriodMapWindow.cpp
  PeriodLayoutWindow.cpp
  EDIImportDialog.cpp
  SurveyReport.cpp
  mainwindow.cpp
  mainwindow.h
  ExportGOFEMDialog.h
  include/MTDataPlot.h
  include/OsmTiles.h
  mainwindow.ui
  ExportGOFEMDialog.ui
)
set_target_properties(editools_ui PROPERTIES AUTOMOC ON AUTOUIC ON)
target_link_libraries(editools_ui PUBLIC editools_data editools_qcustomplot Qt5::Network)

# Keep resources on the executables so static-library linking cannot discard them.
add_executable(main main.cpp Icons/icons.qrc)
set_target_properties(main PROPERTIES AUTORCC ON)
target_link_libraries(main PRIVATE editools_ui)

include(CTest)
if(BUILD_TESTING)
  add_executable(mt_spectra_test tests/mt_spectra_test.cc)
  target_link_libraries(mt_spectra_test PRIVATE editools_data)
  target_compile_definitions(mt_spectra_test PRIVATE
    SPECTRA_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
    SPECTRA_TEST_OUTPUT="${CMAKE_CURRENT_BINARY_DIR}")
  add_test(NAME mt_spectra COMMAND mt_spectra_test)

  add_executable(native_mt_test tests/native_mt_test.cc)
  target_link_libraries(native_mt_test PRIVATE editools_data)
  target_compile_definitions(native_mt_test PRIVATE NATIVE_MT_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
  add_test(NAME native_mt_format COMMAND native_mt_test)

  add_executable(native_mt_dialog_test tests/native_mt_dialog_test.cc)
  target_link_libraries(native_mt_dialog_test PRIVATE editools_ui)
  add_test(NAME native_mt_dialog COMMAND native_mt_dialog_test)

  add_executable(osm_test tests/osm_test.cc)
  target_link_libraries(osm_test PRIVATE editools_ui)
  add_test(NAME osm_basemap COMMAND osm_test)

  add_executable(response_plot_test tests/response_plot_test.cc tests/period_resampling_test.cc Icons/icons.qrc)
  set_target_properties(response_plot_test PROPERTIES AUTORCC ON)
  target_link_libraries(response_plot_test PRIVATE editools_ui)
  add_test(NAME response_plots COMMAND response_plot_test)

  set_tests_properties(native_mt_dialog osm_basemap response_plots PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
  set_tests_properties(native_mt_dialog PROPERTIES TIMEOUT 20)
  set_tests_properties(osm_basemap PROPERTIES TIMEOUT 30)
  set_tests_properties(response_plots PROPERTIES TIMEOUT 60)
endif()
