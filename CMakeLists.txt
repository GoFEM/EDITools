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

set(TARGET main)
project(${TARGET} CXX Fortran)

find_package(Qt5Widgets)
find_package(PROJ REQUIRED CONFIG)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Set your files and resources here
set( Srcs ${QCUSTOMPLOT_PATH}/qcustomplot.cpp
          src/mtcomp.f
          src/mtcomp_interface.f90
          src/MTStationData.cc
          src/NativeMT.cc
          src/MTResponseData.cc
          src/FitStatistics.cc
          src/MTMapData.cc
          src/PeriodResampling.cc
          src/SurveyCoordinates.cc
          src/MTSurveyData.cc
          src/EDIFileReader.cc
          src/MTDataPlot.cc
          src/ApparentResistivityPlot.cc
          src/PhasePlot.cc
          src/TipperPlot.cc
          src/PhaseTensorPlot.cc
          src/MapPlot.cc
          ExportGOFEMDialog.cpp
          ExportNativeMTDialog.cpp
          SurveyCoordinatesDialog.cpp
          FitStatisticsWindow.cpp
          PeriodMapWindow.cpp
          PeriodLayoutWindow.cpp
          SurveyReport.cpp
          mainwindow.cpp
          ${TARGET}.cpp)

set( Hdrs mainwindow.h ExportGOFEMDialog.h include/MTDataPlot.h )

set( MOC_Hdrs mainwindow.h ExportGOFEMDialog.h )

set( UIs mainwindow.ui ExportGOFEMDialog.ui )

set( QRCs Icons/icons.qrc )


# The rest should just work (sure...)
# We have ui files, this will bring in the macro: qt5_wrap_ui
find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QRC_Srcs ${QRCs} )

  source_group("Resources" FILES
    ${UIs}
    ${QRCs}
  )

  source_group("Generated" FILES
    ${UI_Srcs}
    ${MOC_Srcs}
    ${QRC_Srcs}
  )

find_package(Boost REQUIRED COMPONENTS serialization filesystem system)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)

include_directories(${QCUSTOMPLOT_PATH} ${Boost_INCLUDE_DIRS} ${EIGEN_PATH} "/include")

#SET(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(${TARGET} ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs}
    include/datum.h
    README.md)
qt5_use_modules(${TARGET} Core Gui Widgets PrintSupport)
target_link_libraries(${TARGET} ${Boost_LIBRARIES})
target_link_libraries(${TARGET} PROJ::proj)

include(CTest)
if(BUILD_TESTING)
  add_executable(native_mt_test tests/native_mt_test.cc src/NativeMT.cc src/MTResponseData.cc src/MTStationData.cc)
  set_target_properties(native_mt_test PROPERTIES AUTOMOC OFF AUTOUIC OFF)
  target_include_directories(native_mt_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(native_mt_test ${Boost_LIBRARIES})
  target_compile_definitions(native_mt_test PRIVATE NATIVE_MT_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
  add_test(NAME native_mt_format COMMAND native_mt_test)
  add_executable(native_mt_dialog_test tests/native_mt_dialog_test.cc ExportNativeMTDialog.cpp
    SurveyCoordinatesDialog.cpp src/SurveyCoordinates.cc
    src/NativeMT.cc src/MTResponseData.cc src/MTStationData.cc src/MTSurveyData.cc src/EDIFileReader.cc
    src/mtcomp.f src/mtcomp_interface.f90)
  set_target_properties(native_mt_dialog_test PROPERTIES AUTOMOC OFF AUTOUIC OFF
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/native_mt_dialog_modules")
  target_include_directories(native_mt_dialog_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(native_mt_dialog_test Qt5::Widgets ${Boost_LIBRARIES} PROJ::proj)
  add_test(NAME native_mt_dialog COMMAND native_mt_dialog_test)
  set_tests_properties(native_mt_dialog PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 20)
  set(ResponseTestSrcs ${Srcs})
  list(REMOVE_ITEM ResponseTestSrcs ${TARGET}.cpp)
  add_executable(response_plot_test tests/response_plot_test.cc tests/period_resampling_test.cc ${ResponseTestSrcs}
    ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs})
  set_target_properties(response_plot_test PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/response_plot_modules")
  qt5_use_modules(response_plot_test Core Gui Widgets PrintSupport)
  target_link_libraries(response_plot_test ${Boost_LIBRARIES} PROJ::proj)
  add_test(NAME response_plots COMMAND response_plot_test)
  set_tests_properties(response_plots PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()
