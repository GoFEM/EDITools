
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET main)
project(${TARGET} CXX Fortran)

find_package(Qt5Widgets)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Set your files and resources here
set( Srcs ${QCUSTOMPLOT_PATH}/qcustomplot.cpp
          mtcomp.f
          mtcomp_interface.f90
          MTStationData.cc
          MTSurveyData.cc
          EDIFileReader.cc
          MTDataPlot.cc
          ApparentResistivityPlot.cc
          PhasePlot.cc
          TipperPlot.cc
          PhaseTensorPlot.cc
          MapPlot.cc
          ExportGOFEMDialog.cpp
          mainwindow.cpp
          ${TARGET}.cpp)

set( Hdrs mainwindow.h ExportGOFEMDialog.h )

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

include_directories(${QCUSTOMPLOT_PATH} ${Boost_INCLUDE_DIRS} ${EIGEN_PATH})

#SET(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(${TARGET} ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs}
    datum.h
    README.md)
qt5_use_modules(${TARGET} Core Gui Widgets PrintSupport)
target_link_libraries(${TARGET} ${Boost_LIBRARIES})
