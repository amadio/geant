#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB lsources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB lheaders ${PROJECT_SOURCE_DIR}/include/TM*.h ${PROJECT_SOURCE_DIR}/include/TP*.h)
file(GLOB chdrs   ${PROJECT_SOURCE_DIR}/include/T*.h)
set(cdict TXsecCint)
ROOT_GENERATE_DICTIONARY(${cdict} ${chdrs} LINKDEF ${PROJECT_SOURCE_DIR}/include/TXsecLinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(TXsec SHARED ${lsources} ${cdict}.cxx ${lheaders} ${cdict}.h)
target_link_libraries(TXsec ${ROOT_LIBRARIES})
add_executable(tabxsec tabxsec.cc ${sources} ${headers})
target_link_libraries(tabxsec ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -L. -lTXsec)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  tabxsec.in
  tabxsec.out
  init.mac
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()




