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
include_directories(inc ${PROJECT_SOURCE_DIR}/xsec/inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cc)
file(GLOB headers inc/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(tabxsec tabxsec.cc ${sources} ${headers})
target_link_libraries(tabxsec ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lEg -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(_macros ${CMAKE_CURRENT_SOURCE_DIR}/macros)

set(TABXSEC_SCRIPTS
  ${_macros}/tabxsec.in
  ${_macros}/tabxsec.out
  ${_macros}/init.mac
  ${_macros}/init_vis.mac
  ${_macros}/run1.mac
  ${_macros}/run2.mac
  ${_macros}/vis.mac
  )

foreach(_script ${TABXSEC_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ExtracXS DEPENDS tabxsec)
