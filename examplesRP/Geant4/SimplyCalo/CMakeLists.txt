#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SimplyCalo)

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
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(OUTPUT bin/examples/Calorimeter/Geant4)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})
add_executable(SimplyCalo ${CMAKE_SOURCE_DIR}/examplesRP/Geant4/SimplyCalo/SimplyCalo.cc ${sources} ${headers})
target_link_libraries(SimplyCalo ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#
set(SimplyCalo_SCRIPTS
    SimplyCalo_ExN03.mac
    SimplyCalo_ATLASbar.mac
    README
   )

 foreach(_script ${SimplyCalo_SCRIPTS})
   configure_file(
     ${_script}
      ${CMAKE_INSTALL_PREFIX}/${OUTPUT}/${_script}
     COPYONLY
     )
#     configure_file(
#       ${_script}
#       ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
#       COPYONLY
#     )
 endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SimplyCalo DESTINATION ${OUTPUT})
