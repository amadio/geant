#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TestEm5)

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
set(OUTPUT bin/examples/TestEm5/Geant4)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})
add_executable(TestEm5 TestEm5.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(TestEm5 ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#
set(TestEm5_SCRIPTS
    hanson19_G4.mac
    README
   )

 foreach(_script ${TestEm5_SCRIPTS})
   configure_file(
     ${_script}
      ${CMAKE_INSTALL_PREFIX}/${OUTPUT}/${_script}
     COPYONLY
     )
#     configure_file(
#       ${_script}
#       ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${OUTPUT}/${_script}
#       COPYONLY
#     )
 endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestEm5 DESTINATION ${OUTPUT})
