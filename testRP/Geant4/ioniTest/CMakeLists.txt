#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ioniTest_G4)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
#
set(OUTPUT bin/testRP/ioniTest/Geant4)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})
add_executable(ioniTest_G4 ioniTest_G4.cc ${sources} ${headers})
target_link_libraries(ioniTest_G4 ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#
set(ioniTest_SCRIPTS
    README
   )

foreach(_script ${ioniTest_SCRIPTS})
  configure_file(
    ${_script}
     ${CMAKE_INSTALL_PREFIX}/${OUTPUT}/${_script}
    COPYONLY
    )
    configure_file(
      ${_script}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
      COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ioniTest_G4 DESTINATION ${OUTPUT})
