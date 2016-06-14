#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6.3 FATAL_ERROR)
project(G01)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED gdml vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
find_package(ROOT REQUIRED)

include_directories(AFTER SYSTEM ${ROOT_INCLUDE_DIR})

if( USE_ROOT )
   add_definitions( "-DUSE_ROOT")
endif()

if( SNEXTG4 )
  add_definitions("-DSNEXTG4")
endif()

if ( NAVTEST ) 
  add_definitions("-DNAVTEST")
endif()

#----------------------------------------------------------------------------
# Find Pythia8 and setup the libraries
#

find_package(Pythia8 REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Define compiler flags
#

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(load_gdml load_gdml.cc ${sources} ${headers})
target_link_libraries(load_gdml ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lGeom -lEG -lTreePlayer -lGui
  ${PYTHIA8_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G01_SCRIPTS
    g01.in event.in vis.mac cms.gdml gungeantino.in materials.xml scanevent.C vis.mac)

foreach(_script ${G01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(G01 DEPENDS load_gdml)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS load_gdml DESTINATION bin)

