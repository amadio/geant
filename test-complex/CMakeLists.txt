#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TESTCOMPLEX)

#----------------------------------------------------------------------------
# Find Geant4 package, batch mode only
#
find_package(Geant4 REQUIRED gdml)

#----------------------------------------------------------------------------
# Find HepMC (required package)
#
find_package(HepMC REQUIRED)
if(HepMC_FOUND)
  message(STATUS "HepMC includes: ${HEPMC_INCLUDE_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
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
# Define compiler flags : already done in the main CMakeLists
#
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(testcomplex testcomplex.cc ${sources} ${headers})
target_link_libraries(testcomplex ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${XSEC_LINK_LIBS} -lGeom -lEG )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TESTCOMPLEX. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(TESTCOMPLEX_SCRIPTS
    g4macro_TPHYS.mac g4macro_G4.mac testcomplex.sh 
   )

foreach(_script ${TESTCOMPLEX_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    #${PROJECT_BINARY_DIR}/${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
#
add_custom_target(TESTCOMPLEX DEPENDS testcomplex)

#----------------------------------------------------------------------------
# Setting Geant4 version compiler definition  
#
if (Geant4_VERSION_MAJOR LESS 10)
  set_target_properties(testcomplex PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=1") 
else()
  set_target_properties(testcomplex PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=0") 
endif()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS testcomplex DESTINATION bin)

