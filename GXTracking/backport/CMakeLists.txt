#------------------------------------------------------------------------------
# - CMakeLists.txt for tracking
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# - Project definition
#------------------------------------------------------------------------------
project(tracking)

#------------------------------------------------------------------------------
# - Find the package named 'Geant4'
#------------------------------------------------------------------------------
find_package(Geant4 REQUIRED)

#------------------------------------------------------------------------------
# - Find the package named 'Threads'
#------------------------------------------------------------------------------
#find_package(Threads REQUIRED)

#------------------------------------------------------------------------------
# - Include the file listed in the Geant4_USE_FILE variable
#------------------------------------------------------------------------------
include(${Geant4_USE_FILE})

#------------------------------------------------------------------------------
# - collates all the sources and headers in and below the current directory
#------------------------------------------------------------------------------
GEANT4_COLLATE_APPLICATION_SOURCES(TRACKING_SOURCES)

add_executable(tracking ${TRACKING_SOURCES})
target_link_libraries(tracking ${Geant4_LIBRARIES} pthread dl )

#
# Set CXXFLAGS
#
set(CMAKE_LD_FLAGS -O2)
set(CMAKE_CXX_FLAGS -O2)

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
