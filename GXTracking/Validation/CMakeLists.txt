#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
Project(Validation)

#----------------------------------------------------------------------------
# setup #define GPUDEBUG 
#
option(GPUPLOTS      "Turn GPU validation plots ON/OFF" ON)
option(GPUDEBUG      "Turn GPU validation printout ON/OFF" OFF)
option(GPUNONRANDOM  "Force the use of NonRandomEngine ON/OFF" OFF)

#.. Assume user always wants validation histograms
if(DEFINED GPUNONRANDOM AND GPUNONRANDOM STREQUAL "ON")
    #set_target_properties(G4validation PROPERTIES COMPILE_DEFINITIONS "GPUDEBUG=1 GPUNONRANDOM=1")
    add_definitions(-DGPUNONRANDOM)
endif(DEFINED GPUNONRANDOM AND GPUNONRANDOM STREQUAL "ON")

if(DEFINED GPUDEBUG AND GPUDEBUG STREQUAL "ON")
  #set_target_properties(G4validation PROPERTIES COMPILE_DEFINITIONS "GPUDEBUG=1")
  add_definitions(-DGPUDEBUG)
endif(DEFINED GPUDEBUG AND GPUDEBUG STREQUAL "ON")

if(DEFINED GPUPLOTS AND GPUPLOTS STREQUAL "ON")
  add_definitions(-DGPUPLOTS)
endif(DEFINED GPUPLOTS AND GPUPLOTS STREQUAL "ON")

#### Setting up RPATH for both build-time and installed binaries.  ###########
#### See http://www.cmake.org/Wiki/CMake_RPATH_handling for details.  ########
#
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
#
#### end of RPATH-related commands  #####################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../cmake/Modules)
#message(STATUS "From Validation: CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

find_package(ROOT REQUIRED)
#message(STATUS "From Validation: ROOT_INCLUDE_DIR = " ${ROOT_INCLUDE_DIR})

### Process subdirectories' CMakeLists.txt
add_subdirectory(utils)
add_subdirectory(g4app)
add_subdirectory(prototype)
