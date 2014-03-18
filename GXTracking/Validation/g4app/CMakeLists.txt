#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4validation)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

#message( STATUS "From geant4 dir: CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
#message( STATUS "From geant4 dir: CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR} )
#message( STATUS "From geant4 dir: PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#find_package(ROOT REQUIRED)
message(STATUS "From geant4 CMakeList: ROOT_INCLUDE_DIR=" ${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# setup #define GPUDEBUG 
#
#option(GPUDEBUG      "Turn GPU validation printout ON/OFF" OFF)
#option(GPUNONRANDOM  "Force the use of NonRandomEngine ON/OFF" OFF)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
#message(STATUS "From geant4: UTILS_SOURCE_DIR=" ${UTILS_SOURCE_DIR})
include(${Geant4_USE_FILE} ${UTILS_SOURCE_DIR}/include)

#if(${USE_ROOT})
#  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /scratch/lima/geant/geant4.9.6.p01-orig/cmake/modules)
#  find_package(ROOT REQUIRED)
#endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${CMAKE_SOURCE_DIR}/utils/include
		    ${ROOT_INCLUDE_DIR}
)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#add_library( STATIC ${sources})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#add_executable(G4validation G4validation.cc ${sources} ${headers})
add_executable(G4validation G4validation.cc ${sources})

#message(STATUS "From geant4 CMakeLists: CMAKE_BINARY_DIR=" ${CMAKE_BINARY_DIR})
target_link_libraries(G4validation
   ${Geant4_LIBRARIES}
   -L${CMAKE_BINARY_DIR}/utils -lutils
   ${ROOT_LIBRARIES} )

 
#if(DEFINED GPUNONRANDOM AND GPUNONRANDOM STREQUAL "ON")
#    #set_target_properties(G4validation PROPERTIES COMPILE_DEFINITIONS "GPUDEBUG=1 GPUNONRANDOM=1")
#    add_definitions(-DGPUNONRANDOM)
#endif(DEFINED GPUNONRANDOM AND GPUNONRANDOM STREQUAL "ON")
#
#if(DEFINED GPUDEBUG AND GPUDEBUG STREQUAL "ON")
#  #set_target_properties(G4validation PROPERTIES COMPILE_DEFINITIONS "GPUDEBUG=1")
#  add_definitions(-DGPUDEBUG)
#endif(DEFINED GPUDEBUG AND GPUDEBUG STREQUAL "ON")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G4validation. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4app_SCRIPTS
    macros/G4validation.in macros/G4validation.large_N.in macros/exo1.mac macros/exo2.mac macros/gui.mac macros/newgeom.mac macros/rndmSeed.mac macros/run1.mac macros/run2.mac macros/vis.mac macros/icons.mac
  )

foreach(_script ${G4app_SCRIPTS})
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
add_custom_target(G4app DEPENDS G4validation)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4validation DESTINATION bin)
