#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N03)

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

if(USE_ROOT)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ../cmake/modules)
  find_package(ROOT REQUIRED)
  add_definitions(-DUSE_ROOT)
endif()

option(MAKESTAT "Collect detailed statistics during exampleN03 run!" OFF)
if(MAKESTAT)
  add_definitions(-DMAKESTAT)
  message(STATUS "MAKESTAT ON  : Detailed statistics IS collected during exampleN03 run!")
else()
  message(STATUS "MAKESTAT OFF : Detailed statistics IS NOT collected during exampleN03 run!")
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
message(STATUS "Using Geant4 file ${Geant4_USE_FILE}")
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                    ${PROJECT_SOURCE_DIR}/../xsec/inc 
)

if(CMAKE_BUILD_TYPE MATCHES None)
  set(_cf ${CMAKE_CXX_FLAGS})
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
   set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES Release)
   set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
   set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
elseif(CMAKE_BUILD_TYPE MATCHES MiniSizeRel)
   set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINISIZEREL})
endif()

separate_arguments(_cf)
set_directory_properties(PROPERTIES COMPILE_OPTIONS "${_cf}")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleN03 exampleN03.cc ${sources} ${headers})
set_target_properties(tabxsec PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(exampleN03 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} 
                                 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${XSEC_LINK_LIBS} -lGeant_v -lGeom -lEG -lTreePlayer -lGui)
if(USE_VECGEOM_NAVIGATOR)
   target_link_libraries(exampleN03 ${VECGEOM_LIBRARIES})
endif()
add_dependencies(exampleN03 xseclib Geantlib_v)

set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(N03_SCRIPTS
    exampleN03.in exampleN03.large_N.in exampleN03.out exo1.mac exo2.mac gui.mac newgeom.mac rndmSeed.mac run1.mac run2.mac vis.mac icons.mac
    visTutor/exN03Tree0.mac visTutor/exN03Vis12.loop visTutor/exN03Vis12.mac visTutor/exN03Vis13.loop visTutor/exN03Vis13.mac visTutor/exN03Vis5.mac visTutor/g4Xt.xrm visTutor/heprep.mac
  )

foreach(_script ${N03_SCRIPTS})
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
add_custom_target(N03 DEPENDS exampleN03)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleN03 DESTINATION bin)

