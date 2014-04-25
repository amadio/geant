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

if(${USE_ROOT})
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ../cmake/modules)
  find_package(ROOT REQUIRED)
  add_definitions(-DUSE_ROOT)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                    ${PROJECT_SOURCE_DIR}/../xsec/inc 
)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

get_property(G4_COMPILE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${CMAKE_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${CMAKE_EXE_LINKER_FLAGS}")

if (Geant4_VERSION_MAJOR LESS 10)
   set_source_files_properties(${sources} PROPERTY COMPILE_DEFINITIONS            "${G4_COMPILE_DEFINITIONS}")
else()
   set_source_files_properties(${sources} PROPERTY COMPILE_DEFINITIONS            "${G4_COMPILE_DEFINITIONS} -DG4MULTITHREADED")
endif()

set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${GEANT5_COMPILE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${GEANT5_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANT5_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_NOW_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${GEANT5_EXE_NOW_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleN03 exampleN03.cc ${sources} ${headers})
target_link_libraries(exampleN03 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} 
                                 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v -lGeom)
add_dependencies(exampleN03 xseclib Geantlib_v)

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

