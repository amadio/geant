#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

# Before we use Geant4 we save all previous definitions, otherwise we will completely
# pollute the compilation environment

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
# Difficult call... sanitize compiler flags in case G4 has been built with Clang 
  STRING(REPLACE "-Qunused-arguments" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

include_directories(inc ${PROJECT_SOURCE_DIR}/xsec/inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources tabxsec.cc src/*.cc src/*.c src/*.cxx)
file(GLOB headers inc/*.hh src/*.h)

get_property(G4_COMPILE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "-std=c++11 ${CMAKE_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "-std=c++11 ${CMAKE_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "-std=c++11 ${CMAKE_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "-std=c++11 ${CMAKE_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++11 ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "-std=c++11 ${CMAKE_EXE_LINKER_FLAGS}")

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
add_executable(tabxsec tabxsec.cc ${sources} ${headers})
target_link_libraries(tabxsec ${CALLOC_LIBRARY} ${Geant4_LIBRARIES}
${ROOT_LIBRARIES} -lEG ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${XSEC_LINK_LIBS} -lGeant_v)
add_dependencies(tabxsec xseclib)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build tabxsec. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(_macros ${CMAKE_CURRENT_SOURCE_DIR}/macros)

set(TABXSEC_SCRIPTS
  ${_macros}/tabxsec.in
  ${_macros}/tabxsec.out
  ${_macros}/init.mac
  ${_macros}/init_vis.mac
  ${_macros}/run1.mac
  ${_macros}/run2.mac
  ${_macros}/vis.mac
  )

foreach(_script ${TABXSEC_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ExtracXS DEPENDS tabxsec)
install(TARGETS tabxsec DESTINATION bin)
