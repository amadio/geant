#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

# Before we use Geant4 we save all previous definitions, otherwise we will completely
# pollute the compilation environment

get_property(SAVE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
set(SAVE_CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}")
set(SAVE_CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
set(SAVE_CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(SAVE_CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
set(SAVE_CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(SAVE_CMAKE_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS}")

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

get_property(G4_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_CXX_FLAGS                ${CMAKE_CXX_FLAGS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_CXX_FLAGS_DEBUG          ${CMAKE_CXX_FLAGS_DEBUG})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    ${CMAKE_CXX_FLAGS_MINISIZEREL})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_CXX_FLAGS_RELEASE        ${CMAKE_CXX_FLAGS_RELEASE})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY CMAKE_EXEC_LINKER_FLAGS        ${CMAKE_EXE_LINKER_FLAGS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS            ${G4_DEFINITIONS})

set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${SAVE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${SAVE_CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${SAVE_CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${SAVE_CMAKE_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${SAVE_CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${SAVE_CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${SAVE_CMAKE_EXE_LINKER_FLAGS}")

include_directories(inc ${PROJECT_SOURCE_DIR}/xsec/inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cc)
file(GLOB headers inc/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(tabxsec tabxsec.cc ${sources} ${headers})
target_link_libraries(tabxsec ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lEG -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
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
