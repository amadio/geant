#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GEANTV)


#---Set pathes where to put the libraries, executables and headers------------------------------
set(OUTDIR $ENV{GEANT_OUTPUT_DIR})
if(NOT OUTDIR)
  set(OUTDIR ${PROJECT_SOURCE_DIR})
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}/bin)

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Load macros
#
include(cmake/modules/GeantMacros.cmake)

#----------------------------------------------------------------------------
# Define compiler flags
#
get_property(GEANT5_COMPILE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
get_property(GEANT5_COMPILE_FLAGS_NOW GLOBAL PROPERTY COMPILE_FLAGS)
set(GEANT5_CXX_NOW_FLAGS                "${CMAKE_CXX_FLAGS}")
set(GEANT5_CXX_NOW_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
set(GEANT5_CXX_NOW_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(GEANT5_CXX_NOW_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
set(GEANT5_CXX_NOW_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(GEANT5_EXE_NOW_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS}")
set(GEANT5_WARNING_FLAGS "-W -Wall -Weffc++ -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")

set(GEANT5_CXX_FLAGS                "${CMAKE_CXX_FLAGS} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_COMPILE_FLAGS            "${GEANT5_COMPILE_FLAGS_NOW} ${GEANT5_WARNING_FLAGS}")

#---------------------------------------------------------------------------
if(CUSTOM_ALLOCATOR)
  find_package(CAllocator)
endif()

#----------------------------------------------------------------------------
# Add subdirectories
#
add_custom_target(_headers ALL)

add_subdirectory (xsec)
add_subdirectory (tabxsec)
add_subdirectory (vecprot)
add_subdirectory (micromc)
add_subdirectory (geom_vec_tests)

