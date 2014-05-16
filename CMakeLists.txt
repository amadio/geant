#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GEANTV)

#----------------------------------------------------------------------------
# Additional options

option(CUDA "Enable compilation for CUDA." OFF)

#---Set pathes where to put the libraries, executables and headers------------------------------
if( GEANT_OUTPUT_DIR )
    set(OUTDIR ${GEANT_OUTPUT_DIR})
else()
    set(OUTDIR $ENV{GEANT_OUTPUT_DIR})
endif()
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
# Find VecGeom geometry headers library 
#
find_package(VecGeom REQUIRED)
if( VECGEOM_FOUND )
    find_package(Vc)
    if ( VC_FOUND )
    include_directories(${VC_INCLUDE_DIR})
    endif()
if( USE_VECGEOM_NAVIGATOR )
    add_definitions( "-DUSE_VECGEOM_NAVIGATOR" )
    set(GEANT5_CXX11_ENABLE "ON")
endif() 
   add_definitions( "-DVECGEOM_ROOT" )
   include_directories(${VECGEOM_INCLUDE_DIR})
endif()

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
set(CMAKE_CXX_FLAGS_DEBUG               "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-2")
set(GEANT5_CXX_NOW_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
set(GEANT5_CXX_NOW_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(GEANT5_CXX_NOW_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
set(GEANT5_CXX_NOW_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(GEANT5_EXE_NOW_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS}")
set(GEANT5_WARNING_FLAGS "-W -Wall -Weffc++ -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")


if(GEANT5_CXX11_ENABLE)
   set(CXXSTD "-std=c++11")
else()
   set(CXXSTD "")
endif()

set(GEANT5_CXX_FLAGS                "${CMAKE_CXX_FLAGS} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_COMPILE_FLAGS            "${GEANT5_COMPILE_FLAGS_NOW} ${CXXSTD} ${GEANT5_WARNING_FLAGS}")

#---------------------------------------------------------------------------
if(CUSTOM_ALLOCATOR)
  find_package(CAllocator)
endif()
if (CUDA)
  find_package(CUDA REQUIRED)
  include(FindCUDA)

  if (NOT CUDA_ARCH)
    set(CUDA_ARCH 35)
  endif()
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")

  #set(CUDA_NVCC_FLAGS --ptxas-options=-v)
  set(
    CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}
    -O2 --use_fast_math -Xcudafe "--diag_suppress=code_is_unreachable"
    -Xcudafe "--diag_suppress=initialization_not_reachable"
  )
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -g -G)
  endif()
  # set(CUDA_SEPARABLE_COMPILATION ON)
  # set(CUDA_VERBOSE_BUILD ON)

  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")

endif()

#---------------------------------------------------------------------------
# Include directories
include_directories("${PROJECT_SOURCE_DIR}/vecprot_v2/inc")

#----------------------------------------------------------------------------
# Add subdirectories
#
add_custom_target(_headers ALL)

add_subdirectory (xsec)
add_subdirectory (vecprot_v2)
add_subdirectory (tabxsec)
add_subdirectory (micromc)
add_subdirectory (test-small)

if (CUDA)
  add_subdirectory (vecprot_v2/cuda)
endif()

