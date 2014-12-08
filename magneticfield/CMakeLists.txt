#-----------------------------------------------------------------------------
# - Define CMake requirements and Project
#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8)

set(project VecMagField)
project(VecMagField)

#------------------------------------------------------------------------------
# Configuration options
#------------------------------------------------------------------------------

option(VC_ACCELERATION "Internal vectorization using Vc." OFF)  ## Not available Dec 2014
option(CUDA "Enable compilation for CUDA." OFF)

set(VECFIELD_VECTOR "" CACHE STRING "Vector instruction set to use (sse/avx/..)")

if (NOT BACKEND)
  message(STATUS "Backend not specified. Defaulting to Scalar.")
  set(BACKEND "Scalar")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
                      ${PROJECT_SOURCE_DIR}/cmake/modules)

#------------------------------------------------------------------------------
# Determine backend
#------------------------------------------------------------------------------

set(Vc)
set(Cilk)
set(Scalar)

if (BACKEND STREQUAL "Vc")
  set(Vc TRUE)
endif()
if (BACKEND STREQUAL "Cilk")
  set(Cilk TRUE)
endif()
if (BACKEND STREQUAL "Scalar")
  set(Scalar TRUE)
endif()
if (NOT Vc AND NOT Cilk AND NOT Scalar)
  message(FATAL_ERROR "Invalid backend specified: \"${BACKEND}\"")
endif()

message(STATUS "Configuring with vector backend ${BACKEND}.")

#------------------------------------------------------------------------------
# Determine C++ compiler being used
#------------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(Clang TRUE)
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(GNU TRUE)
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(Intel TRUE)
endif()
if (NOT GNU AND NOT Clang AND NOT Intel)
  message(WARNING "Unsupported compiler. Build will likely fail.")
endif()

## Flags segment copied from VECPHYS CMakeLists.txt -- which was derived from VecGeom

#------------------------------------------------------------------------------
# Set C++ compiler flags. Will also be forwarded by CUDA when compiling C++.
#------------------------------------------------------------------------------

set(VECFIELD_COMPILATION_FLAGS)
set(VECFIELD_COMPILATION_FLAGS_GNU)
set(VECFIELD_COMPILATION_FLAGS_CLANG)
set(VECFIELD_COMPILATION_FLAGS_INTEL)
set(VECFIELD_NVCC_COMPILATION_FLAGS)

set(VECFIELD_COMPILATION_FLAGS
    "${VECFIELD_COMPILATION_FLAGS} -Wall -fPIC -DVECFIELD")
if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++ -ggdb")
  if (APPLE OR UNIX) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

#------------------------------------------------------------------------------
# Compiler specific optimization flags
#------------------------------------------------------------------------------

set(VECFIELD_COMPILATION_FLAGS_GNU
    "${VECFIELD_COMPILATION_FLAGS_GNU} -finline-limit=10000000")
set(VECFIELD_COMPILATION_FLAGS_GNU
    "${VECFIELD_COMPILATION_FLAGS_GNU} -ffast-math -ftree-vectorize")
set(VECFIELD_COMPILATION_FLAGS_CLANG
    "${VECFIELD_COMPILATION_FLAGS_CLANG} -ffast-math -ftree-vectorize")
set(VECFIELD_COMPILATION_FLAGS_INTEL
    "${VECFIELD_COMPILATION_FLAGS_INTEL} -xHost")
if (VECFIELD_VECTOR)
  message(STATUS "Compiling with vector instruction set ${VECFIELD_VECTOR}")
  set(VECFIELD_COMPILATION_FLAGS_GNU
      "${VECFIELD_COMPILATION_FLAGS_GNU} -m${VECFIELD_VECTOR}")
  set(VECFIELD_COMPILATION_FLAGS_CLANG
      "${VECFIELD_COMPILATION_FLAGS_CLANG} -m${VECFIELD_VECTOR}")
endif()

#------------------------------------------------------------------------------
# Backend specific
#------------------------------------------------------------------------------

if (Scalar)
  set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_SCALAR")
endif()
if (Vc)
  set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_VC")
  set(VECFIELD_COMPILATION_FLAGS_GNU
      "${VECFIELD_COMPILATION_FLAGS_GNU} -fabi-version=6")
endif()
if (Cilk)
  set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_CILK")
endif()
if (CUDA)
  set(CUDA_ARCH 20 CACHE STRING "CUDA device architecture.")
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")
  set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_CUDA")
  if (NOT CUDA_SPECIALIZATION)
    set(VECFIELD_COMPILATION_FLAGS
        "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_CUDA_NO_SPECIALIZATION")
  endif()
  set(VECFIELD_NVCC_COMPILATION_FLAGS ${VECFIELD_NVCC_COMPILATION_FLAGS}
      -O2 --use_fast_math -Xcudafe "--diag_suppress=code_is_unreachable"
      -Xcudafe "--diag_suppress=initialization_not_reachable")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -g -G)
  endif()
endif()

#------------------------------------------------------------------------------
# Forward CMake options
#------------------------------------------------------------------------------

if (VC_ACCELERATION)
  set(VECFIELD_COMPILATION_FLAGS
      "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_VC_ACCELERATION")
endif()

if (ROOT)
  set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -DVECFIELD_ROOT")
endif()

#------------------------------------------------------------------------------
# Pass flags to compilers
#------------------------------------------------------------------------------

if (GNU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECFIELD_COMPILATION_FLAGS_GNU}")
endif()
if (Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECFIELD_COMPILATION_FLAGS_CLANG}")
endif()
if (Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECFIELD_COMPILATION_FLAGS_INTEL}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECFIELD_COMPILATION_FLAGS}")
set(VECFIELD_COMPILATION_FLAGS "${VECFIELD_COMPILATION_FLAGS} -O2 -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${VECFIELD_NVCC_COMPILATION_FLAGS})

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (CUDA)
  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Backend configuration
#------------------------------------------------------------------------------

find_package(VecGeom REQUIRED)
if( VECGEOM_FOUND )
  include_directories(${VECGEOM_INCLUDE_DIR})

  if (VC_ACCELERATION OR Vc)

  find_package(Vc REQUIRED)
  set(VECFIELD_LIBRARIES_EXTERNAL ${VECFIELD_LIBRARIES_EXTERNAL} ${Vc_LIBRARIES})
  include_directories(AFTER SYSTEM ${Vc_INCLUDE_DIR})

  if (Vc)
    set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/vc/Backend.cpp)
  endif()

  endif()
endif()

if (Cilk)
  if (NOT Intel)
    message(FATAL_ERROR "Must use Intel C++ compiler (icc) for Cilk backend.")
  endif()
  set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/cilk/Backend.cpp)
endif()

if (CUDA)
  find_package(CUDA REQUIRED)
  include(FindCUDA)

  set(CUDA_SEPARABLE_COMPILATION ON)

  set(SRC_EXTERNAL ${SRC_EXTERNAL} source/CudaManager.cpp)
  set(SRC_CUDA ${SRC_CUDA}
      source/CudaManager.cu source/backend/cuda/Interface.cu)
endif()

#------------------------------------------------------------------------------
# Include and source
#------------------------------------------------------------------------------

include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/inc
)

set(HEADERS 
    inc/ConstFieldHelixStepper.h
    inc/GVIntegratorDriver.h
    inc/TClassicalRK4.h
    inc/TIntegratorDriver.h
    inc/TMagFieldEquation.h
    inc/TSimpleRunge.h
    inc/TUniformMagField.h
)

set(SOURCES
  src/FieldSetup.cc
  src/circlestepper.cxx
  src/TestConstBzStepper.cxx
)

foreach(SRC ${SOURCES})
  set(SRC_CPP ${SRC_CPP} ${CMAKE_SOURCE_DIR}/${SRC})
endforeach()

# Copy all source files to .cu-files in order for NVCC to compile them as CUDA
# code and not regular C++ files.

if (CUDA)

  foreach(SRC_FILE ${SRC_CPP})

    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME_WE)

    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC_FILE}
          ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      DEPENDS ${SRC_FILE}
    )

    set(
      SRC_CUDA ${SRC_CUDA}
      ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
    )

  endforeach()

endif()

# Compile cpp-files for C++11

set_source_files_properties(${SRC_CPP} PROPERTIES COMPILE_FLAGS -std=c++11)

#------------------------------------------------------------------------------
# Build libraries
#------------------------------------------------------------------------------

add_library(magfield ${SRC_CPP})
target_link_libraries(magfield ${VECFIELD_LIBRARIES}
                      ${VECFIELD_LIBRARIES_EXTERNAL})
set(VECFIELD_LIBRARIES ${VECFIELD_LIBRARIES} magfield)
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(VECFIELD_LIBRARIES_EXTERNAL ${VECFIELD_LIBRARIES_EXTERNAL} rt)
endif()

if (CUDA)
  cuda_add_library(
    magfieldcuda
    ${SRC_CUDA}
    OPTIONS ${CUDA_ARCH}
  )
  target_link_libraries(magfieldcuda ${VECFIELD_LIBRARIES}
                        ${VECFIELD_LIBRARIES_EXTERNAL})
  set(VECFIELD_LIBRARIES ${VECFIELD_LIBRARIES} magfieldcuda)
endif()

#------------------------------------------------------------------------------
# Define executables
#------------------------------------------------------------------------------

set(EXECUTABLES
  ${CMAKE_SOURCE_DIR}/test/TestVecMagField
)

set_source_files_properties(
  ${EXECUTABLES}
  PROPERTIES COMPILE_FLAGS -std=c++11
)

message(STATUS "Project name: " ${Project})
message(STATUS "Executables = " ${EXECUTABLES})
message(STATUS "Sources     = " ${SOURCES})
message(STATUS "Headers     = " ${HEADERS})

add_executable(${Project} ${EXECUTABLES} ${${SOURCES}} ${HEADERS})


add_executable(${Project} ${EXECUTABLES} ${${SOURCES}} ${HEADERS})
target_link_libraries(${Project} ${VECFIELD_LIBRARIES} 
                                 ${VECFIELD_LIBRARIES_EXTERNAL})

# add_executable(load_gdml load_gdml.cc ${SOURCES} ${headers})
#   target_link_libraries(load_gdml ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lGeom
#   ${PYTHIA8_LIBRARIES})

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
