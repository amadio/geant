#-----------------------------------------------------------------------------
# - Define CMake requirements and Project
#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8)
project(${Project})

#------------------------------------------------------------------------------
# Configuration options
#------------------------------------------------------------------------------

option(VC_ACCELERATION "Internal vectorization using Vc." ON)

option(CUDA "Enable compilation for CUDA." OFF)

set(VECPHYS_VECTOR "" CACHE STRING "Vector instruction set to use (sse/avx/..)")

if (NOT BACKEND)
  message(STATUS "Backend not specified. Defaulting to Vc.")
  set(BACKEND "Vc")
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

#------------------------------------------------------------------------------
# Set C++ compiler flags. Will also be forwarded by CUDA when compiling C++.
#------------------------------------------------------------------------------

set(VECPHYS_COMPILATION_FLAGS)
set(VECPHYS_COMPILATION_FLAGS_GNU)
set(VECPHYS_COMPILATION_FLAGS_CLANG)
set(VECPHYS_COMPILATION_FLAGS_INTEL)
set(VECPHYS_NVCC_COMPILATION_FLAGS)

set(VECPHYS_COMPILATION_FLAGS
    "${VECPHYS_COMPILATION_FLAGS} -Wall -fPIC -DVECPHYS")
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

set(VECPHYS_COMPILATION_FLAGS_GNU
    "${VECPHYS_COMPILATION_FLAGS_GNU} -finline-limit=10000000")
set(VECPHYS_COMPILATION_FLAGS_GNU
    "${VECPHYS_COMPILATION_FLAGS_GNU} -ffast-math -ftree-vectorize")
set(VECPHYS_COMPILATION_FLAGS_CLANG
    "${VECPHYS_COMPILATION_FLAGS_CLANG} -ffast-math -ftree-vectorize")
set(VECPHYS_COMPILATION_FLAGS_INTEL
    "${VECPHYS_COMPILATION_FLAGS_INTEL} -xHost")
if (VECPHYS_VECTOR)
  message(STATUS "Compiling with vector instruction set ${VECPHYS_VECTOR}")
  set(VECPHYS_COMPILATION_FLAGS_GNU
      "${VECPHYS_COMPILATION_FLAGS_GNU} -m${VECPHYS_VECTOR}")
  set(VECPHYS_COMPILATION_FLAGS_CLANG
      "${VECPHYS_COMPILATION_FLAGS_CLANG} -m${VECPHYS_VECTOR}")
endif()

#------------------------------------------------------------------------------
# Backend specific
#------------------------------------------------------------------------------

if (Scalar)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_SCALAR")
endif()
if (Vc)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_VC")
  set(VECPHYS_COMPILATION_FLAGS_GNU
      "${VECPHYS_COMPILATION_FLAGS_GNU} -fabi-version=6")
endif()
if (Cilk)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CILK")
endif()
if (CUDA)
  set(CUDA_ARCH 20 CACHE STRING "CUDA device architecture.")
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CUDA")
  if (NOT CUDA_SPECIALIZATION)
    set(VECPHYS_COMPILATION_FLAGS
        "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CUDA_NO_SPECIALIZATION")
  endif()
  set(VECPHYS_NVCC_COMPILATION_FLAGS ${VECPHYS_NVCC_COMPILATION_FLAGS}
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
  set(VECPHYS_COMPILATION_FLAGS
      "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_VC_ACCELERATION")
endif()

if (ROOT)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_ROOT")
endif()

#------------------------------------------------------------------------------
# Pass flags to compilers
#------------------------------------------------------------------------------

if (GNU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECPHYS_COMPILATION_FLAGS_GNU}")
endif()
if (Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECPHYS_COMPILATION_FLAGS_CLANG}")
endif()
if (Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECPHYS_COMPILATION_FLAGS_INTEL}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECPHYS_COMPILATION_FLAGS}")
set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -O2 -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${VECPHYS_NVCC_COMPILATION_FLAGS})

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (CUDA)
  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Backend configuration
#------------------------------------------------------------------------------

if (VC_ACCELERATION OR Vc)

  find_package(Vc REQUIRED)
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} ${Vc_LIBRARIES})
  include_directories(AFTER SYSTEM ${Vc_INCLUDE_DIR})

  if (Vc)
    set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/vc/Backend.cpp)
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

if (ROOT)
  
  find_package(ROOT REQUIRED)

  include_directories(AFTER SYSTEM ${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})
  link_directories(${ROOT_LIBRARY_DIR}/root) # For natively installed ROOT
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL}
      Cint Core MathCore RIO Thread Geom Gui Net Hist Graf Graf3d Gpad Tree
      Rint Postscript Matrix Physics m dl)
endif()

#------------------------------------------------------------------------------
# Include and source
#------------------------------------------------------------------------------

include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/inc
)

set(SRC_CPP_RELATIVE
  src/GUTrackHandler.cc
  src/GUComptonKleinNishina.cc
  src/GUAliasSampler.cc
)
foreach(SRC ${SRC_CPP_RELATIVE})
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

add_library(vecphys ${SRC_CPP})
target_link_libraries(vecphys ${VECPHYS_LIBRARIES}
                      ${VECPHYS_LIBRARIES_EXTERNAL})
set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphys)
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} rt)
endif()

if (CUDA)
  cuda_add_library(
    vecphyscuda
    ${SRC_CUDA}
    OPTIONS ${CUDA_ARCH}
  )
  target_link_libraries(vecphyscuda ${VECPHYS_LIBRARIES}
                        ${VECPHYS_LIBRARIES_EXTERNAL})
  set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphyscuda)
endif()

#------------------------------------------------------------------------------
# Define executables
#------------------------------------------------------------------------------

set(EXECUTABLES
  ${CMAKE_SOURCE_DIR}/test/${Project}.cc
)

set_source_files_properties(
  ${EXECUTABLES}
  PROPERTIES COMPILE_FLAGS -std=c++11
)

add_executable(${Project} ${EXECUTABLES})
target_link_libraries(${Project} ${VECPHYS_LIBRARIES} 
                                 ${VECPHYS_LIBRARIES_EXTERNAL})

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
