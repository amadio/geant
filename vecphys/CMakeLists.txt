cmake_minimum_required(VERSION 2.8.8)

################################################################################

# Configuration options

option(VC_ACCELERATION "Internal vectorization using Vc." ON)

option(MIC "Enable compilation for MIC." OFF)

option(GUBENCHMARK "Enable performance comparisons." OFF)

option(DISTANCE_DEBUG "Enable comparision of calculcated distances againt ROOT/Geant4 behind the scenes" OFF)

option(VECPHYS_STANDALONE "Build VECPHYS standalone" OFF)

set(VECPHYS_VECTOR "" CACHE STRING
    "Vector instruction set to use (sse/avx/...).")

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries.")
set(INSTALL_INCLUDE_DIR include CACHE PATH
    "Installation directory for header files.")
set(INSTALL_CMAKE_DIR lib/CMake/VecPhys CACHE PATH
    "Installation directory for CMake files.")

# Make relative paths absolute
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

if (NOT BACKEND)
  message(STATUS "Backend not specified. Defaulting to Vc.")
  set(BACKEND "Vc")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

################################################################################
SET(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda-7.5")

# Determine backend

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


# Determine C++ compiler being used

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(Clang TRUE)
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(GNU TRUE)
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(Intel TRUE)
endif()
if (NOT GNU AND NOT Clang AND NOT Intel)
  message(WARNING "Unsupported compiler. Build will likely fail.")
endif()

################################################################################

# Set C++ compiler flags. Will also be forwarded by CUDA when compiling C++.

set(VECPHYS_COMPILATION_FLAGS)
set(VECPHYS_COMPILATION_FLAGS_GNU)
set(VECPHYS_COMPILATION_FLAGS_CLANG)
set(VECPHYS_COMPILATION_FLAGS_INTEL)
set(VECPHYS_NVCC_COMPILATION_FLAGS)

set(VECPHYS_COMPILATION_FLAGS
    "${VECPHYS_COMPILATION_FLAGS} -Wall -Wextra -fPIC")
if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++ -ggdb")
  if (APPLE OR UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

# Compiler specific optimization flags
if (NOT CMAKE_BUILD_TYPE MATCHES Debug)
set(VECPHYS_COMPILATION_FLAGS_GNU
    "${VECPHYS_COMPILATION_FLAGS_GNU} -finline-limit=10000000")
set(VECPHYS_COMPILATION_FLAGS_GNU
    "${VECPHYS_COMPILATION_FLAGS_GNU} -ffast-math -ftree-vectorize")
set(VECPHYS_COMPILATION_FLAGS_CLANG
    "${VECPHYS_COMPILATION_FLAGS_CLANG} -ffast-math -ftree-vectorize")
endif()
set(VECPHYS_COMPILATION_FLAGS_INTEL
    "${VECPHYS_COMPILATION_FLAGS_INTEL} ")
if (VECPHYS_VECTOR)
  message(STATUS "Compiling with vector instruction set ${VECPHYS_VECTOR}")
  set(VECPHYS_COMPILATION_FLAGS_GNU
      "${VECPHYS_COMPILATION_FLAGS_GNU} -m${VECPHYS_VECTOR}")
  set(VECPHYS_COMPILATION_FLAGS_CLANG
      "${VECPHYS_COMPILATION_FLAGS_CLANG} -m${VECPHYS_VECTOR}")
endif()

# Backend specific
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
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CUDA")
  if (CUDA_VOLUME_SPECIALIZATION)
    set(VECPHYS_COMPILATION_FLAGS
        "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CUDA_VOLUME_SPECIALIZATION")
  else()
    set(VECPHYS_COMPILATION_FLAGS
        "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_CUDA_NO_VOLUME_SPECIALIZATION")
  endif()
  set(VECPHYS_NVCC_COMPILATION_FLAGS ${VECPHYS_NVCC_COMPILATION_FLAGS}
      -Xcudafe "--diag_suppress=code_is_unreachable"
      -Xcudafe "--diag_suppress=initialization_not_reachable" 
      --expt-relaxed-constexpr
  )
  if(VECPHYS_STANDALONE) 
    set(CUDA_ARCH 35 CACHE STRING "CUDA device architecture.")
    set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")
    set(VECPHYS_NVCC_COMPILATION_FLAGS ${VECPHYS_NVCC_COMPILATION_FLAGS}
       -std=c++11)
    set(VECPHYS_NVCC_CONFIGS
        DEBUG -g -G
        RELEASE  -O2 --use_fast_math
        RELWITHDEBINFO -g -G  -O2 --use_fast_math
        MINSIZEREL --use_fast_math )
  endif()
endif()

# Enable pretty diagnostics coloring if gcc 4.9+ is being used
if (GNU)
  set(VECPHYS_ERROR_LIMIT 20 CACHE STRING
      "Limit number of errors output by diagnostics.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=${VECPHYS_ERROR_LIMIT}")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=auto")
  endif()
endif()

# Forward CMake options
if (VC_ACCELERATION)
  set(VECPHYS_COMPILATION_FLAGS
      "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_VC_ACCELERATION")
endif()
if (GUBENCHMARK)
  set(VECPHYS_COMPILATION_FLAGS
      "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_BENCHMARK")
endif()
if (ROOT)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_ROOT")
endif()

# Pass flags to compilers
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
#set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -O2")

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (CUDA)
  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

################################################################################

# Backend configuration

if (VC_ACCELERATION OR Vc)

  find_package(Vc REQUIRED)
  if (MIC)
    set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} ${Vc_MIC_LIBRARIES})
  else()
    set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} ${Vc_LIBRARIES})
  endif()
  include_directories(AFTER ${Vc_INCLUDE_DIR})

  if (Vc)
    set(SRC_EXTERNAL ${SRC_EXTERNAL} src/backend/vc/Backend.cpp)
  endif()

endif()

if (Cilk)

  if (NOT Intel)
    message(FATAL_ERROR "Must use Intel C++ compiler (icc) for Cilk backend.")
  endif()

  set(SRC_EXTERNAL ${SRC_EXTERNAL} src/backend/cilk/Backend.cpp)

endif()

if (CUDA)

  find_package(CUDA REQUIRED)
  include(FindCUDA)

  set(CUDA_SEPARABLE_COMPILATION ON)

  set(SRC_EXTERNAL ${SRC_EXTERNAL} src/CudaManager.cpp src/backend/cuda/Interface.cpp)
#  set(SRC_CUDA ${SRC_CUDA}
#      source/CudaManager.cu)

  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${VECPHYS_NVCC_COMPILATION_FLAGS})
endif()

if (ROOT)
  find_package(ROOT REQUIRED)

  include_directories(AFTER ${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL}
      ${ROOT_LIBRARIES} )
endif()

#----------------------------------------------------------------------------
# Find VecGeom geometry headers and library
#
################################################################################
find_package(VecGeom REQUIRED)
if( VECGEOM_FOUND )
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} ${VECGEOM_LIBRARIES} )
else()
   message(FATAL_ERROR "Required VecGeom library not found")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc
                    ${CMAKE_CURRENT_SOURCE_DIR}/inc/base
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/inc
                    ${VECGEOM_INCLUDE_DIR}
)

set(SRC_CPP_RELATIVE
  src/GUAliasSampler.cc
  src/ComptonKleinNishina.cc
  src/ConversionBetheHeitler.cc
  src/PhotoElectronSauterGavrila.cc
  src/IonisationMoller.cc
  src/BremSeltzerBerger.cc
  src/UrbanWentzelVI.cc
  src/GUAliasTable.cc
  src/GUAliasTableManager.cc
  src/Physics2DVector.cc
  src/EmModelBase.cc
  src/MaterialHandler.cc
  src/EmProcess.cc
  src/PhotonProcess.cc
  src/ElectronProcess.cc
)
foreach(SRC ${SRC_CPP_RELATIVE})
  set(SRC_CPP ${SRC_CPP} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC})
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

# Add files that won't be compiled for CUDA

if (GUBENCHMARK) 
  if (ROOT) 
    set(SRC_CPP ${SRC_CPP}
      src/GUHistogram.cc
    )
  endif()

  set(SRC_CPP ${SRC_CPP}
    src/GUTrackHandler.cc
  )
  set(SRC_CUDA ${SRC_CUDA}
    test/src/GUCurand.cu

  )

  # build benchmark executables separately
  if (GUBENCHMARK STREQUAL "Model") 
    set(SRC_CPP ${SRC_CPP}
      test/src/GUBenchmarker.cc
      test/src/GUBenchmarker_cpu.cc
    )
    set(SRC_CUDA ${SRC_CUDA}
      test/src/GUBenchmarker.cu
      test/src/GUBenchmarker_gpu.cu
    )
  elseif(GUBENCHMARK STREQUAL "Process")
    set(SRC_CPP ${SRC_CPP}
      test/src/ProcessBenchmarker.cc
      test/src/ProcessBenchmarker_cpu.cc
    )
    set(SRC_CUDA ${SRC_CUDA}
      test/src/ProcessBenchmarker.cu
      test/src/ProcessBenchmarker_gpu.cu
    )
  else()
    message(FATAL_ERROR "Invalid GUBENCHMARK=${GUBENCHMARK}: 
    Use GUBENCHMARK=[Model|Process]")
  endif()
endif()

# Compile cpp-files for C++11

set_source_files_properties(${SRC_CPP} PROPERTIES COMPILE_FLAGS -std=c++11)

################################################################################

# Build libraries

add_library(vecphys ${SRC_CPP})
target_link_libraries(vecphys ${VECPHYS_LIBRARIES}
                      ${VECPHYS_LIBRARIES_EXTERNAL})
set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphys)
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} rt)
endif()

set_source_files_properties(
  ${SRC_EXECUTABLES}
  PROPERTIES COMPILE_FLAGS -std=c++11
)

# build the CUDA version of the library
if (CUDA)
  cuda_add_library(
    vecphyscuda
    ${SRC_CUDA}
    OPTIONS ${CUDA_ARCH}
    ${VECPHYS_NVCC_CONFIGS}
  )
  target_link_libraries(vecphyscuda ${VECPHYS_LIBRARIES}
                        ${VECPHYS_LIBRARIES_EXTERNAL})
  set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphyscuda)
endif()

################################################################################

# Define executables

if (GUBENCHMARK STREQUAL "Model")
  set(TEST_EXECUTABLES_BENCHMARK
      ${CMAKE_CURRENT_SOURCE_DIR}/test/GUModelBenchmark.cc
  )
elseif (GUBENCHMARK STREQUAL "Process")
  set(TEST_EXECUTABLES_BENCHMARK
      ${CMAKE_CURRENT_SOURCE_DIR}/test/ProcessBenchmark.cc
  )
endif()

set(TEST_EXECUTABLES_CUDA
  # Files go here
)

set(TEST_EXECUTABLES
  ${TEST_EXECUTABLES_BENCHMARK}
  ${TEST_EXECUTABLES_CUDA}
)
set_source_files_properties(
  ${TEST_EXECUTABLES}
  PROPERTIES COMPILE_FLAGS -std=c++11
)

# Build executables

set(VECPHYS_EXECUTABLES)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} ${VECPHYS_LIBRARIES}
                          ${VECPHYS_LIBRARIES_EXTERNAL}
                          ${USERKERNELLIB})
    set(VECPHYS_EXECUTABLES ${VECPHYS_EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      add_test(${TARGET_NAME} ${TARGET_NAME})
    endif()
  endforeach()
endfunction()

build_executables("${TEST_EXECUTABLES_CORE}")
add_to_ctest("${TEST_EXECUTABLES_CORE}")
if (GUBENCHMARK)
  build_executables("${TEST_EXECUTABLES_BENCHMARK}")
endif()
if (CUDA)
  build_executables("${TEST_EXECUTABLES_CUDA}")
  add_to_ctest("${TEST_EXECUTABLES_CUDA}")
endif()

################################################################################

# Installation

export(PACKAGE VecPhys)

file(RELATIVE_PATH INSTALL_INCLUDE_DIR_RELATIVE
     "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}/.")
file(RELATIVE_PATH INSTALL_LIB_DIR_RELATIVE
     "${INSTALL_CMAKE_DIR}" "${INSTALL_LIB_DIR}/.")
# Build
set(CONF_TYPE "build")
set(CONF_INCLUDE_DIR "${CMAKE_SOURCE_DIR}")
set(CONF_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
set(CONF_LIBRARIES ${VECPHYS_LIBRARIES_EXTERNAL}
    "${CONF_LIBRARY_DIR}/libvecphys.a")
configure_file(VecPhysConfig.cmake.in
               "${PROJECT_BINARY_DIR}/VecPhysConfig.cmake" @ONLY)
# Installation
set(CONF_TYPE "install")
set(CONF_INCLUDE_DIR "\${THIS_DIR}/${INSTALL_INCLUDE_DIR_RELATIVE}")
set(CONF_LIBRARY_DIR "\${THIS_DIR}/${INSTALL_LIB_DIR_RELATIVE}")
set(CONF_LIBRARIES ${VECPHYS_LIBRARIES_EXTERNAL}
    "${CONF_LIBRARY_DIR}/libvecphys.a")
configure_file(VecPhysConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/VecPhysConfig.cmake" @ONLY)

# Install the VecPhysConfig.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/VecPhysConfig.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install headers and libraries
file(GLOB _headers LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/inc inc/*)
foreach(_file ${_headers})
  install(FILES inc/${_file} DESTINATION ${INSTALL_INCLUDE_DIR}/vecphys)
endforeach()

file(GLOB _headers_test LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/test/inc test/inc/*)
foreach(_file ${_headers_test})
  install(FILES test/inc/${_file} DESTINATION ${INSTALL_INCLUDE_DIR}/vecphys/test)
endforeach()
install(TARGETS vecphys DESTINATION ${INSTALL_LIB_DIR})
if (CUDA)
  install(TARGETS vecphyscuda DESTINATION ${INSTALL_LIB_DIR})
endif()

