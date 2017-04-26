#----------------------------------------------------------------------------------------------
# Configuration options
#
set(GUBENCHMARK "" CACHE STRING "Enable performance comparisons.")
option(DISTANCE_DEBUG "Enable comparision of calculcated distances againt ROOT/Geant4 behind the scenes" OFF)
option(VECPHYS_STANDALONE "Build VECPHYS standalone" OFF)
set(VECPHYS_VECTOR "" CACHE STRING "Vector instruction set to use (sse/avx/...).")

#----------------------------------------------------------------------------------------------
# Make relative paths absolute
#
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

#----------------------------------------------------------------------------------------------
# Determine backend
#
set(Vc)
set(Cilk)
set(Scalar)

if (NOT BACKEND)
  message(STATUS "Backend not specified. Defaulting to Vc.")
  set(BACKEND "Vc")
endif()

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
  message(FATAL_ERROR "vecphys: Invalid backend specified: \"${BACKEND}\"")
endif()
message(STATUS "vecphys: Configuring with vector backend ${BACKEND}.")

#----------------------------------------------------------------------------------------------
# Determine C++ compiler being used
#
if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
    OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
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

#----------------------------------------------------------------------------------------------
# Set C++ compiler flags. Will also be forwarded by CUDA when compiling C++.
#
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

#----------------------------------------------------------------------------------------------
# Compiler specific optimization flags
#
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

#----------------------------------------------------------------------------------------------
# Backend specific
#
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

#----------------------------------------------------------------------------------------------
# CUDA specific
#
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

#----------------------------------------------------------------------------------------------
# Enable pretty diagnostics coloring if gcc 4.9+ is being used
#
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

#----------------------------------------------------------------------------------------------
# Forward CMake options
#
if (VC_ACCELERATION)
  set(VECPHYS_COMPILATION_FLAGS
      "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_VC_ACCELERATION")
endif()
if (GUBENCHMARK)
  set(VECPHYS_COMPILATION_FLAGS
      "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_BENCHMARK")
endif()
if (USE_ROOT)
  set(VECPHYS_COMPILATION_FLAGS "${VECPHYS_COMPILATION_FLAGS} -DVECPHYS_ROOT")
endif()

#----------------------------------------------------------------------------------------------
# Pass flags to compilers
#
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

message(STATUS "vecphys: Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (CUDA)
  message(STATUS "vecphys: Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

#----------------------------------------------------------------------------------------------
# Backend configuration
#
if (VC_ACCELERATION OR Vc)
  if (KNC)
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
    message(FATAL_ERROR "vecphys: Must use Intel C++ compiler (icc) for Cilk backend.")
  endif()
  set(SRC_EXTERNAL ${SRC_EXTERNAL} src/backend/cilk/Backend.cpp)
endif()

if (CUDA)
  set(SRC_EXTERNAL ${SRC_EXTERNAL} src/CudaManager.cpp src/backend/cuda/Interface.cpp)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${VECPHYS_NVCC_COMPILATION_FLAGS})
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(inc
  inc/base
  test/inc
  ${CMAKE_SOURCE_DIR}/physics/particle/inc
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

set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} ${VECGEOM_LIBRARIES} )

if (CUDA)
#----------------------------------------------------------------------------------------------
# Copy all source files to .cu-files in order for NVCC to compile them as CUDA
# code and not regular C++ files.
#
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

#----------------------------------------------------------------------------------------------
# Add files that won't be compiled for CUDA
#
if (GUBENCHMARK) 
  if (USE_ROOT) 
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
    message(FATAL_ERROR "vecphys: Invalid GUBENCHMARK=${GUBENCHMARK}: 
    Use GUBENCHMARK=[Model|Process]")
  endif()
endif()

#----------------------------------------------------------------------------------------------
# Build libraries
#
add_library(vecphys ${SRC_CPP})
target_link_libraries(vecphys -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECPHYS_LIBRARIES}
  ${VECPHYS_LIBRARIES_EXTERNAL})
set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphys)

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(VECPHYS_LIBRARIES_EXTERNAL ${VECPHYS_LIBRARIES_EXTERNAL} rt)
endif()

set_source_files_properties(
  ${SRC_EXECUTABLES}
  PROPERTIES COMPILE_FLAGS -std=c++11
  )

#----------------------------------------------------------------------------------------------
# Build the CUDA version of the library
#
if (CUDA)
  cuda_add_library(
    vecphyscuda
    ${SRC_CUDA}
    SHARED
    OPTIONS ${CUDA_ARCH}
    ${VECPHYS_NVCC_CONFIGS}
    )

  target_link_libraries(vecphyscuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECPHYS_LIBRARIES}
    ${VECPHYS_LIBRARIES_EXTERNAL})
  set(VECPHYS_LIBRARIES ${VECPHYS_LIBRARIES} vecphyscuda)
endif()

#----------------------------------------------------------------------------------------------
# Define executables
#
if (GUBENCHMARK STREQUAL "Model")
  set(TEST_EXECUTABLES_BENCHMARK
    test/GUModelBenchmark.cc
    )
elseif (GUBENCHMARK STREQUAL "Process")
  set(TEST_EXECUTABLES_BENCHMARK
    test/ProcessBenchmark.cc
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

#----------------------------------------------------------------------------------------------
# Build executables
#
set(VECPHYS_EXECUTABLES)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECPHYS_LIBRARIES}
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

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
