# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

PROJECT(GPTracking)

# cuda
SET(CMAKE_MODULE_PATH 
  ${GPTracking_SOURCE_DIR}/CMakeModule 
  ${CMAKE_MODULE_PATH}
  )

# disable warning on the FindCUDA module
CMAKE_POLICY(SET CMP0017 OLD)

SET(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")

FIND_PACKAGE(CUDA REQUIRED)

SET(CUDA_NVCC_FLAGS "-arch=sm_20 --ptxas-options=-v")

# Use CUDA separate compilation and linking
SET(CUDA_SEPARATE_COMPILATION ON)

# include/library dirs
INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/include 
                     Common/include 
                     EMPhysics/include 
                     Geometry/include 
                     Manager/include 
                     Material/include 
                     Transportation/include 
                     kernel/include
)

LINK_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/lib64 )

#------------------------------------------------------------------------------
# - input data
#------------------------------------------------------------------------------

file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})

FILE( GLOB src_files */src/*.cu )

# executables
cuda_add_executable ( tracking
                      test/tracking.cc
                      ${src_files}
)

#target_link_libraries(test ${CUDA_LIBRARIES} )
