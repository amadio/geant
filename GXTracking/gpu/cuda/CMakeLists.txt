#------------------------------------------------------------------------------
# - CMakeLists.txt for transportation
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# - Project definition
#------------------------------------------------------------------------------
project(${Project})

SET(PROJECT_TOP_DIR "./../../")
#------------------------------------------------------------------------------
# - Find the package named 'Threads'
#------------------------------------------------------------------------------
FIND_PACKAGE(CUDA REQUIRED)

if(${USE_ROOT})
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_TOP_DIR}/cmake/Modules)
  find_package(ROOT REQUIRED)
endif()

#------------------------------------------------------------------------------
# - cuda SDK and libs 
#------------------------------------------------------------------------------
SET(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
#------------------------------------------------------------------------------
# - Build Options: Release (default --optimize 3) or Debug (-G -lineinfo)
#------------------------------------------------------------------------------
SET(CMAKE_BUILD_TYPE Release)
SET(CUDA_NVCC_FLAGS -arch=sm_35 --use_fast_math --ptxas-options=-v)
#------------------------------------------------------------------------------
# - Incude GPTracking
#------------------------------------------------------------------------------
INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/include 
                     ${ROOT_INCLUDE_DIR}
                     ${PROJECT_TOP_DIR}/GPTracking/Common/include 
                     ${PROJECT_TOP_DIR}/GPTracking/EMPhysics/include 
                     ${PROJECT_TOP_DIR}/GPTracking/Geometry/include 
                     ${PROJECT_TOP_DIR}/GPTracking/Manager/include 
                     ${PROJECT_TOP_DIR}/GPTracking/Material/include 
                     ${PROJECT_TOP_DIR}/GPTracking/Transportation/include 
                     ${PROJECT_TOP_DIR}/GPTracking/Common/src 
                     ${PROJECT_TOP_DIR}/GPTracking/EMPhysics/src
                     ${PROJECT_TOP_DIR}/GPTracking/Geometry/src 
                     ${PROJECT_TOP_DIR}/GPTracking/Manager/src 
                     ${PROJECT_TOP_DIR}/GPTracking/Material/src 
                     ${PROJECT_TOP_DIR}/GPTracking/Transportation/src 
                     ${PROJECT_TOP_DIR}/SimCore/Common/include 
                     ${PROJECT_TOP_DIR}/SimCore/EMPhysics/include 
                     ${PROJECT_TOP_DIR}/SimCore/Geometry/include 
                     ${PROJECT_TOP_DIR}/SimCore/Manager/include
                     ${PROJECT_TOP_DIR}/SimCore/Transportation/include
                     ${PROJECT_TOP_DIR}/SimCore/Common/src
                     ${PROJECT_TOP_DIR}/SimCore/EMPhysics/src 
                     ${PROJECT_TOP_DIR}/SimCore/Geometry/src
                     ${PROJECT_TOP_DIR}/SimCore/Manager/src
                     ${PROJECT_TOP_DIR}/SimCore/Transportation/src 
                     ${PROJECT_TOP_DIR}/FWCore/Application/include 
                     ${PROJECT_TOP_DIR}/FWCore/FrameWork/include 
                     ${PROJECT_TOP_DIR}/FWCore/Utility/include 
                     ${PROJECT_TOP_DIR}/FWCore/Application/src 
                     ${PROJECT_TOP_DIR}/FWCore/FrameWork/src 
                     kernel/include
                     kernel/src
)

LINK_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/lib64 )

#------------------------------------------------------------------------------
# - compile and library 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
# - executable and link 
#------------------------------------------------------------------------------
cuda_add_executable ( ${Project}
                      test/${Project}.cc
                      kernel/src/${Project}_kernel.cu
)

target_link_libraries(${Project} ${CUDA_LIBRARIES} ${ROOT_LIBRARIES} )

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
