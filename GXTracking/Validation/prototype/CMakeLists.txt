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
project(electronTest)

#------------------------------------------------------------------------------
# - Find the package named 'Threads'
#------------------------------------------------------------------------------
SET(GXTracking_DIR "${CMAKE_SOURCE_DIR}/..")
SET(GPTracking_DIR "${GXTracking_DIR}/GPTracking")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${GXTracking_DIR}/cmake/Modules)
FIND_PACKAGE(CUDA REQUIRED)
#find_package(ROOT REQUIRED)

#------------------------------------------------------------------------------
# - cuda SDK and libs
#------------------------------------------------------------------------------
SET(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
SET(CUDA_NVCC_FLAGS -arch=sm_20 --optimize 2 --use_fast_math --ptxas-options=-v)

#message(STATUS "From prototype/CMakeLists.txt: GPTracking_DIR=" ${GPTracking_DIR})
#message(STATUS "From prototype/CMakeLists.txt: ROOT_INCLUDE_DIR=" ${ROOT_INCLUDE_DIR})

#------------------------------------------------------------------------------
# - Incude GPTracking
#------------------------------------------------------------------------------
INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/include 
                     ${GXTracking_DIR}/SimCore/Common/include
                     ${GXTracking_DIR}/gpu/cuda/kernel/include 
                     ${GXTracking_DIR}/gpu/cuda/kernel/src
                     ${GPTracking_DIR}/Common/include 
                     ${GPTracking_DIR}/EMPhysics/include 
                     ${GPTracking_DIR}/Geometry/include 
                     ${GPTracking_DIR}/Manager/include 
                     ${GPTracking_DIR}/Material/include 
                     ${GPTracking_DIR}/Transportation/include 
                     ${GPTracking_DIR}/Common/src 
                     ${GPTracking_DIR}/EMPhysics/src
                     ${GPTracking_DIR}/Geometry/src 
                     ${GPTracking_DIR}/Manager/src
                     ${GPTracking_DIR}/Material/src 
                     ${GPTracking_DIR}/Transportation/src 
		     ${GXTracking_DIR}/Validation/utils/include
                     ${ROOT_INCLUDE_DIR}
)

LINK_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR}/lib64 )

#------------------------------------------------------------------------------
# - compile and library 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
# - executable and link 
#------------------------------------------------------------------------------

#message(STATUS "from prototype CMakeLists: PROJECT_NAME=" ${PROJECT_NAME})
cuda_add_executable ( ${PROJECT_NAME}
                      ${PROJECT_NAME}.cc
                      src/${PROJECT_NAME}_kernel.cu
                      include/${PROJECT_NAME}_kernel.h
		      )

target_link_libraries(${PROJECT_NAME} -L${CMAKE_BINARY_DIR}/utils -lutils ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
