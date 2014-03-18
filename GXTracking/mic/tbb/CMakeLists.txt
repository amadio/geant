#------------------------------------------------------------------------------
# - CMakeLists.txt for transportation
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# - Project definition
#------------------------------------------------------------------------------
project(${Project})

SET(PROJECT_TOP_DIR "./../../")
#------------------------------------------------------------------------------
# - Intel(R) Manycore Platform Software Stack 
#------------------------------------------------------------------------------
SET(MIC_TOOLKIT_ROOT_DIR "/usr/linux-k1om-4.7")
SET(CMAKE_TOOLCHAIN_FILE "../../cmake/Templates/cmake_mic.cmake")
add_definitions(-DUSE_MIC)
#------------------------------------------------------------------------------
# - Incude GPTracking
#------------------------------------------------------------------------------
INCLUDE_DIRECTORIES( ${MIC_TOOLKIT_ROOT_DIR}/include 
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

LINK_DIRECTORIES( ${MIC_TOOLKIT_ROOT_DIR}/lib64 )

#------------------------------------------------------------------------------
# - executable and link 
#------------------------------------------------------------------------------
add_executable ( ${Project}
                 test/${Project}.cc
                 kernel/src/${Project}_mic.cc
)

target_link_libraries(${Project})

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
