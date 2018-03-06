cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(elasticTest_GV)

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/tests/physics/elasticTest/GeantV/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/tests/physics/elasticTest/GeantV/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests/physics/elasticTest/GeantV/inc/
  ${VECGEOM_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/core/numa/inc/
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc/
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc/
  ${CMAKE_SOURCE_DIR}/run/scheduler/inc/
  ${CMAKE_SOURCE_DIR}/core/base/inc/
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/models/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
)

#----------------------------------------------------------------------------------------------
# Executable
#
# add the clear model-level version of the test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/physics/elasticTest/GeantV)
add_executable(elasticTest_GV ${CMAKE_SOURCE_DIR}/tests/physics/elasticTest/GeantV/elasticTest_GV.cc ${sources})
target_link_libraries(elasticTest_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics Geant_v)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS elasticTest_GV DESTINATION bin/tests/physics/elasticTest_GV)
