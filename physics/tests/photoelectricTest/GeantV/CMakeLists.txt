#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/physics/tests/photoelectricTest/GeantV/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/physics/tests/photoelectricTest/GeantV/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/physics/tests/photoelectricTest/GeantV/inc/
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
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
  )

#----------------------------------------------------------------------------------------------
# Executable
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/physics/tests/photoelectricTest/GeantV)
add_executable(photoelectricTest_GV ${CMAKE_SOURCE_DIR}/physics/tests/photoelectricTest/GeantV/photoelectricTest_GV.cc ${sources})
target_link_libraries(photoelectricTest_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics Geant_v)
