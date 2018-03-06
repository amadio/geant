#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/tests/physics/TestHad0/GeantV/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/tests/physics/TestHad0/GeantV/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests/physics/TestHad0/GeantV/inc/
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
  ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
  )

#----------------------------------------------------------------------------------------------
# Executable
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/physics/TestHad0/GeantV)
add_executable(TestHad0_GV ${CMAKE_SOURCE_DIR}/tests/physics/TestHad0/GeantV/TestHad0_GV.cc ${sources})
target_link_libraries(TestHad0_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics)
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestHad0_GV DESTINATION bin/tests/physics/TestHad0_GV)
