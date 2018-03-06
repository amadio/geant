#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/tests/physics/xsectionsTest/GeantV/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/tests/physics/xsectionsTest/GeantV/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests/physics/xsectionsTest/GeantV/inc/
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
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/inc
  )

#----------------------------------------------------------------------------------------------
# Executable
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/physics/xsectionsTest/GeantV)
add_executable(xsectionsTest ${CMAKE_SOURCE_DIR}/tests/physics/xsectionsTest/GeantV/xsectionsTest.cc ${sources})
target_link_libraries(xsectionsTest -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics Geant_v ${VECGEOM_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS xsectionsTest DESTINATION bin/tests/physics/xsectionsTest)
