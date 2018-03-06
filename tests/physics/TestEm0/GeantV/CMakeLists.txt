#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/tests/physics/GeantV/TestEm0/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/tests/physics/GeantV/TestEm0/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/tests/physics/GeantV/TestEm0/inc/
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
set(OUTPUT bin/tests/physics/TestEm0/GeantV)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})
add_executable(TestEm0_GV ${CMAKE_SOURCE_DIR}/tests/physics/GeantV/TestEm0/TestEm0_GV.cc ${sources})
target_link_libraries(TestEm0_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics)

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#
set(TestEm0_SCRIPTS
    README
   )

foreach(_script ${TestEm0_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestEm0_GV DESTINATION ${OUTPUT})
install(FILES ${TestEm5_SCRIPTS} DESTINATION ${OUTPUT})
