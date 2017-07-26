#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/testRP/GeantV/bremTest/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/testRP/GeantV/bremTest/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/testRP/GeantV/bremTest/inc/
  ${VECGEOM_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/numa/inc/
  ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/
  ${CMAKE_SOURCE_DIR}/base/inc/
  ${CMAKE_SOURCE_DIR}/material/inc
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/utils/inc
  )

#----------------------------------------------------------------------------------------------
# Executable
#
set(OUTPUT bin/testRP/bremTest/GeantV)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${OUTPUT})
add_executable(bremTest_GV ${CMAKE_SOURCE_DIR}/testRP/GeantV/bremTest/bremTest_GV.cc ${sources})
target_link_libraries(bremTest_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics Geant_v)

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#
set(bremTest_SCRIPTS
    README
   )

foreach(_script ${bremTest_SCRIPTS})
  configure_file(
    ${_script}
     ${CMAKE_INSTALL_PREFIX}/${OUTPUT}/${_script}
    COPYONLY
    )
    configure_file(
      ${_script}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
      COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS bremTest_GV DESTINATION ${OUTPUT})
