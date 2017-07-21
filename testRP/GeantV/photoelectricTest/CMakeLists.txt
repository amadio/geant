#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources ${CMAKE_SOURCE_DIR}/testRP/GeantV/photoelectricTest/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/testRP/GeantV/photoelectricTest/inc/*.h)

include_directories(
  ${CMAKE_SOURCE_DIR}/testRP/GeantV/photoelectricTest/inc/
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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/testRP/photoelectricTest/GeantV)
add_executable(photoelectricTest_GV ${CMAKE_SOURCE_DIR}/testRP/GeantV/photoelectricTest/photoelectricTest_GV.cc ${sources})
target_link_libraries(photoelectricTest_GV ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics Geant_v)
