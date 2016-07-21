#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources )  # Common source  files
file(GLOB headers )  # Common include files

include_directories(${CMAKE_SOURCE_DIR}/magneticfield/inc)

#----------------------------------------------------------------------------------------------
# Redefine compilation flags since Geant4 is overriding them
#
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${GEANT5_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${GEANT5_EXE_LINKER_FLAGS}")
set_source_files_properties(${sources} PROPERTY COMPILE_FLAGS                  "${GEANT5_COMPILE_FLAGS}")

#----------------------------------------------------------------------------------------------
# Add the executables
#

add_executable(testMagFieldEquation testMagFieldEquation.cxx) ## ${sources} ${headers})
target_link_libraries(testMagFieldEquation -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testMagFieldEquation")

# add_executable(testTClassicalRK4 testTClassicalRK4.cxx) ## ${sources} ${headers})
# target_link_libraries(testTClassicalRK4 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
# message(STATUS "Created target for testTClassicalRK4")

add_executable(testStepperFixed testStepperFixed.cxx) ## ${sources} ${headers})
target_link_libraries(testStepperFixed -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testStepperFixed")

add_executable(testIntegrationDriver testIntegrationDriver.cxx) ## ${sources} ${headers})
target_link_libraries(testIntegrationDriver -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testIntegrationDriver")

add_executable(testGUPropagateField testGUPropagateField.cxx) ## ${sources} ${headers})
target_link_libraries(testGUPropagateField -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testGUPropagateField")

## add_executable(testFieldPropagator TestGUFieldPropagator.cxx) ## ${sources} ${headers})
## target_link_libraries(testFieldPropagator Vmagfield  -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
## # ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Xsec Geant_v)
## # add_dependencies(micromc xseclib Geantlib_v)
## message(STATUS "Created target for testFieldPropagator")

#add_custom_target(testFP DEPENDS testFieldPropagator)
add_custom_target(testMagEq DEPENDS testMagFieldEquation)

#----------------------------------------------------------------------------------------------
# Install the files
#
#install(TARGETS testFieldPropagator DESTINATION bin)
install(TARGETS testMagFieldEquation DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

