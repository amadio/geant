#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources )  # Common source  files
file(GLOB headers )  # Common include files
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${GEANT5_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${GEANT5_EXE_LINKER_FLAGS}")
set_source_files_properties(${sources} PROPERTY COMPILE_FLAGS                  "${GEANT5_COMPILE_FLAGS}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
include_directories(${PROJECT_SOURCE_DIR}/magneticfield/inc)

add_executable(testMagFieldEquation testMagFieldEquation.cxx) ## ${sources} ${headers})
target_link_libraries(testMagFieldEquation Vmagfield ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testMagFieldEquation")

add_executable(testTClassicalRK4 testTClassicalRK4.cxx) ## ${sources} ${headers})
target_link_libraries(testTClassicalRK4 Vmagfield ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testTClassicalRK4")

## add_executable(testFieldPropagator TestGUFieldPropagator.cxx) ## ${sources} ${headers})
## target_link_libraries(testFieldPropagator Vmagfield  ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
## # ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lGeom -lEG -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v)
## # add_dependencies(micromc xseclib Geantlib_v)
## message(STATUS "Created target for testFieldPropagator")

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
#add_custom_target(testFP DEPENDS testFieldPropagator)
#install(TARGETS testFieldPropagator DESTINATION bin)

add_custom_target(testMagEq DEPENDS testMagFieldEquation)
install(TARGETS testMagFieldEquation DESTINATION bin)
