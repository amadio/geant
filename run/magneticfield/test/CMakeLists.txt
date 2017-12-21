#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources )  # Common source  files
file(GLOB headers )  # Common include files

include_directories(${CMAKE_SOURCE_DIR}/magneticfield/inc ${CMAKE_SOURCE_DIR}/VecMagFieldRoutine/include)
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
# include_directories(${PROJECT_SOURCE_DIR}/magneticfield/inc)
# include_directories(${PROJECT_SOURCE_DIR}/magneticfield/inc ${PROJECT_SOURCE_DIR}/VecMagFieldRoutine/include)


## Need -lGeant_v in order to obtain the function/template Geant::cxx::ErrorHandlerImpl

#add_executable(testMagFieldEquation testMagFieldEquation.cxx) ## ${sources} ${headers})

#target_link_libraries(testMagFieldEquation Vmagfield Geant_v ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
# target_link_libraries(testMagFieldEquation Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testMagFieldEquation")

#add_executable(testVectorMagFieldEquation testVectorMagFieldEquation.cxx)
#target_link_libraries(testVectorMagFieldEquation Vmagfield Geant_v ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testVectorMagFieldEquation")

## Test for 'simple' scalar class
#add_executable(testScalarCashKarp testScalarCashKarp.cxx)
#target_link_libraries(testScalarCashKarp Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testScalarCashKarp")

## Test for 'vector' Stepper class  ( renamed from testNewStepperFixedCashKarp )
#add_executable(testVectorCashKarp testVectorCashKarp.cxx)
#target_link_libraries(testVectorCashKarp Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testVectorCashKarp")

## New test for 'vector' Stepper class  ( mid-version towards flexible test below )
#add_executable(templatedTestVectorCashKarp templatedTestVectorCashKarp.cxx)
#target_link_libraries(templatedTestVectorCashKarp Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for templatedTestVectorCashKarp")

## Test newest Cash Karp, that can integrate scalar or Double_v
add_executable(testFlexibleCashKarp testFlexibleCashKarp.cxx)
target_link_libraries(testFlexibleCashKarp Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testFlexibleCashKarp")
## --------  End of Stepper Tests

## ########  Start of Integration Driver Tests #############
## 
add_executable(testScalarIntegrationDriver testScalarIntegrationDriver.cxx)
target_link_libraries(testScalarIntegrationDriver Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testScalarIntegrationDriver")

add_executable(testVectorIntegrationDriver testVectorIntegrationDriver.cxx)
target_link_libraries(testVectorIntegrationDriver Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testVectorIntegrationDriver")
## --------  End   of Integration Driver Tests -------------

# add_executable(testToyClass testToyClass.cxx)
# target_link_libraries(testToyClass Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
# message(STATUS "Created target for testToyClass")

# add_executable(BenchmarkIntegrationDriver BenchmarkIntegrationDriver.cxx) ## ${sources} ${headers})
# target_link_libraries(BenchmarkIntegrationDriver Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
# message(STATUS "Created target for BenchmarkIntegrationDriver")

#add_executable(TestConstBzStepper TestConstBzStepper.cxx) ## ${sources} ${headers})
#target_link_libraries(TestConstBzStepper Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

#add_executable(ReproduceError ReproduceError.cxx)
#target_link_libraries(ReproduceError Vmagfield Geant_v ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#message(STATUS "Created target for testVectorMagFieldEquation")

# add_executable(testTClassicalRK4 testTClassicalRK4.cxx) ## ${sources} ${headers})
# target_link_libraries(testTClassicalRK4 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield ${VECGEOM_LIBRARIES})
# message(STATUS "Created target for testTClassicalRK4")

#add_executable(testStepperFixed testStepperFixed.cxx)
#target_link_libraries(testStepperFixed -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield Geant_v ${VECGEOM_LIBRARIES})
#message(STATUS "Created target for testStepperFixed")

### Old test(s) ------
# add_executable(testIntegrationDriver testIntegrationDriver.cxx) ## ${sources} ${headers})
# target_link_libraries(testIntegrationDriver -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield Geant_v ${VECGEOM_LIBRARIES})
# message(STATUS "Created target for testIntegrationDriver")

#add_executable(testGUPropagateField testGUPropagateField.cxx) ## ${sources} ${headers})
#target_link_libraries(testGUPropagateField -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield Geant_v ${VECGEOM_LIBRARIES})
#message(STATUS "Created target for testGUPropagateField")

add_executable(testConstVecFieldStepper testConstVecFieldStepper.cxx) ## ${sources} ${headers})
target_link_libraries(testConstVecFieldStepper Vmagfield Geant_v ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Created target for testConstVecFieldStepper")

## add_executable(testFieldPropagator TestGUFieldPropagator.cxx) ## ${sources} ${headers})
## target_link_libraries(testFieldPropagator -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield Geant_v ${VECGEOM_LIBRARIES})
## # ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Xsec Geant_v)
## # add_dependencies(micromc xseclib Geantlib_v)
## message(STATUS "Created target for testFieldPropagator")

#add_custom_target(testFP DEPENDS testFieldPropagator)
add_custom_target(testMagEq DEPENDS testMagFieldEquation)

#----------------------------------------------------------------------------------------------
# Install the files
#
#install(TARGETS testFieldPropagator DESTINATION bin)
#install(TARGETS testMagFieldEquation DESTINATION bin)
#install(TARGETS testVectorMagFieldEquation DESTINATION bin)
# Stepper tests
#install(TARGETS testScalarCashKarp DESTINATION bin)
#install(TARGETS testVectorCashKarp DESTINATION bin)
#install(TARGETS templatedTestVectorCashKarp DESTINATION bin)
install(TARGETS testFlexibleCashKarp DESTINATION bin)
# Driver  tests
## install(TARGETS testFlexibleIntegrationDriver DESTINATION bin)
install(TARGETS testVectorIntegrationDriver DESTINATION bin)
# Old tests - steppers
#install(TARGETS testStepperFixed DESTINATION bin)
install(TARGETS testConstVecFieldStepper DESTINATION bin)
# Old tests - drivers
install(TARGETS testScalarIntegrationDriver DESTINATION bin)
#install(TARGETS testGUPropagateField DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

