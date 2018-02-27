#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
set(TEST_EXECUTABLES
  executable_test.C
  )

set(headers
   inc/GeantThreadData.h
   inc/PhysicsProcess.h
   inc/GeantBasket.h
   inc/GeantScheduler.h
   inc/Geant/Propagator.h
   inc/Geant/Event.h
   inc/Geant/Track.h
   inc/GeantOutput.h
   inc/Geant/WorkloadManager.h
   inc/Geant/UserApplication.h
   inc/TrackStat.h
   inc/Geant/MyHit.h
   ../examples/inc/ExN03Application.h
   ../inc/Geant/PrimaryGenerator.h
   ../xsec/inc/GunGenerator.h
   ../xsec/inc/TTabPhysProcess.h
   )

include_directories(inc ../magneticfield/inc)

#----------------------------------------------------------------------------------------------
# Make sure to use GEANTV and not GEANT4 flags, in case they have been overwritten
# Not sure this is really necessary
#
set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${GEANTV_COMPILE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_FLAGS}")
# Bug in VC (?), cannot compile without optimisation
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_LINKER_FLAGS}")

ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF inc/Geant/LinkDef.h)

#----------------------------------------------------------------------------------------------
# CTEST section
#
if(CTEST)
  message(STATUS "Testing with CTest enabled.")
  enable_testing()
endif()

#----------------------------------------------------------------------------------------------
# Add the test executable
#
set(GEANTV_EXECUTABLES)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    if(VECGEOM_FOUND)
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} Xsec Geant_v GeantExamples)
    else()
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${ROOT_LIBRARIES} Xsec Geant_v GeantExamples)
    endif()
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "executable_test")
      	add_test(NAME ${TARGET_NAME}
	  COMMAND  $<TARGET_FILE:${TARGET_NAME}> 4 true tests/data/ExN03.root tests/data/xsec_FTFP_BERT.root tests/data/fstate_FTFP_BERT.root true)
      else()
      	add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
      set_tests_properties(${TARGET_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED")
    endif()
  endforeach()
endfunction()

build_executables("${TEST_EXECUTABLES}")
add_to_ctest("${TEST_EXECUTABLES}")
install(TARGETS executable_test DESTINATION tests)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
