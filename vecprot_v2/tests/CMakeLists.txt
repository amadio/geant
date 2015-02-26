#----------------------------------------------------------------------------
# Setup the project

include_directories(inc ../magneticfield/inc)

set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${GEANT5_COMPILE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG} -O0")
endif()
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANT5_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${GEANT5_EXE_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(TEST_EXECUTABLES
   ${CMAKE_CURRENT_SOURCE_DIR}/executable_test.C
)     
     
set(headers
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantThreadData.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/PhysicsProcess.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantBasket.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantScheduler.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantPropagator.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantEvent.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantTrack.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantOutput.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/WorkloadManager.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantVApplication.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/GeantTrackStat.h
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/MyApplication.h 
   ${CMAKE_SOURCE_DIR}/vecprot_v2/inc/MyHit.h 
   ${CMAKE_SOURCE_DIR}/inc/PrimaryGenerator.h
   ${CMAKE_SOURCE_DIR}/xsec/inc/GunGenerator.h
   ${CMAKE_SOURCE_DIR}/xsec/inc/TTabPhysProcess.h
)   

ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE libGeant_v LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
if(CTEST)
  message(STATUS "Testing with CTest enabled.")
  enable_testing()
endif()

set(GEANTV_EXECUTABLES)
link_directories(${PROJECT_SOURCE_DIR}/lib)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    if(VECGEOM_FOUND)
      target_link_libraries(${TARGET_NAME} ${ROOT_LIBRARIES} -lEG -lGeom ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} -lXsec -lGeant_v)
    else()
      target_link_libraries(${TARGET_NAME} ${ROOT_LIBRARIES} -lEG -lGeom -lXsec -lGeant_v)
    endif()
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      add_test(NAME ${TARGET_NAME} COMMAND $<TARGET_FILE:${TARGET_NAME}> 4 false ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/ExN03.root ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/xsec_FTFP_BERT.root ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/fstate_FTFP_BERT.root)
      set_tests_properties(${TARGET_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED")
    endif()
  endforeach()
endfunction()

#if(TEST_EXECUTABLES STREQUAL "executable_test")
#  set(OPTIONS "4 false ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/ExN03.root ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/xsec_FTFP_BERT.root ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/data/fstate_FTFP_BERT.root")
#endif()

build_executables("${TEST_EXECUTABLES}")
add_to_ctest("${TEST_EXECUTABLES}")
install(TARGETS executable_test DESTINATION ${CMAKE_SOURCE_DIR}/vecprot_v2/tests)
