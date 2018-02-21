#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(
  inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/run/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/physics/generators/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm5/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm3/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullCMS/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullLHCb/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  )

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

# Bug in VC (?), cannot compile without optimisation
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")

set(sources
   ${CMAKE_SOURCE_DIR}/core/base/src/Error.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTaskData.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/GeantTrack.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrackGeo.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Basket.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Handler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SimulationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldLookup.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LinearPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PropagationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsStage.cxx
   ${CMAKE_SOURCE_DIR}/core/concurrency/src/sync_objects.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantPropagator.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/GeantEvent.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/EventSet.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/PhysicsInterface.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/WorkloadManager.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/GeantVApplication.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/GeantVDetectorConstruction.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/GeantFactoryStore.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/MyHit.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/StdApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TransportManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LocalityManager.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/MCTruthMgr.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantRunManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantEventServer.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackStat.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Region.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVG.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVGM.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/VectorNavInterface.cxx
)

set(headers
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Config.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Error.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/globals.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/GeantFwd.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/dcqueue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/mpmc_bounded_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/priority_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/condition_locker.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTaskData.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/PhysicsInterface.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Basket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Handler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Basketizer.h
 #  ${CMAKE_CURRENT_SOURCE_DIR}/inc/StackLikeBuffer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SimulationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PreStepHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeomQueryHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/FieldPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/FieldLookup.h   
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinearPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SteppingActionsHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PreStepStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeomQueryStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PropagationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SteppingActionsStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantPropagator.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantTrack.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrackGeo.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/GeantVApplication.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantVDetectorConstruction.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/MyHit.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantFactory.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantFactoryStore.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/StdApplication.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/PrimaryGenerator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackManager.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TransportManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/LocalityManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/MCTruthMgr.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/cuckoohash_map.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantRunManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantEvent.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/EventSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantEventServer.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantConfig.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackStat.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Region.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/ExternalFramework.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/ScalarNavInterfaceVG.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/ScalarNavInterfaceVGM.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/VectorNavInterface.h
)

if (USE_ROOT)
  if (USE_UMESIMD)
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h OPTIONS "-DVECGEOM_UMESIMD=1")
   else()
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
   endif()
endif()

#----------------------------------------------------------------------------
#

set(GEANTV_HEADERS ${headers} )
set(GEANTV_SOURCES ${sources} )

install(FILES ${GEANTV_HEADERS} DESTINATION inc)

#----------------------------------------------------------------------------------------------
# CUDA section
#
if (CUDA)
  # See vecprot_v2/cuda/CMakeLists.txt
  add_subdirectory (cuda)
else()
  if(USE_ROOT)
    set(sources ${sources} G__Geant.cxx)
  endif()
  add_library(Geant_v SHARED ${sources})
  target_link_libraries(Geant_v -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield GeantNuma ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} Material)

  if (VECGEOM_GEANT4)
    include(${Geant4_USE_FILE})
    include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIRS})
    target_link_libraries(Geant_v ${Geant4_LIBRARIES})
  endif()

#----------------------------------------------------------------------------------------------
# Install files
#
  add_custom_target(Geantlib_v DEPENDS Geant_v)
  add_dependencies(Geantlib_v GeantNuma)

  if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
    install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
  else()
    install(TARGETS Geant_v DESTINATION lib)
  endif()
  if (USE_ROOT)
    install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v.rootmap
      ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v_rdict.pcm DESTINATION lib)
  endif()

endif()

#----------------------------------------------------------------------------------------------
# Build executables
#

set(TEST_EXECUTABLES
  examples/TestEm3_GV.cc
  examples/TestEm5_GV.cc
  )
set(NOTEST_EXECUTABLES
  examples/fullCMS.cc
  examples/fullLHCb.cc
  )

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics GeantExamplesRP ${ROOT_LIBRARIES} -lpthread ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
#    if (USE_TBB)
#      target_link_libraries(${TARGET_NAME} Geant_tbb  ${TBB_LIBS})
#    endif()
    add_dependencies(${TARGET_NAME} Geant_v RealPhysics GeantExamplesRP GeantNuma)
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "TestEm3_GV")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      elseif(TARGET_NAME STREQUAL "TestEm5_GV")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      else()
        add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
    endif()
  endforeach()
endfunction()

link_directories(${PROJECT_BINARY_DIR}/lib)

build_executables("${TEST_EXECUTABLES}")
build_executables("${NOTEST_EXECUTABLES}")
if(CTEST)
  add_to_ctest("${TEST_EXECUTABLES}")
endif()


#----------------------------------------------------------------------------
# Copy all scripts and install the executable(s)
#
### TestEm5_GV ###
set(TestEm5_GV_OUTPUT bin/examples/TestEm5/GeantV)
set(TestEm5_GV_SCRIPTS examples/TestEm5_GV_hanson19.mac examples/hanson_19.dat examples/README_TestEm5_GV)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TestEm5_GV_OUTPUT})
foreach(_script ${TestEm5_GV_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${TestEm5_GV_OUTPUT}/${_script}
    COPYONLY
    )
endforeach()
install(TARGETS TestEm5_GV DESTINATION ${TestEm5_GV_OUTPUT})
install(FILES ${TestEm5_GV_SCRIPTS} DESTINATION ${TestEm5_GV_OUTPUT})
### fullCMS ###
set(fullCMS_OUTPUT bin/examples/FullCMS/GeantV)
set(fullCMS_SCRIPTS examples/fullCMS.mac examples/README_fullCMS ${CMAKE_SOURCE_DIR}/examples/physics/Geant4/FullCMS/cms.gdml)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${fullCMS_OUTPUT})
foreach(_script ${fullCMS_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${fullCMS_OUTPUT}/${_script}
    COPYONLY
    )
endforeach()
install(TARGETS fullCMS DESTINATION ${fullCMS_OUTPUT})
install(FILES ${fullCMS_SCRIPTS} DESTINATION ${fullCMS_OUTPUT})
### fullLHCb ###
set(fullLHCb_OUTPUT bin/examples/FullLHCb/GeantV)
set(fullLHCb_SCRIPTS examples/fullLHCb.mac examples/README_fullLHCb ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullLHCb/LHCb_201603.gdml)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${fullLHCb_OUTPUT})
foreach(_script ${fullLHCb_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${fullLHCb_OUTPUT}/${_script}
    COPYONLY
    )
endforeach()
install(TARGETS fullLHCb DESTINATION ${fullLHCb_OUTPUT})
install(FILES ${fullLHCb_SCRIPTS} DESTINATION ${fullLHCb_OUTPUT})

### TestEm3_GV ###
set(TestEm3_GV_OUTPUT bin/examples/Calorimeter/GeantV)
set(TestEm3_GV_SCRIPTS examples/TestEm3_GV_ATLASbar.mac examples/README_TestEm3_GV)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TestEm3_GV_OUTPUT})
foreach(_script ${TestEm3_GV_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${TestEm3_GV_OUTPUT}/${_script}
    COPYONLY
    )
endforeach()
install(TARGETS TestEm3_GV DESTINATION ${TestEm3_GV_OUTPUT})
install(FILES ${TestEm3_GV_SCRIPTS} DESTINATION ${TestEm3_GV_OUTPUT})

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
