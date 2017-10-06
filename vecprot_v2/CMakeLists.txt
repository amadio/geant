#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if (USE_REAL_PHYSICS)
  include_directories(
    inc
    ${CMAKE_SOURCE_DIR}/numa/inc
    ${CMAKE_SOURCE_DIR}/magneticfield/inc
    ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
    ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc
    ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc
    ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc
    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/ExN03RP/inc
    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3RP/inc
    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc
    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/CaloRP/inc
    )
else()
  include_directories(
    inc
    ${CMAKE_SOURCE_DIR}/numa/inc
    ${CMAKE_SOURCE_DIR}/magneticfield/inc
    ${CMAKE_SOURCE_DIR}/physics/particle/inc
    ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
    )
endif()

if (USE_VECGEOM_NAVIGATOR)
  include_directories(${CMAKE_SOURCE_DIR}/material/inc)
endif()

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

set(sources
   ${CMAKE_SOURCE_DIR}/base/src/Error.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTaskData.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrack.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrackVec.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrackGeo.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantBasket.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Basket.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Handler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SimulationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LinearPropagationHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsHandler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PreStepStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeomQueryStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PropagationStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/SteppingActionsStage.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/sync_objects.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantOutput.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantPropagator.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantEvent.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantScheduler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PhysicsInterface.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PhysicsProcessOld.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/WorkloadManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantVApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantVDetectorConstruction.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantFactoryStore.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyHit.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/StdApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TThreadMergingFile.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TThreadMergingServer.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TransportManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LocalityManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MCTruthMgr.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantRunManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantEventServer.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackStat.cxx
)

set(headers
   ${CMAKE_SOURCE_DIR}/base/inc/Geant/Config.h
   ${CMAKE_SOURCE_DIR}/base/inc/Geant/Error.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTaskData.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PhysicsInterface.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PhysicsProcessOld.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Basket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Handler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/StackLikeBuffer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SimulationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PreStepHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeomQueryHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/FieldPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinearPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SteppingActionsHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PreStepStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeomQueryStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PropagationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/SteppingActionsStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantBasket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantPropagator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrack.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrackVec.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrackGeo.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantOutput.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantVApplication.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantVDetectorConstruction.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyHit.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantFactory.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/StdApplication.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PrimaryGenerator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TThreadMergingFile.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TThreadMergingServer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TransportManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/LocalityManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MCTruthMgr.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/cuckoohash_map.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantRunManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantEventServer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantConfig.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackStat.h
)

# Add codes that describes some stages (and handlers) that are different when real-physics is used
if (NOT USE_REAL_PHYSICS)
  set (sources ${sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/XSecSamplingStage.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/XSecSamplingHandler.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ContinuousProcStage.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ContinuousProcHandler.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DiscreteProcStage.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DiscreteProcHandler.cxx
    )
  set (headers ${headers}
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/XSecSamplingStage.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/XSecSamplingHandler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ContinuousProcStage.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ContinuousProcHandler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/DiscreteProcStage.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/DiscreteProcHandler.h
    )
endif ()

if(USE_ROOT AND (NOT USE_VECGEOM_NAVIGATOR))
   set(sources ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/src/ScalarNavInterfaceTGeo.cxx)
   set(headers ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/inc/ScalarNavInterfaceTGeo.h)
endif()

if(USE_VECGEOM_NAVIGATOR)
  set(sources ${sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ScalarNavInterfaceVG.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ScalarNavInterfaceVGM.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src/VectorNavInterface.cxx
    )

  set(headers ${headers}
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ScalarNavInterfaceVG.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ScalarNavInterfaceVGM.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/VectorNavInterface.h
    )
endif()

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
  if (NOT USE_VECGEOM_NAVIGATOR)
    target_link_libraries(Geant_v -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield GeantNuma ${ROOT_LIBRARIES} )
  else()
    target_link_libraries(Geant_v -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield GeantNuma ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} Material)
  endif()
  if( VECGEOM_FOUND )
    target_link_libraries(Geant_v ${VC_LIBRARIES})

    if (VECGEOM_GEANT4)
      include(${Geant4_USE_FILE})
      include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIRS})
      target_link_libraries(Geant_v ${Geant4_LIBRARIES})
    endif()
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
# -------  START == NOT WITH REAL PHYSICS -------------------------------------------------------
if(NOT (USE_REAL_PHYSICS))
# -------
if(CUDA)
  set(TEST_EXECUTABLES
    CMSGPUApp.cc
    runGPUApp.cc
#   testBasketizer.cc
)
else()
  set(TEST_EXECUTABLES
    CMSApp.cc
    runApp.cc
    testBasketizer.cc
)
endif()

if(USE_ROOT)
  set(TEST_EXECUTABLES
    ${TEST_EXECUTABLES}
    testLocality.cc
    )
endif()

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    if(USE_VECGEOM_NAVIGATOR)
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Xsec Geant_v GeantExamples ${ROOT_LIBRARIES} -lpthread ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
    else()
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v Xsec GeantExamples ${ROOT_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
    endif()
    if (USE_TBB)
      target_link_libraries(${TARGET_NAME} Geant_tbb  ${TBB_LIBS})
    endif()
    add_dependencies(${TARGET_NAME} Geant_v Xsec GeantExamples GeantNuma)
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "runApp")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ExN03.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT.root)
      elseif(TARGET_NAME STREQUAL "runGPUApp")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ExN03.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT.root)
      elseif(TARGET_NAME STREQUAL "CMSGPUApp")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/cms2015.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT_G496p02_1mev.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT_G496p02_1mev.root -E ${CMAKE_SOURCE_DIR}/data/pp14TeVminbias.root)
      elseif(TARGET_NAME STREQUAL "CMSApp")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/cms2015.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT_G496p02_1mev.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT_G496p02_1mev.root -E ${CMAKE_SOURCE_DIR}/data/pp14TeVminbias.root)
      elseif(TARGET_NAME STREQUAL "testBasketizer")
        add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> 8)
      elseif(TARGET_NAME STREQUAL "testLocality")
        add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> 8 ${CMAKE_SOURCE_DIR}/data/ExN03.root)
      else()
        add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
    endif()
  endforeach()
endfunction()

link_directories(${PROJECT_BINARY_DIR}/lib)

build_executables("${TEST_EXECUTABLES}")

if(CTEST)
  add_to_ctest("${TEST_EXECUTABLES}")
endif()

if(CUDA)
  install(TARGETS CMSGPUApp runGPUApp DESTINATION bin)
else()
  install(TARGETS CMSApp runApp testBasketizer DESTINATION bin)
endif()

if(USE_ROOT)
  install(TARGETS testLocality DESTINATION bin)
endif()

if(CUDA)
  add_dependencies(runGPUApp Geant_v Xsec GeantExamples)
  add_dependencies(CMSGPUApp Geant_v Xsec GeantExamples)
else()
  add_dependencies(runApp Geant_v Xsec GeantExamples)
  add_dependencies(CMSApp Geant_v Xsec GeantExamples)
endif()
# -------  END == NOT WITH REAL PHYSICS -------------------------------------------------------


else()
# -------  START == WITH REAL PHYSICS -------------------------------------------------------
set(TEST_EXECUTABLES
  caloAppRP.cc
  runAppRP.cc
  runAppTestEM3RP.cc
  runTestEm5.cc
  )

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics GeantExamplesRP ${ROOT_LIBRARIES} -lpthread ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
    if (USE_TBB)
      target_link_libraries(${TARGET_NAME} Geant_tbb  ${TBB_LIBS})
    endif()
    add_dependencies(${TARGET_NAME} Geant_v RealPhysics GeantExamplesRP GeantNuma)
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "runAppRP")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ExN03.root)
      elseif(TARGET_NAME STREQUAL "caloAppRP")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      elseif(TARGET_NAME STREQUAL "runAppTestEM3RP")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ATLASbar.root)
      elseif(TARGET_NAME STREQUAL "runTestEm5")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      else()
        add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
    endif()
  endforeach()
endfunction()

link_directories(${PROJECT_BINARY_DIR}/lib)

build_executables("${TEST_EXECUTABLES}")

if(CTEST)
  add_to_ctest("${TEST_EXECUTABLES}")
endif()


#----------------------------------------------------------------------------
# Copy all scripts and install the executable(s)
#
### TestEm5 ###
set(TestEm5_OUTPUT bin/examples/TestEm5/GeantV)
set(TestEm5_SCRIPTS hanson19_GV.mac hanson_19.dat)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TestEm5_OUTPUT})
foreach(_script ${TestEm5_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_INSTALL_PREFIX}/${TestEm5_OUTPUT}/${_script}
    COPYONLY
    )
endforeach()
install(TARGETS runTestEm5 DESTINATION ${TestEm5_OUTPUT})
### The remaining ### (calo applications; runAppRP and runAppTestEM3RP will be removed soon)
set(caloApp_OUTPUT bin/examples/Calorimeter/GeantV)
set(caloApp_SCRIPTS caloApp.mac)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${caloApp_OUTPUT})
  configure_file(
    caloApp.mac
    ${CMAKE_INSTALL_PREFIX}/${caloApp_OUTPUT}/caloApp.mac
    COPYONLY
    )
install(TARGETS caloAppRP runAppRP runAppTestEM3RP DESTINATION ${caloApp_OUTPUT})
# -------  END == WITH REAL PHYSICS -------------------------------------------------------
endif()


#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
