#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if (USE_REAL_PHYSICS)
  include_directories(
    inc
    ../magneticfield/inc
    ../realphysics/dummies/inc
    ../realphysics/physics/management/inc
    ../examplesRP/GeantV/ExN03RP/inc
    ../examplesRP/GeantV/TestEm3RP/inc
    ../numa/inc
    )
else()
  include_directories(
    inc
    ../magneticfield/inc
    ../physics/particle/inc
    ../numa/inc
    )
endif()

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

set(sources
   ../base/src/Error.cxx
   src/GeantTaskData.cxx
   src/GeantTrack.cxx
   src/GeantTrackVec.cxx
   src/GeantTrackGeo.cxx
   src/GeantBasket.cxx
   src/Basket.cxx
   src/Handler.cxx
   src/PreStepHandler.cxx
   src/XSecSamplingHandler.cxx
   src/GeomQueryHandler.cxx
   src/FieldPropagationHandler.cxx
   src/LinearPropagationHandler.cxx
   src/ContinuousProcHandler.cxx
   src/DiscreteProcHandler.cxx
   src/SteppingActionsHandler.cxx
   src/SimulationStage.cxx
   src/PreStepStage.cxx
   src/XSecSamplingStage.cxx
   src/GeomQueryStage.cxx
   src/PropagationStage.cxx
   src/ContinuousProcStage.cxx
   src/DiscreteProcStage.cxx
   src/SteppingActionsStage.cxx
   src/sync_objects.cxx
   src/GeantOutput.cxx
   src/GeantPropagator.cxx
   src/GeantEvent.cxx
   src/GeantScheduler.cxx
   src/PhysicsInterface.cxx
   src/PhysicsProcessOld.cxx
   src/WorkloadManager.cxx
   src/GeantVApplication.cxx
   src/GeantFactoryStore.cxx
   src/GeantTrackStat.cxx
   src/MyHit.cxx
   src/StdApplication.cxx
   src/TThreadMergingFile.cxx
   src/TThreadMergingServer.cxx
   src/TransportManager.cxx
   src/LocalityManager.cxx
   src/MCTruthMgr.cxx
   src/GeantRunManager.cxx
   src/GeantEventServer.cxx
   src/TrackStat.cxx
)

set(headers
   ../base/inc/Geant/Config.h
   ../base/inc/Geant/Error.h
   inc/GeantTaskData.h
   inc/PhysicsInterface.h
   inc/PhysicsProcessOld.h
   inc/Basket.h
   inc/Handler.h
   inc/StackLikeBuffer.h
   inc/PreStepHandler.h
   inc/XSecSamplingHandler.h
   inc/GeomQueryHandler.h
   inc/FieldPropagationHandler.h
   inc/LinearPropagationHandler.h
   inc/ContinuousProcHandler.h
   inc/DiscreteProcHandler.h
   inc/SteppingActionsHandler.h
   inc/SimulationStage.h
   inc/PreStepStage.h
   inc/XSecSamplingStage.h
   inc/GeomQueryStage.h
   inc/PropagationStage.h
   inc/ContinuousProcStage.h
   inc/DiscreteProcStage.h
   inc/SteppingActionsStage.h
   inc/GeantBasket.h
   inc/GeantPropagator.h
   inc/GeantTrack.h
   inc/GeantTrackVec.h
   inc/GeantTrackGeo.h
   inc/GeantOutput.h
   inc/GeantVApplication.h
   inc/GeantTrackStat.h
   inc/MyHit.h
   inc/GeantFactory.h
   inc/StdApplication.h
   inc/PrimaryGenerator.h
   inc/TThreadMergingFile.h
   inc/TThreadMergingServer.h
   inc/TrackManager.h
   inc/TransportManager.h
   inc/LocalityManager.h
   inc/MCTruthMgr.h
   inc/cuckoohash_map.hh
   inc/GeantRunManager.h
   inc/GeantEventServer.h
   inc/GeantConfig.h
   inc/TrackStat.h
)

if(USE_ROOT AND (NOT USE_VECGEOM_NAVIGATOR))
   set(sources ${sources} src/ScalarNavInterfaceTGeo.cxx)
   set(headers ${headers} inc/ScalarNavInterfaceTGeo.h)
endif()

if(USE_VECGEOM_NAVIGATOR)
  set(sources ${sources}
    src/ScalarNavInterfaceVG.cxx
    src/ScalarNavInterfaceVGM.cxx
    src/VectorNavInterface.cxx
    )

  set(headers ${headers}
    inc/ScalarNavInterfaceVG.h
    inc/ScalarNavInterfaceVGM.h
    inc/VectorNavInterface.h
    )
endif()

if (USE_ROOT)
  if (USE_UMESIMD)
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF inc/LinkDef.h OPTIONS "-DVECGEOM_UMESIMD=1")
   else() 
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
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
  target_link_libraries(Geant_v -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield GeantNuma ${ROOT_LIBRARIES} )

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
  runAppRP.cc
  runAppTestEM3RP.cc
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
      elseif(TARGET_NAME STREQUAL "runAppTestEM3RP")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ATLASbar.root)
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

install(TARGETS runAppRP runAppTestEM3RP DESTINATION bin)
add_dependencies(runAppRP Geant_v RealPhysics GeantExamplesRP)
# -------  END == WITH REAL PHYSICS -------------------------------------------------------
endif()


#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
