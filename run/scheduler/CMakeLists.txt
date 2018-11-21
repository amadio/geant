#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(
  inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/physics/generators/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  )

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

# Bug in VC (?), cannot compile without optimisation
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")

set(sources
   ${CMAKE_SOURCE_DIR}/core/base/src/Error.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TaskData.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Track.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackGeo.cxx
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
   ${CMAKE_CURRENT_SOURCE_DIR}/src/Propagator.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Event.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/EventSet.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/PhysicsInterface.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/WorkloadManager.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/UserApplication.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/UserDetectorConstruction.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/FactoryStore.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/MyHit.cxx
   ${CMAKE_SOURCE_DIR}/run/userapp/src/StdApplication.cxx
#  ${CMAKE_CURRENT_SOURCE_DIR}/src/TransportManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LocalityManager.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/MCTruthMgr.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/RunManager.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/EventServer.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackStat.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/Region.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVG.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/ScalarNavInterfaceVGM.cxx
   ${CMAKE_SOURCE_DIR}/core/interfaces/src/VectorNavInterface.cxx
)

set(headers
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Config.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Error.h
   ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/Fwd.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/Geant/dcqueue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/mpmc_bounded_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/priority_queue.h
   ${CMAKE_SOURCE_DIR}/core/concurrency//inc/Geant/condition_locker.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/TaskData.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/RngWrapper.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/PhysicsInterface.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Basket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Handler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Basketizer.h
 #  ${CMAKE_CURRENT_SOURCE_DIR}/inc/StackLikeBuffer.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SimulationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PreStepHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/GeomQueryHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/FieldPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/FieldLookup.h   
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LinearPropagationHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SteppingActionsHandler.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PreStepStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/GeomQueryStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/PropagationStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/SteppingActionsStage.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/Propagator.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Track.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/TrackGeo.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/UserApplication.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/UserDetectorConstruction.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/MyHit.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Factory.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/FactoryStore.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/StdApplication.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/PrimaryGenerator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/TrackManager.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TransportManager.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LocalityManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/MCTruthMgr.h
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc/Geant/cuckoohash_map.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/RunManager.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Event.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/EventSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/EventServer.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/GeantConfig.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TrackStat.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/Region.h
   ${CMAKE_SOURCE_DIR}/run/userapp/inc/Geant/ExternalFramework.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/ScalarNavInterfaceVG.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/ScalarNavInterfaceVGM.h
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc/Geant/VectorNavInterface.h
)

if (USE_ROOT)
  if (USE_UMESIMD)
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LinkDef.h OPTIONS "-DVECGEOM_UMESIMD=1")
   else()
     ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE Geant_v LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/LinkDef.h )
   endif()
endif()

#----------------------------------------------------------------------------
#

set(GEANTV_HEADERS ${headers} )
set(GEANTV_SOURCES ${sources} )

install(FILES ${GEANTV_HEADERS} DESTINATION inc/Geant)

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
  target_link_libraries(Geant_v -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Vmagfield GeantNuma ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES}
          Material ${VecCoreLib_LIBRARIES})

  if (VECGEOM_GEANT4)
    #include(${Geant4_USE_FILE})
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
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
