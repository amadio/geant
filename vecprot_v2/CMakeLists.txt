#----------------------------------------------------------------------------
# Setup the project

include_directories(inc ../magneticfield/inc)

set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${GEANT5_COMPILE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG} -O1")
endif()
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANT5_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${GEANT5_EXE_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantThreadData.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrack.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantBasket.cxx  
   ${CMAKE_CURRENT_SOURCE_DIR}/src/sync_objects.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantOutput.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantPropagator.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantEvent.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantScheduler.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/PhysicsProcess.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/WorkloadManager.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantVApplication.cxx 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantFactoryStore.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTrackStat.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyHit.cxx
)     
     
set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantThreadData.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PhysicsProcess.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantBasket.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantScheduler.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantPropagator.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantEvent.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrack.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantOutput.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/WorkloadManager.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantVApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrackStat.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyHit.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PrimaryGenerator.h 
)   

#file(GLOB sources src/*.cxx)
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/GeantFactory.cxx" "" sources "${sources}")
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/GeantHitBlock.cxx" "" sources "${sources}")
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTaskManager.cxx" "" sources "${sources}")
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/GeantVRunnable.cxx" "" sources "${sources}")

#file(GLOB headers inc/*.h)
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantEMGamma.h" "" headers "${headers}")
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/ComptonCrossSection.h" "" headers "${headers}")
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h" "" headers "${headers}")

#ROOT_GENERATE_DICTIONARY(G__Geant ${headers} LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h)
ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE libGeant_v LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
#

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED

add_library(Geant_v SHARED ${sources} G__Geant.cxx)
if( VECGEOM_FOUND )
target_link_libraries(Geant_v ${ROOT_LIBRARIES} -lEG -lGeom ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
# Bug in VC (?), cannot compile without optimisation
#set(CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG} -O1")
else()
target_link_libraries(Geant_v ${ROOT_LIBRARIES} -lEG -lGeom )
endif()

add_custom_target(Geantlib_v DEPENDS Geant_v)
