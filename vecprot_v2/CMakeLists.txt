#----------------------------------------------------------------------------
# Setup the project

include_directories(inc ../magneticfield/inc)

if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_CURRENT_SOURCE_DIR}/src/GeantTaskData.cxx
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
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyHit.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/StdApplication.cxx
)     
     
set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTaskData.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PhysicsProcess.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantBasket.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantScheduler.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantPropagator.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrack.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantOutput.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/WorkloadManager.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantVApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrackStat.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyHit.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/StdApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PrimaryGenerator.h 
)   

ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE libGeant_v LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
#

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED

if (CUDA)
  CUDA_ADD_LIBRARY(Geant_v
                 ${sources} G__Geant.cxx
                 SHARED
                 OPTIONS ${CUDA_ARCH}
  )
else()
  add_library(Geant_v SHARED ${sources} G__Geant.cxx)
endif()

if( VECGEOM_FOUND )
  target_link_libraries(Geant_v ${ROOT_LIBRARIES} -lEG -lGeom ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
  # Bug in VC (?), cannot compile without optimisation
  #set(CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG} -O0")
else()
  target_link_libraries(Geant_v ${ROOT_LIBRARIES} -lEG -lGeom )
endif()

add_custom_target(Geantlib_v DEPENDS Geant_v)

find_package(HepMC REQUIRED)
include_directories("${HEPMC_INCLUDE_DIR}")

add_executable(CMSApp CMSApp.cc)

target_link_libraries(CMSApp -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v -lXsec -lGeantExamples
       ${ROOT_LIBRARIES} -lGeom ${HEPMC_LIBRARIES})

add_dependencies(CMSApp Geant_v Xsec GeantExamples)
