#----------------------------------------------------------------------------
# Setup the project

include_directories(inc)
include_directories("${PROJECT_SOURCE_DIR}/magneticfield/inc")
## else: include_directories(inc ../magneticfield/inc)

if( USE_VECGEOM_NAVIGATOR )
# Bug in VC (?), cannot compile without optimisation
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_HOME_DIRECTORY}/base/src/Error.cxx
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
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TThreadMergingFile.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/TThreadMergingServer.cxx
)     

set(headers
   ${CMAKE_HOME_DIRECTORY}/base/inc/Geant/Config.h
   ${CMAKE_HOME_DIRECTORY}/base/inc/Geant/Error.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTaskData.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PhysicsProcess.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantBasket.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantPropagator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrack.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantOutput.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantVApplication.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantTrackStat.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyHit.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantFactory.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/StdApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/PrimaryGenerator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TThreadMergingFile.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/TThreadMergingServer.h
)   

ROOT_GENERATE_DICTIONARY(G__Geant ${headers} MODULE libGeant_v LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")

#----------------------------------------------------------------------------
#

set(GEANTV_HEADERS ${headers} )
set(GEANTV_SOURCES ${sources} )

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED

if (CUDA)
  # See vecprot_v2/cuda/CMakeLists.txt
  add_subdirectory (cuda)
else()
  add_library(Geant_v SHARED ${sources} G__Geant.cxx)

target_link_libraries(Geant_v Vmagfield ${ROOT_LIBRARIES} -lEG -lGeom )

if( VECGEOM_FOUND )
  target_link_libraries(Geant_v Vmagfield ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
  
   if (VECGEOM_GEANT4)
      find_package(Geant4 REQUIRED gdml)
      include(${Geant4_USE_FILE})
      message(${Geant4_INCLUDE_DIRS})
      include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIRS})
      target_link_libraries(Geant_v ${Geant4_LIBRARIES})
   endif()
endif()

add_custom_target(Geantlib_v DEPENDS Geant_v)
if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v.dylib DESTINATION lib)
else()
  install(TARGETS Geant_v DESTINATION lib)
endif()
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v.rootmap
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeant_v_rdict.pcm DESTINATION lib)

endif()

find_package(HepMC REQUIRED)
include_directories("${HEPMC_INCLUDE_DIR}")

#
#---------------------- Building executables // ctest ------------------------------------
#

set(TEST_EXECUTABLES
   ${CMAKE_CURRENT_SOURCE_DIR}/CMSApp.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/runApp.cc
)

#----------------------------------------------------------------------------

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    if(USE_VECGEOM_NAVIGATOR)
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v -lGeantExamples ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
    else()
      target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v -lXsec -lGeantExamples ${ROOT_LIBRARIES} -lGeom ${HEPMC_LIBRARIES})
    endif()
    add_dependencies(${TARGET_NAME} Geant_v Xsec GeantExamples)
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

#----------------------------------------------------------------------------

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "runApp")
        if(CUDA)
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ExN03.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT.root -r true)
        else()
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/ExN03.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT.root)
        endif()
      elseif(TARGET_NAME STREQUAL "CMSApp")
        if(CUDA)
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/cms2015.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT_G496p02_1mev.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT_G496p02_1mev.root -E ${CMAKE_SOURCE_DIR}/data/pp14TeVminbias.root -r true)
        else()
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}> -g ${CMAKE_SOURCE_DIR}/data/cms2015.root -x ${CMAKE_SOURCE_DIR}/data/xsec_FTFP_BERT_G496p02_1mev.root -f ${CMAKE_SOURCE_DIR}/data/fstate_FTFP_BERT_G496p02_1mev.root -E ${CMAKE_SOURCE_DIR}/data/pp14TeVminbias.root)
        endif()
      else()
        add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
    endif()
  endforeach()
endfunction()

build_executables("${TEST_EXECUTABLES}")

if(CTEST)
  add_to_ctest("${TEST_EXECUTABLES}")
endif()

install(TARGETS CMSApp runApp DESTINATION bin)
