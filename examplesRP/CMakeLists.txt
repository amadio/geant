#----------------------------------------------------------------------------------------------
# Add include directories, sources and headers
#
include_directories(
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/HAD/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/HAD/inc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/cross_sections/HAD/inc
  ${CMAKE_SOURCE_DIR}/material/inc
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/vecprot_v2/inc
  ${CMAKE_SOURCE_DIR}/numa/inc
)

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

set(sources
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/src/CMSDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/src/CMSFullApp.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/src/CMSData.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/src/CMSParticleGun.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/src/TestEm3App.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/src/TestEm3DetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/src/TestEm3Data.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/src/TestEm3PrimaryGenerator.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/src/TestEm3PhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/TestEm5.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/TestEm5Data.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserPrimaryGenerator.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/UserPhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/StepMaxProcess.cxx
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/src/Hist.cxx
  )

set(headers
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/inc/CMSDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/inc/CMSFullApp.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/inc/CMSData.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/FullCMS/inc/CMSParticleGun.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc/TestEm3App.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc/TestEm3DetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc/TestEm3Data.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc/TestEm3PrimaryGenerator.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm3/inc/TestEm3PhysicsList.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/TestEm5.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/TestEm5Data.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserPrimaryGenerator.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/UserPhysicsList.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/StepMaxProcess.h
  ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/Hist.h
  )

#if (USE_ROOT)
#  ROOT_GENERATE_DICTIONARY(G__ExamplesRP ${headers} MODULE GeantExamplesRP
#    LINKDEF
#    ${CMAKE_SOURCE_DIR}/examplesRP/GeantV/TestEm5/inc/LinkDef.h
#    OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
#    set(sources ${sources} G__ExamplesRP.cxx)
#endif()
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# Library and dependencies
#
add_library(GeantExamplesRP SHARED ${sources})

target_link_libraries(GeantExamplesRP -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics Material
  ${TBB_LIBRARIES})

if (USE_ROOT)
  target_link_libraries(GeantExamplesRP ${ROOT_LIBRARIES})
endif()

if (KNC)
  target_link_libraries(GeantExamplesRP ${Vc_MIC_LIBRARIES})
endif()

add_dependencies(GeantExamplesRP Geant_v RealPhysics Material)
add_custom_target(ExamplesRPlib DEPENDS GeantExamplesRP)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
  # for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP.so DESTINATION lib)
else()
  install(TARGETS GeantExamplesRP DESTINATION lib)
endif()

#if(USE_ROOT)
#  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP.rootmap
#    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP_rdict.pcm DESTINATION lib)
#endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------------------------
# The Geant4 version of the example applications (only if build with Geant4 was requested)
#
if (WITH_GEANT4)
  add_subdirectory(Geant4/TestEm3)
  add_subdirectory(Geant4/TestEm5)
  add_subdirectory(Geant4/FullCMS)
endif()
