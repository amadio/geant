#----------------------------------------------------------------------------------------------
# Add include directories, sources and headers
#
include_directories(
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBTestModules
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/models/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
  ${CMAKE_SOURCE_DIR}/core/base/inc
  ${CMAKE_SOURCE_DIR}/run/scheduler/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
)

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

set(sources
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSFieldConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSFullApp.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSData.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSParticleGun.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/src/CMSPhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/src/LHCbDetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/src/LHCbFullApp.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/src/LHCbData.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/src/LHCbParticleGun.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/src/LHCbPhysicsList.cxx   
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3App.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3DetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3Data.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3PrimaryGenerator.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3PhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/TestEm3FastSimProcess.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/src/HepMCTruth.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/TestEm5.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/TestEm5Data.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/TestEm5DetectorConstruction.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/TestEm5PrimaryGenerator.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/TestEm5PhysicsList.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/StepMaxProcess.cxx
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/src/Hist.cxx
  )

set(headers
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSFieldConstruction.h  
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSFullApp.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSData.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSParticleGun.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSPhysicsList.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/GeantV/inc/CMSmagneticField.h
  
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc/LHCbDetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc/LHCbFullApp.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc/LHCbData.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc/LHCbParticleGun.h
  ${CMAKE_SOURCE_DIR}/examples/physics/FullLHCb/GeantV/inc/LHCbPhysicsList.h

  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3App.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3DetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3Data.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3PrimaryGenerator.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3PhysicsList.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/TestEm3FastSimProcess.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm3/GeantV/inc/HepMCTruth.h

  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/TestEm5.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/TestEm5Data.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/TestEm5DetectorConstruction.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/TestEm5PrimaryGenerator.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/TestEm5PhysicsList.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/StepMaxProcess.h
  ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/Hist.h

  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Event.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Filter.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Module.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Producer.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Schedule.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/Path.cpp
  ${CMAKE_SOURCE_DIR}/examples/physics/cmsToyGV/TBBProcessingDemo/TBBFrameworkCore/EventProcessor.cpp
  )

#if (USE_ROOT)
#  ROOT_GENERATE_DICTIONARY(G__ExamplesRP ${headers} MODULE GeantExamplesRP
#    LINKDEF
#    ${CMAKE_SOURCE_DIR}/examples/physics/TestEm5/GeantV/inc/Geant/LinkDef.h
#    )
#    set(sources ${sources} G__ExamplesRP.cxx)
#endif()
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# Library and dependencies
#
add_library(GeantExamplesRP SHARED ${sources})

target_link_libraries(GeantExamplesRP -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics Material ${HEPMC_LIBRARIES}
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
# The GeantV version of the example applications
#

add_subdirectory(TestEm3/GeantV)
add_subdirectory(TestEm5/GeantV)
if (USE_ROOT)
  add_subdirectory(FullCMS/GeantV)
  add_subdirectory(FullLHCb/GeantV)
endif()

#----------------------------------------------------------------------------------------------
# The Geant4 version of the example applications (only if build with Geant4 was requested)
#
if (WITH_GEANT4)
  add_subdirectory(TestEm3/Geant4)
  add_subdirectory(TestEm5/Geant4)
  add_subdirectory(FullCMS/Geant4)
endif()

if(USE_TBB)
  add_subdirectory(cmsToyGV)
endif()
