#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(GeantV/ExN03RP/inc GeantV/TestEm3RP/inc ../vecprot_v2/inc ../numa/inc ../realphysics/dummies/inc ../realphysics/physics/management/inc)
if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

set(sources
  GeantV/ExN03RP/src/ExN03ApplicationRP.cxx
  GeantV/TestEm3RP/src/TestEM3ApplicationRP.cxx
  )

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
  GeantV/ExN03RP/src/ExN03ApplicationRP.cxx
  GeantV/TestEm3RP/src/TestEM3ApplicationRP.cxx
  )

set(headers
  GeantV/ExN03RP/inc/ExN03ApplicationRP.h
  GeantV/TestEm3RP/inc/TestEM3ApplicationRP.h
  )
if (USE_ROOT)
  ROOT_GENERATE_DICTIONARY(G__ExamplesRP ${headers} MODULE GeantExamplesRP
    LINKDEF   GeantV/ExN03RP/inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
    set(sources ${sources} G__ExamplesRP.cxx)
endif()

#----------------------------------------------------------------------------------------------
# Library and dependencies
#
add_library(GeantExamplesRP SHARED ${sources})

target_link_libraries(GeantExamplesRP -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics
  ${TBB_LIBRARIES})

if (USE_ROOT)
  target_link_libraries(GeantExamplesRP ${ROOT_LIBRARIES})
endif()

if (KNC)
  target_link_libraries(GeantExamplesRP ${Vc_MIC_LIBRARIES})
endif()

add_dependencies(GeantExamplesRP Geant_v RealPhysics)
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

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamplesRP_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------------------------
# The Geant4 version of the examples
#
add_subdirectory(Geant4/SimplyCalo)
