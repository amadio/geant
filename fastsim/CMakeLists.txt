#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc
  ${CMAKE_SOURCE_DIR}/physics/particle/inc
  ${CMAKE_SOURCE_DIR}/numa/inc)

file(GLOB sources src/*.cxx src/GunGenerator.cxx)
file(GLOB headers inc/*.h inc/GunGenerator.h)
string(REGEX REPLACE "[^ ]*/fastsimLinkDef.h" "" headers "${headers}")

#----------------------------------------------------------------------------------------------
# Add the library
#
if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(fastsimDict ${headers} MODULE Fastsim LINKDEF inc/fastsimLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  add_library(Fastsim SHARED ${sources} fastsimDict.cxx) 
  target_link_libraries(Fastsim -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v NewPhysics
    ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})

else()
  add_library(Fastsim SHARED ${sources})
  target_link_libraries(Fastsim -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v NewPhysics
    ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
endif()

set(FASTSIM_LINK_LIBS "Fastsim" ${HEPMC_LIBRARIES} PARENT_SCOPE)
add_dependencies(Fastsim Geantlib_v NewPhysics)
add_headers("${headers}")
add_custom_target(fastsimlib DEPENDS Fastsim)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Fastsim DESTINATION lib)
endif()
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim.rootmap 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
