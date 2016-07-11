# Add this directpory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx src/GunGenerator.cxx)
file(GLOB headers inc/*.h inc/GunGenerator.h)
string(REGEX REPLACE "[^ ]*/fastsimLinkDef.h" "" headers "${headers}")

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)
if(USE_ROOT)
ROOT_GENERATE_DICTIONARY(fastsimDict ${headers} MODULE Fastsim LINKDEF inc/fastsimLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
endif()
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
if(USE_ROOT)
add_library(Fastsim SHARED ${sources} fastsimDict.cxx)

target_link_libraries(Fastsim -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
else()
add_library(Fastsim SHARED ${sources})

target_link_libraries(Fastsim -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
endif()
set(FASTSIM_LINK_LIBS "-lFastsim" ${HEPMC_LIBRARIES} PARENT_SCOPE)

add_dependencies(Fastsim Geantlib_v)

add_headers("${headers}")

add_custom_target(fastsimlib DEPENDS Fastsim)
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Fastsim DESTINATION lib)
endif()
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim_rdict.pcm DESTINATION lib)

