#----------------------------------------------------------------------------
# Setup the project


#----------------------------------------------------------------------------
# Find HepMC (required package)
#
find_package(HepMC QUIET)
if(HepMC_FOUND)
  message(STATUS "HepMC includes: ${HEPMC_INCLUDE_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
endif()

#set(HepMC_FOUND NO)

#----------------------------------------------------------------------------
# Setup include directory for this project
#


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx src/GunGenerator.cxx)
file(GLOB headers inc/*.h inc/GunGenerator.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)

ROOT_GENERATE_DICTIONARY(fastsimDict ${headers} MODULE libFastsim LINKDEF inc/fastsimLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Fastsim SHARED ${sources} fastsimDict.cxx)

set(ROOT_LIBS "")
if(USE_ROOT OR NOT USE_VECGEOM_NAVIGATOR)
  set(ROOT_LIBS "${ROOT_LIBRARIES} -lEG -lGui -lGeom")
endif()

set(HEPMC_LIBS "")
if(HepMC_FOUND)
  set(HEPMC_LIBS "${HEPMC_LIBRARIES}")
endif()

target_link_libraries(Fastsim -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})

set(FASTSIM_LINK_LIBS "-lFastsim" ${HEPMC_LIBRARIES} PARENT_SCOPE)

add_dependencies(Fastsim Geantlib_v)

add_headers("${headers}")

add_custom_target(fastsimlib DEPENDS Fastsim)
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim.dylib DESTINATION lib)
else()
  install(TARGETS Fastsim DESTINATION lib)
endif()
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFastsim_rdict.pcm DESTINATION lib)

