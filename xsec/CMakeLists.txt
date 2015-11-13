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
file(GLOB sources src/T*.cxx src/GunGenerator.cxx src/TestProcess.cxx)
file(GLOB headers inc/T*.h inc/GunGenerator.h inc/TestProcess.h)

#set(XSEC_HEADERS ${headers} PARENT_SCOPE)
#set(XSEC_SOURCES ${sources} PARENT_SCOPE)

#set(headers ${headers} ${VC_INCLUDE_DIR}/base/Vector.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)


if (NOT MIC)
  ROOT_GENERATE_DICTIONARY(xsecDict ${headers} MODULE libXsec LINKDEF inc/xsecLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  #----------------------------------------------------------------------------
  # Add the executable, and link it to the Geant4 libraries
  #
  add_library(Xsec SHARED ${sources} xsecDict.cxx)

  set(ROOT_LIBS "")
  if(USE_ROOT OR NOT USE_VECGEOM_NAVIGATOR)
    set(ROOT_LIBS "${ROOT_LIBRARIES} -lEG -lGui -lGeom")
  endif()
else()
  add_library(Xsec SHARED ${sources})
endif()

set(HEPMC_LIBS "")
if(HepMC_FOUND)
  set(HEPMC_LIBS "${HEPMC_LIBRARIES}")
endif()

if (NOT MIC)
target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})
else()
target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v 
     ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})
endif()

set(XSEC_LINK_LIBS "-lXsec" ${HEPMC_LIBRARIES} PARENT_SCOPE)

add_dependencies(Xsec Geantlib_v)

add_headers("${headers}")

add_custom_target(xseclib DEPENDS Xsec)
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Xsec DESTINATION lib)
endif()
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec_rdict.pcm DESTINATION lib)

