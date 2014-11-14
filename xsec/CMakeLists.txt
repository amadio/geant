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


#----------------------------------------------------------------------------
# Setup include directory for this project
#


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/T*.cxx src/GunGenerator.cxx)
file(GLOB headers inc/T*.h inc/GunGenerator.h)

#set(headers ${headers} ${VC_INCLUDE_DIR}/base/Vector.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  set(sources "${sources} src/HepMCGenerator.cxx")
  set(headers "${headers} inc/HepMCGenerator.h")
endif()

set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${GEANT5_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${GEANT5_EXE_LINKER_FLAGS}")
set_source_files_properties(${sources} PROPERTY COMPILE_FLAGS                  "${GEANT5_COMPILE_FLAGS}")

include_directories(inc)

if(HepMC_FOUND)
  ROOT_GENERATE_DICTIONARY(xsecCint ${headers} LINKDEF inc/xsmcLinkDef.h)
else()
  ROOT_GENERATE_DICTIONARY(xsecCint ${headers} LINKDEF inc/xsecLinkDef.h)
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Xsec SHARED ${sources} xsecCint.cxx ${headers} xsecCint.h)

if(HepMC_FOUND)
  target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
    ${ROOT_LIBRARIES} -lEG -lGui -lGeom ${HEPMC_LIBRARIES})
else()
  target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
    ${ROOT_LIBRARIES} -lEG -lGui -lGeom)
endif()

add_dependencies(Xsec Geantlib_v)

add_headers("${headers}")

add_custom_target(xseclib DEPENDS Xsec)

