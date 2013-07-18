#----------------------------------------------------------------------------
# Setup the project

include_directories(inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx)
string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/GeantEMGamma.cxx" "" sources "${sources}")
string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/src/ComptonCrossSection.cxx" "" sources "${sources}")

file(GLOB headers inc/*.h)
string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantEMGamma.h" "" headers "${headers}")
string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/ComptonCrossSection.h" "" headers "${headers}")
string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h" "" headers "${headers}")

#ROOT_GENERATE_DICTIONARY(geantCint ${headers} LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h)
ROOT_GENERATE_DICTIONARY(geantCint ${headers} LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Geant SHARED ${sources} geantCint.cxx ${headers} geantCint.h)
target_link_libraries(Geant ${ROOT_LIBRARIES} -lEG -lGeom)

add_custom_target(Geantlib DEPENDS Geant)
