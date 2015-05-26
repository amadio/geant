#----------------------------------------------------------------------------
# Setup the project

#----------------------------------------------------------------------------
# Setup include directory for this project
#

include_directories(inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/Geant/[A-Z]*.h)

# ROOT_GENERATE_DICTIONARY(baseDict ${headers} LINKDEF src/baseLinkDef.h)

#----------------------------------------------------------------------------
# Add library
#

#add_library(Xsec SHARED ${sources} xsecCint.cxx ${headers} xsecCint.h)
#target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
#                           ${ROOT_LIBRARIES} -lEG -lGui -lGeom)
#add_dependencies(Xsec Geantlib_v)

add_headers("${headers}")

#add_custom_target(xseclib DEPENDS Xsec)
