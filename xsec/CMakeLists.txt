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
file(GLOB headers inc/T*.h)
include_directories(inc)
ROOT_GENERATE_DICTIONARY(xsecCint ${headers} LINKDEF inc/xsecLinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Xsec SHARED ${sources} xsecCint.cxx ${headers} xsecCint.h)
target_link_libraries(Xsec ${ROOT_LIBRARIES} -lEG)

add_headers("${headers}")


