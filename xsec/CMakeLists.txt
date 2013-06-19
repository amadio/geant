#----------------------------------------------------------------------------
# Setup the project

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
#find_package(ROOT REQUIRED)
#include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB lsources ${PROJECT_SOURCE_DIR}/xsec/src/*.cxx)
file(GLOB lheaders ${PROJECT_SOURCE_DIR}/xsec/inc/T*.h)
include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)
ROOT_GENERATE_DICTIONARY(xsecCint ${lheaders} LINKDEF ${PROJECT_SOURCE_DIR}/xsec/inc/xsecLinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Xsec SHARED ${lsources} xsecCint.cxx ${lheaders} xsecCint.h)
target_link_libraries(Xsec ${ROOT_LIBRARIES} -lEg)



