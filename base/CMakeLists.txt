# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()


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

add_headers("${headers}")
