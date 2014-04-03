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

set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_DEBUG          "${GEANT5_CXX_FLAGS_DEBUG}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_MINISIZEREL    "${GEANT5_CXX_FLAGS_MINISIZEREL}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set_source_files_properties(${sources} PROPERTY CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set_source_files_properties(${sources} PROPERTY CMAKE_EXEC_LINKER_FLAGS        "${GEANT5_EXE_LINKER_FLAGS}")
set_source_files_properties(${sources} PROPERTY COMPILE_FLAGS                  "${GEANT5_COMPILE_FLAGS}")

include_directories(inc)
ROOT_GENERATE_DICTIONARY(xsecCint ${headers} LINKDEF inc/xsecLinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Xsec SHARED ${sources} xsecCint.cxx ${headers} xsecCint.h)
target_link_libraries(Xsec ${ROOT_LIBRARIES} -lEG -lGui -lGeom)

add_headers("${headers}")

add_custom_target(xseclib DEPENDS Xsec)

