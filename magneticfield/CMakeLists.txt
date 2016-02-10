#----------------------------------------------------------------------------
# Find VecPhys (required )
#----------------------------------------------------------------------------
find_package(VecGeom REQUIRED)
if(VECGEOM_FOUND)
  include_directories(${VECGEOM_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Required library VecGeom was not found!")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB sources src/*.cxx) # (T*.cxx src/GU*.cxx)
file(GLOB headers inc/*.h)   # (T*.h inc/GU*.h)

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
#include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

# print include directories
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include directories: ${inc_dirs}")

#----------------------------------------------------------------------------
# Generate the root dictionary: this library will be used in run.C root macro.  
#----------------------------------------------------------------------------
## ROOT_GENERATE_DICTIONARY(vMagFieldDict ${headers} MODULE libVphysproc LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(Vmagfield SHARED ${sources} ) ## vMagFieldDict.cxx)
target_link_libraries(Vmagfield  ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_headers("${headers}")
## add_dependencies

add_subdirectory(test)

install(TARGETS Vmagfield LIBRARY DESTINATION lib)
