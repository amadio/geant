# Adding to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB sources src/*.cxx) # (T*.cxx src/GU*.cxx)
file(GLOB headers inc/*.h)   # (T*.h inc/GU*.h)

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(inc)

# print include directories
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include directories: ${inc_dirs}")

#----------------------------------------------------------------------------
# Generate the root dictionary: this library will be used in run.C root macro.  
#----------------------------------------------------------------------------
## ROOT_GENERATE_DICTIONARY(vMagFieldDict ${headers} MODULE Vphysproc LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(Vmagfield SHARED ${sources} ) ## vMagFieldDict.cxx)
target_link_libraries(Vmagfield  ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_headers("${headers}")
## add_dependencies

add_subdirectory(test)

install(TARGETS Vmagfield LIBRARY DESTINATION lib)

