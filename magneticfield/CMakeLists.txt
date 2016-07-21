#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cxx) # (T*.cxx src/GU*.cxx)
file(GLOB headers inc/*.h)   # (T*.h inc/GU*.h)

include_directories(inc)

#----------------------------------------------------------------------------------------------
# Add the library
#
# ROOT_GENERATE_DICTIONARY(vMagFieldDict ${headers} MODULE Vphysproc LINKDEF inc/LinkDef.h)
add_library(Vmagfield SHARED ${sources} ) ## vMagFieldDict.cxx)
target_link_libraries(Vmagfield  -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
add_headers("${headers}")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS Vmagfield LIBRARY DESTINATION lib)

add_subdirectory(test)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
