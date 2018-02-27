#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/Geant/T*.h)

include_directories(inc)

ROOT_GENERATE_DICTIONARY(NudyDict ${headers} MODULE Nudy LINKDEF inc/Geant/LinkDef.h)

#----------------------------------------------------------------------------------------------
# Add the executable
#
add_library(Nudy SHARED ${sources} NudyDict.cxx)
target_link_libraries(Nudy -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${ROOT_LIBRARIES} MathMore ${VC_LIBRARIES})

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Nudy DESTINATION lib)
endif()

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy.rootmap 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

