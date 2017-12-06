#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/Geant/*.h)
add_headers("${headers}")

include_directories(inc 
                    ${CMAKE_SOURCE_DIR}/numa/inc)

install(FILES ${headers} DESTINATION inc/Geant)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
