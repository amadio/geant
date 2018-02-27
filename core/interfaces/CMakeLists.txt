#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories("${PROJECT_SOURCE_DIR}/core/interfaces/inc")

file(GLOB sources src/*.cxx)
file(GLOB headers inc/Geant/*.h)

# library building currently handled by run/scheduler/CMakeLists.txt

install(FILES ${headers} DESTINATION inc/Geant)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
