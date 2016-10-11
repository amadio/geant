#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(inc ../include ${TBB_INCLUDE_DIR})

file(GLOB SOURCES src/*.cxx)
# Find header files
file(GLOB HEADERS inc/*.h)

# Create shared library
add_library(Geant_tbb SHARED ${SOURCES})

target_link_libraries(Geant_tbb -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v ${TBB_LIBRARIES})

# Install library
install(TARGETS Geant_tbb DESTINATION lib)

# Install library headers
install(FILES ${HEADERS} DESTINATION inc/${PROJECT_NAME})
