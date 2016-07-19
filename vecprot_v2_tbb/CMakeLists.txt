#cmake_minimum_required(VERSION 2.6.3 FATAL_ERROR)


#include_directories("${PROJECT_SOURCE_DIR}/vecprot_v2/inc")

# Find source files
file(GLOB SOURCES src/*.cxx)
# Find header files
file(GLOB HEADERS inc/*.h)

# Include header files
include_directories(inc)

# Include TBB
include_directories(${TBB_INCLUDE_DIR})

# Create shared library
add_library(Geant_tbb SHARED ${SOURCES})

target_link_libraries(Geant_tbb -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v
     ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS} ${TBB_LIBRARIES})

add_dependencies(Xsec Geantlib_v)

# Install library
install(TARGETS Geant_tbb DESTINATION lib)

# Install library headers
install(FILES ${HEADERS} DESTINATION inc/${PROJECT_NAME})
