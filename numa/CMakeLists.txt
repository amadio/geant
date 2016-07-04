#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
set(sources 
${CMAKE_CURRENT_SOURCE_DIR}/src/NumaNode.cxx
${CMAKE_CURRENT_SOURCE_DIR}/src/NumaPolicy.cxx
${CMAKE_CURRENT_SOURCE_DIR}/src/NumaTopology.cxx
${CMAKE_CURRENT_SOURCE_DIR}/src/NumaUtils.cxx
)
set(headers
${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantNuma.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaBlock.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaBlockMgr.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaNode.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaPolicy.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaTopology.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaUtils.h
)

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(GeantNuma SHARED ${sources})
target_link_libraries(GeantNuma  ${NUMA_LIBRARIES})

add_headers("${headers}")
install(TARGETS GeantNuma DESTINATION lib)
install(FILES DESTINATION lib)
