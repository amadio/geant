#----------------------------------------------------------------------------
# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
set(sources 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/NumaNode.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/NumaCore.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/NumaPolicy.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/NumaTopology.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/NumaUtils.cxx
)

set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantNuma.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaBlock.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaBlockMgr.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaNode.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaCore.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaPolicy.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaTopology.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaUtils.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/Geant/NumaAllocator.h
)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(GeantNuma SHARED ${sources})
target_link_libraries(GeantNuma  ${NUMA_LIBRARIES} ${VECGEOM_LIBRARIES})

install(FILES ${headers} DESTINATION inc/Geant)
install(TARGETS GeantNuma DESTINATION lib)
