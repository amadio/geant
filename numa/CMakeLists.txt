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
# Find if the system is NUMA aware
#----------------------------------------------------------------------------
set(NUMA_ALL_LIBRARIES "-lGeantNuma")
if (USE_NUMA)
  find_package(Numa)
  if (NUMA_FOUND)
    message(STATUS "Found Numa: "  ${NUMA_INCLUDE_DIR})
    # Check the presence of hwloc library
    find_package(hwloc)
    if (HWLOC_FOUND)
      add_definitions( "-DUSE_NUMA" )
      include_directories( AFTER SYSTEM ${NUMA_INCLUDE_DIR})
      set(NUMA_LIBRARIES -L${NUMA_LIBRARY_DIR} -lnuma ${HWLOC_LIBRARIES})
      set(NUMA_ALL_LIBRARIES "-lGeantNuma ${NUMA_LIBRARIES}")    
    else()
      message(STATUS "Hardware locality (hwloc) libraries not found. Disabling NUMA.")
      set(USE_NUMA OFF)
    endif()
  else()
    message(STATUS "Numa libraries (libnuma, libhwloc) not found.")
    set(USE_NUMA OFF)
  endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
set(sources 
   src/NumaNode.cxx
   src/NumaPolicy.cxx
   src/NumaTopology.cxx
   src/NumaUtils.cxx
)

set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/GeantNuma.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaBlock.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaBlockMgr.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaNode.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaPolicy.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaTopology.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaUtils.h
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/NumaAllocator.h
)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(GeantNuma SHARED ${sources})
target_link_libraries(GeantNuma  ${NUMA_LIBRARIES})

#----------------------------------------------------------------------------
# Headers need to be exported to include directory, explicit paths needed
#----------------------------------------------------------------------------
add_headers("${headers}")

install(TARGETS GeantNuma DESTINATION lib)
