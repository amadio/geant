#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cc)
file(GLOB headers inc/*.h ${CMAKE_SOURCE_DIR}/base/inc/PhysicalConstants.h ${CMAKE_SOURCE_DIR}/base/inc/SystemOfUnits.h)

include_directories(inc ${CMAKE_SOURCE_DIR}/base/inc ${VECGEOM_INCLUDE_DIR})

set(_extraLibs "")

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(materialDict ${headers} MODULE Material LINKDEF materialLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  set(sources ${sources} materialDict.cxx)
  set(_extraLibs "${ROOT_LIBRARIES}")
endif()

add_library(Material SHARED ${sources})

# will depend on VecGeom through vecgeom::Vector and vecgeom::map
target_link_libraries(Material ${VECGEOM_LIBRARIES} ${_extraLibs})

add_headers("${headers}")
add_custom_target(materiallib DEPENDS Material)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libMaterial${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Material DESTINATION lib)
endif()
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libMaterial.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libMaterial_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
  set(DOXYFILE_CITE_BIB_FILES "${DOXYFILE_CITE_BIB_FILES} ${CMANE_CURRENT_SOURCE_DIR}/material/doc/biblio" PARENT_SCOPE)
endif()
