# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
  set(DOXYFILE_CITE_BIB_FILES "${DOXYFILE_CITE_BIB_FILES} ${CMANE_CURRENT_SOURCE_DIR}/doc/biblio" PARENT_SCOPE)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers
#
file(GLOB sources management/src/*.cxx material/src/*.cxx particle/src/*.cxx)
file(GLOB headers management/inc/*.h material/inc/*.h particle/inc/*.h)

include_directories(management/inc material/inc particle/inc)

set(_extraLibs "")

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(physicsDict ${headers} MODULE NewPhysics LINKDEF physicsLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  set(sources ${sources} physicsDict.cxx)
  set(_extraLibs "${ROOT_LIBRARIES}")
endif()

add_library(NewPhysics SHARED ${sources})

target_link_libraries(NewPhysics -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} 
     ${_extra_libs} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})

add_headers("${headers}")

add_custom_target(physlib DEPENDS NewPhysics)

if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS NewPhysics DESTINATION lib)
endif()

install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics_rdict.pcm DESTINATION lib)


