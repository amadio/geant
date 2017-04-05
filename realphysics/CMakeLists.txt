#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources dummies/src/*.cc physics/management/src/*.cxx physics/models/EM/src/*.cc
                  physics/processes/EM/src/*.cc physics/utils/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/base/inc/PhysicalConstants.h ${CMAKE_SOURCE_DIR}/base/inc/SystemOfUnits.h
                  dummies/inc/*.h physics/management/inc/*.h physics/models/EM/inc/*.h
                  physics/processes/EM/inc/*.h physics/utils/inc/*.h)

include_directories(${CMAKE_SOURCE_DIR}/material/inc ${CMAKE_SOURCE_DIR}/base/inc dummies/inc physics/management/inc physics/models/EM/inc physics/processes/EM/inc physics/utils/inc ${CMAKE_SOURCE_DIR}/numa/inc)

set(_extraLibs "")

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(realphysicsDict ${headers} MODULE RealPhysics LINKDEF realphysicsLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  set(sources ${sources} realphysicsDict.cxx)
  set(_extraLibs "${ROOT_LIBRARIES}")
endif()

add_library(RealPhysics SHARED ${sources})

target_link_libraries(RealPhysics -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material Geant_v ${_extra_libs}
                      ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})

add_headers("${headers}")
add_custom_target(realphyslib DEPENDS RealPhysics)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS RealPhysics DESTINATION lib)
endif()
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
  set(DOXYFILE_CITE_BIB_FILES "${DOXYFILE_CITE_BIB_FILES} ${CMANE_CURRENT_SOURCE_DIR}/doc/biblio" PARENT_SCOPE)
endif()
