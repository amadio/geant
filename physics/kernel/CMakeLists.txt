#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources particle/src/*.cxx particle/src/generated/*.cxx)
file(GLOB headers particle/inc/*.h)

include_directories(particle/inc)
include_directories(${CMAKE_SOURCE_DIR}/numa/inc)

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

#----------------------------------------------------------------------------------------------
# CUDA section
#
if (CUDA)
  foreach(SRC_FILE ${sources})
    if (NOT "${SRC_FILE}" MATCHES "Dict.cxx$")
      get_filename_component(SRC_FILENAME ${SRC_FILE} NAME_WE)

      ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC_FILE}
        ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
        DEPENDS ${SRC_FILE}
        )
      set(
        SRC_CUDA ${SRC_CUDA}
        ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
	)
    endif()
  endforeach()

  cuda_add_library(NewPhysicsCuda
    ${SRC_CUDA}
    SHARED
    OPTIONS ${CUDA_ARCH}
    )

  cuda_add_library_depend(NewPhysicsCuda ${VECGEOM_CUDA_STATIC_LIBRARY} )
  target_link_libraries(NewPhysicsCuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL})
endif()

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS NewPhysics DESTINATION lib)
endif()
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNewPhysics_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
  set(DOXYFILE_CITE_BIB_FILES "${DOXYFILE_CITE_BIB_FILES} ${CMANE_CURRENT_SOURCE_DIR}/doc/biblio" PARENT_SCOPE)
endif()
