#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources 
  ${CMAKE_SOURCE_DIR}/physics/neutron/nudy/src/*.cxx
)
file(GLOB headers 
  ${CMAKE_SOURCE_DIR}/physics/neutron/nudy/inc/Geant/T*.h
)

include_directories(
  ${CMAKE_SOURCE_DIR}/physics/neutron/nudy/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
)

if(USE_ROOT)
ROOT_GENERATE_DICTIONARY(NudyDict ${headers} MODULE Nudy LINKDEF inc/Geant/LinkDef.h)
set(sources ${sources} NudyDict.cxx)
endif()
#----------------------------------------------------------------------------------------------
# Add the executable
#
add_library(Nudy SHARED ${sources} NudyDict.cxx)
target_link_libraries(Nudy -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${ROOT_LIBRARIES} MathMore 
                      ${VECGEOM_LIBRARIES}  ${Vc_LIBRARIES}
                      ${Geant4_LIBRARIES})
#----------------------------------------------------------------------------------------------
# Install the files
#
install(FILES ${headers} DESTINATION inc/Geant)

if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Nudy DESTINATION lib)
endif()

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy.rootmap 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy_rdict.pcm DESTINATION lib)
endif()

#---------------------------------------------------------------------------------------------
# tests to create root file from endf, reconstruction and retriving data after reconstruction
add_subdirectory(EndfToRoot) 
# compile standalone
#
add_subdirectory(standalone)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

