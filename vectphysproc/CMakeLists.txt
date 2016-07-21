#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
set(sources 
  src/GVectorPhysicsProcess.cc
  src/GVComptonProcess.cc
  )
set(headers
  inc/GVectorPhysicsProcess.h
  inc/GVComptonProcess.h
  )

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(
  inc
  ${CMAKE_SOURCE_DIR}/xsec/inc
  ${CMAKE_SOURCE_DIR}/vecphys/inc
  ${CMAKE_SOURCE_DIR}/physics/particle/inc
  )

#----------------------------------------------------------------------------------------------
# Add the library
#

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(vphysprocDict ${headers} MODULE Vphysproc 
    LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  add_library(Vphysproc SHARED ${sources} vphysprocDict.cxx)
  target_link_libraries(Vphysproc -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VC_LIBRARIES} Xsec Geant_v vecphys ${ROOT_LIBRARIES})
else()
  add_library(Vphysproc SHARED ${sources})
  target_link_libraries(Vphysproc -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VC_LIBRARIES} Xsec Geant_v vecphys)
endif()

add_custom_target(VectPhysProc DEPENDS Vphysproc)

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS Vphysproc DESTINATION lib)
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc.rootmap 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

