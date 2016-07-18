# Add this directrory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

link_libraries(vecphys)
include_directories(${CMAKE_SOURCE_DIR}/vecphys/inc
  ${CMAKE_SOURCE_DIR}/physics/particle/inc)

if (CUDA)
  link_libraries(vecphyscuda)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
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
include_directories(inc)
include_directories(${CMAKE_SOURCE_DIR}/xsec/inc)

#----------------------------------------------------------------------------
# Generate the root dictionary: this library will be used in run.C root macro.  
#----------------------------------------------------------------------------
if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(vphysprocDict ${headers} MODULE Vphysproc LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")

  #----------------------------------------------------------------------------
  # Add traget library and dependencies
  #----------------------------------------------------------------------------
  add_library(Vphysproc SHARED ${sources} vphysprocDict.cxx)
  target_link_libraries(Vphysproc ${VC_LIBRARIES} Xsec Geant_v ${ROOT_LIBRARIES})
else()
  add_library(Vphysproc SHARED ${sources})
  target_link_libraries(Vphysproc ${VC_LIBRARIES} Xsec Geant_v)
endif()

add_custom_target(VectPhysProc DEPENDS Vphysproc)

install(TARGETS Vphysproc DESTINATION lib)

install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc_rdict.pcm DESTINATION lib)

