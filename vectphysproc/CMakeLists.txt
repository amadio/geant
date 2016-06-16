
link_libraries(vecphys)
include_directories(${CMAKE_SOURCE_DIR}/vecphys/inc)

if (CUDA)
  link_libraries(vecphyscuda)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
set(sources 
${CMAKE_CURRENT_SOURCE_DIR}/src/GVectorPhysicsProcess.cc
${CMAKE_CURRENT_SOURCE_DIR}/src/GVComptonProcess.cc
)
set(headers
${CMAKE_CURRENT_SOURCE_DIR}/inc/GVectorPhysicsProcess.h
${CMAKE_CURRENT_SOURCE_DIR}/inc/GVComptonProcess.h
)

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

# print include directories
#get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message("-----  ${inc_dirs}")

#----------------------------------------------------------------------------
# Generate the root dictionary: this library will be used in run.C root macro.  
#----------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(vphysprocDict ${headers} MODULE Vphysproc LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(Vphysproc SHARED ${sources} vphysprocDict.cxx)
target_link_libraries(Vphysproc ${VC_LIBRARIES} Xsec Geant_v ${ROOT_LIBRARIES} EG Geom)

add_custom_target(VectPhysProc DEPENDS Vphysproc)

install(TARGETS Vphysproc DESTINATION lib)
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc_rdict.pcm DESTINATION lib)
