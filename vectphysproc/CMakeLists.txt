
#----------------------------------------------------------------------------
# Find VecPhys (required )
#----------------------------------------------------------------------------
find_package(VecPhys REQUIRED)
if(VECPHYS_FOUND)
  include_directories(${VECPHYS_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Required library vecphys was not found!")
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
ROOT_GENERATE_DICTIONARY(vphysprocDict ${headers} MODULE libVphysproc LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")


#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
add_library(Vphysproc SHARED ${sources} vphysprocDict.cxx)
target_link_libraries(Vphysproc  ${VECPHYS_LIBRARIES} ${VC_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v)

add_custom_target(VectPhysProc DEPENDS Vphysproc)

install(TARGETS Vphysproc DESTINATION lib)
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVphysproc_rdict.pcm DESTINATION lib)
