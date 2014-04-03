#----------------------------------------------------------------------------
# Setup the project

include_directories(inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/MyHit.cxx
)     
     
set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/MyHit.h 
)   

ROOT_GENERATE_DICTIONARY(G__User ${headers} LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
#
include_directories(${PROJECT_SOURCE_DIR}/vecprot_v2/inc)
add_library(User SHARED ${sources} G__User.cxx ${headers} G__User.h)
target_link_libraries(User ${ROOT_LIBRARIES} -lEG -lGeom)
target_link_libraries(User -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant_v)

add_custom_target(Userlib DEPENDS User)
add_dependencies(Userlib Geantlib_v)
