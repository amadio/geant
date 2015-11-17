#----------------------------------------------------------------------------
# Setup the project

include_directories(inc ../vecprot_v2/inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_CURRENT_SOURCE_DIR}/src/ExN03Application.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/CMSApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/LHCbApplication.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/FastSimApplication.cxx
)     
     
set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/ExN03Application.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/CMSApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/LHCbApplication.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/FastSimApplication.h 
)   

#file(GLOB sources src/*.cxx)

#file(GLOB headers inc/*.h)
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h" "" headers "${headers}")

if (NOT MIC)
  ROOT_GENERATE_DICTIONARY(G__Examples ${headers} MODULE libGeantExamples LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
endif()
#----------------------------------------------------------------------------
#

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED
if (NOT MIC)
if (CUDA)
  CUDA_ADD_LIBRARY(GeantExamples
                 ${sources} G__Examples.cxx
                 SHARED
                 OPTIONS ${CUDA_ARCH}
  )
else()
  add_library(GeantExamples SHARED ${sources} G__Examples.cxx)
endif()

target_link_libraries(GeantExamples -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -lGeant_v ${ROOT_LIBRARIES} -lEG -lGeom ${Vc_LIBRARIES})

else()
add_library(GeantExamples SHARED ${sources})
target_link_libraries(GeantExamples -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -lGeant_v)
endif()
add_dependencies(GeantExamples Geant_v)

add_custom_target(Exampleslib DEPENDS GeantExamples)

if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples.so DESTINATION lib)
else()
  install(TARGETS GeantExamples DESTINATION lib)
endif()
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples_rdict.pcm DESTINATION lib)
