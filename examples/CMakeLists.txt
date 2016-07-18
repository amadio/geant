# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

include_directories(inc ../vecprot_v2/inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   src/ExN03Application.cxx
   src/CMSApplication.cxx
   src/LHCbApplication.cxx
   src/FastSimApplication.cxx
)     
     
set(headers
   inc/ExN03Application.h 
   inc/CMSApplication.h 
   inc/LHCbApplication.h 
   inc/FastSimApplication.h 
)   
if (USE_ROOT)
ROOT_GENERATE_DICTIONARY(G__Examples ${headers} MODULE GeantExamples LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
endif()
#----------------------------------------------------------------------------
#

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED
if (NOT KNC)
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
  -lGeant_v ${ROOT_LIBRARIES} ${Vc_LIBRARIES})

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

