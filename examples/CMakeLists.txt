#----------------------------------------------------------------------------
# Setup the project

include_directories(inc ../vecprot_v2/inc)

set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS "${GEANT5_COMPILE_DEFINITIONS}")
set(CMAKE_CXX_FLAGS                "${GEANT5_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANT5_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANT5_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANT5_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${GEANT5_EXE_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
set(sources
   ${CMAKE_CURRENT_SOURCE_DIR}/src/ExN03Application.cxx
   ${CMAKE_CURRENT_SOURCE_DIR}/src/CMSApplication.cxx
)     
     
set(headers
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/ExN03Application.h 
   ${CMAKE_CURRENT_SOURCE_DIR}/inc/CMSApplication.h 
)   

#file(GLOB sources src/*.cxx)

#file(GLOB headers inc/*.h)
#string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h" "" headers "${headers}")

ROOT_GENERATE_DICTIONARY(G__Examples ${headers} MODULE libGeantExamples LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
#

#----------------------------------------------------------------------------
# CUSTOMIZE COMPILATION IF VECGEOM IS USED

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
  -lGeant_v ${ROOT_LIBRARIES} -lEG -lGeom )

add_dependencies(GeantExamples Geantlib_v)

add_custom_target(Exampleslib DEPENDS GeantExamples)
