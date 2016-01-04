#----------------------------------------------------------------------------------------------
# Add source files & include directories

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(inc ../vecprot_v2/inc ../physics/particle/inc ../numa/inc)

if (USE_VECGEOM_NAVIGATOR)
 include_directories(${CMAKE_SOURCE_DIR}/material/inc)
endif()

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()
include_directories(../magneticfield/inc  ../VecMagFieldRoutine/include)

set(sources
  src/ExN03Application.cxx
  src/CMSApplication.cxx
#  src/CMSDetectorConstruction.cxx
  src/LHCbApplication.cxx
  src/FastSimApplication.cxx
  )

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
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/ExN03Application.h
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/CMSApplication.h
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/CMSDetectorConstruction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/LHCbApplication.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/FastSimApplication.h 
  )   
if (USE_ROOT)
  if (USE_UMESIMD)
     ROOT_GENERATE_DICTIONARY(G__Examples ${headers} MODULE GeantExamples
     LINKDEF inc/LinkDef.h OPTIONS "-DVECGEOM_UMESIMD=1")
  else()
     ROOT_GENERATE_DICTIONARY(G__Examples ${headers} MODULE GeantExamples
     LINKDEF inc/LinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  endif()
    set(sources ${sources} G__Examples.cxx)
endif()

#----------------------------------------------------------------------------------------------
# CUDA section
#
if (CUDA)
  CUDA_ADD_LIBRARY(GeantExamples
    ${sources}
    SHARED
    OPTIONS ${CUDA_ARCH}
    )
else()
  add_library(GeantExamples SHARED ${sources})
endif()

target_link_libraries(GeantExamples -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v -lVmagfield
  ${TBB_LIBRARIES})

if (USE_ROOT)
  target_link_libraries(GeantExamples ${ROOT_LIBRARIES})
endif()

if (KNC)
  target_link_libraries(GeantExamples ${Vc_MIC_LIBRARIES})
endif()

add_dependencies(GeantExamples Geant_v)
add_custom_target(Exampleslib DEPENDS GeantExamples)

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE)
  # for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples.so DESTINATION lib)
else()
  install(TARGETS GeantExamples DESTINATION lib)
endif()

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libGeantExamples_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
