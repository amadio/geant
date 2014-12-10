#----------------------------------------------------------------------------
# Find VecPhys (required )
#----------------------------------------------------------------------------
find_package(VecGeom REQUIRED)
if(VECGEOM_FOUND)
  include_directories(${VECGEOM_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Required library VecGeom was not found!")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
set(headers 
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ConstFieldHelixStepper.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/GVIntegratorDriver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/TClassicalRK4.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/TIntegratorDriver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/TMagFieldEquation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/TSimpleRunge.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/TUniformMagField.h
)

set(sources
  ${CMAKE_CURRENT_SOURCE_DIR}/src/FieldSetup.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/circlestepper.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestConstBzStepper.cxx
)

#----------------------------------------------------------------------------
# Set include directories
#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
#include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

# print include directories
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include directories: ${inc_dirs}")

#----------------------------------------------------------------------------
# Generate the root dictionary: this library will be used in run.C root macro.  
#----------------------------------------------------------------------------
## ROOT_GENERATE_DICTIONARY(vMagFieldDict ${headers} MODULE libVphysproc LINKDEF inc/LinkDef.h)

#----------------------------------------------------------------------------
# Add traget library and dependencies
#----------------------------------------------------------------------------
## add_library(Vmagfield SHARED ${sources} vMagFieldDict.cxx)
target_link_libraries(Vmagfield  ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_custom_target(VectMagField DEPENDS Vmagfield)

