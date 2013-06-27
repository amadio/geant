#----------------------------------------------------------------------------
# Cmake for vectorized geometry part 
#
include_directories(inc)
include_directories(${TBB_INCLUDE_DIR})

# depending on whether we use Vc, have to set some compile flags
if( VC_FOUND )
    if( USEAVX )
    	add_definitions(-mavx -std=c++11 -fabi-version=6 -DVEC_EXTENSIONS)
    else()
        add_definitions(-msse4.1 -std=c++11 -fabi-version=6 -DVEC_EXTENSIONS)
    endif()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/TGeoBBox_v.cxx)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)


#----------------------------------------------------------------------------
# Add the library for vectorized Geometries
#
add_library(VecGeom SHARED ${sources} ${headers} ${VC_LIBRARIES})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom)

add_custom_target(VectorizedGeom DEPENDS VecGeom)

#----------------------------------------------------------------------------
# Add some benchmarks if wanted
#
if( VECGEOBENCHMARKS )
   file(GLOB benchsource1 src/BenchmarkFramework/BenchTGeoBBox_v.cxx)
   add_executable(BenchTGeoBBox_v ${benchsource1} ${headers} ${VC_LIBRARIES})
   add_dependencies(BenchTGeoBBox_v VectorizedGeom)
   target_link_libraries(BenchTGeoBBox_v ${ROOT_LIBRARIES} -lGeom ${VC_LIBRARIES} ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom)
endif()


#----------------------------------------------------------------------------
