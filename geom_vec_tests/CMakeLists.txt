#----------------------------------------------------------------------------
# Cmake for vectorized geometry part 
#
include_directories(inc)

find_package(Vc REQUIRED)
include_directories(${VC_INCLUDE_DIR})

find_package(TBB)
include_directories(${TBB_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/TGeoBBox_v.cxx)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)

# depending on whether we use Vc, have to set some compile flags
if( VC_FOUND )
  set_source_files_properties(${sources} PROPERTIES COMPILE_DEFINITIONS VEC_EXTENSIONS)
  set( BASEVCFLAGS "-std=c++11 -fabi-version=6" )
  if( USEAVX )
    	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEVCFLAGS} -mavx")
  else()
    	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEVCFLAGS} -msse4.1")
  endif()
endif()


#----------------------------------------------------------------------------
# Add the library for vectorized Geometries
#
add_library(VecGeom SHARED ${sources} ${headers} ${VC_LIBRARIES})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom -L${VC_LIBRARY_DIR} -lVc)

add_custom_target(VectorizedGeom DEPENDS VecGeom)

#----------------------------------------------------------------------------
# Add some benchmarks if wanted
#
if( VECGEOBENCHMARKS )
   file(GLOB benchsource1 src/BenchmarkFramework/BenchTGeoBBox_v.cxx)
   add_executable(BenchTGeoBBox_v ${benchsource1} ${headers} ${VC_LIBRARIES})
   add_dependencies(BenchTGeoBBox_v VectorizedGeom)
   target_link_libraries(BenchTGeoBBox_v ${ROOT_LIBRARIES} -lGeom ${VC_LIBRARIES} ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom)

# for more complete test from CMS:
   file(GLOB benchsource2 src/BenchmarkFramework/BenchSelectedCMSShapes.cxx)
   add_executable(BenchSelectedCMSShapes ${benchsource2} ${headers})
   target_link_libraries(BenchSelectedCMSShapes ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES}) 
endif()


#----------------------------------------------------------------------------
 
