#----------------------------------------------------------------------------
# Cmake for vectorized geometry part 
#
include_directories(inc)

find_package(Vc)
if ( VC_FOUND )
include_directories(${VC_INCLUDE_DIR})
endif()

find_package(TBB)
include_directories(${TBB_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/TGeoBBox_v.cxx src/TGeoPcon_v.cxx src/TGeoTube_v.cxx src/TGeoCone_v.cxx )
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)

# for optimization studies, we compile this directory with higher Ox
# TODO: check compiler etc
set( BASEFLAGS "-Ofast -ftree-vectorize -ftree-vectorizer-verbose=0" )

# depending on whether we use Vc, have to set some additional compile flags
if( VC_FOUND )
  set_source_files_properties(${sources} PROPERTIES COMPILE_DEFINITIONS VEC_EXTENSIONS)
  set( BASEVCFLAGS "-std=c++11 -fabi-version=6" )
endif()
if( USEAVX )
  	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} ${BASEVCFLAGS} -mavx")
else()
   	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} ${BASEVCFLAGS} -msse4.1")
endif()



#----------------------------------------------------------------------------
# Add the library for vectorized Geometries
#
if( VC_FOUND )
add_library(VecGeom SHARED ${sources} ${headers} ${VC_LIBRARIES})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom -L${VC_LIBRARY_DIR} -lVc)
else()
add_library(VecGeom SHARED ${sources} ${headers})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom)
endif()

add_custom_target(VectorizedGeom DEPENDS VecGeom)

#----------------------------------------------------------------------------
# Add some benchmarks and standard tests if wanted
# TODO: make this general ( with patterns, etc. )
if( VECGEOBENCHMARKS )
    file(GLOB benchsource1 src/BenchmarkFramework/BenchTGeoBBox_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoBBox_v ${benchsource1} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoBBox_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoBBox_v ${benchsource1} ${headers})
        target_link_libraries(BenchTGeoBBox_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoBBox_v VectorizedGeom)

    file(GLOB benchsourcepcon src/BenchmarkFramework/BenchTGeoPcon_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoPcon_v ${benchsourcepcon} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoPcon_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoPcon_v ${benchsourcepcon} ${headers})
        target_link_libraries(BenchTGeoPcon_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoPcon_v VectorizedGeom)

   file(GLOB benchsourceTube src/BenchmarkFramework/BenchTGeoTube_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoTube_v ${benchsourceTube} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoTube_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoTube_v ${benchsourceTube} ${headers})
        target_link_libraries(BenchTGeoTube_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoTube_v VectorizedGeom)

   file(GLOB benchsourceTubeSeg src/BenchmarkFramework/BenchTGeoTubeSeg_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoTubeSeg_v ${benchsourceTubeSeg} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoTubeSeg_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoTubeSeg_v ${benchsourceTubeSeg} ${headers})
        target_link_libraries(BenchTGeoTubeSeg_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoTubeSeg_v VectorizedGeom)

   file(GLOB benchsourceCone src/BenchmarkFramework/BenchTGeoCone_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoCone_v ${benchsourceCone} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoCone_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoCone_v ${benchsourceCone} ${headers})
        target_link_libraries(BenchTGeoCone_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoCone_v VectorizedGeom)


   file(GLOB benchsourceConeSeg src/BenchmarkFramework/BenchTGeoConeSeg_v.cxx)	
    if( VC_FOUND )
        add_executable(BenchTGeoConeSeg_v ${benchsourceConeSeg} ${headers} ${VC_LIBRARIES})
        target_link_libraries(BenchTGeoConeSeg_v ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so ${VC_LIBRARIES} ${TBB_LIBRARIES})
    else()				      
        add_executable(BenchTGeoConeSeg_v ${benchsourceConeSeg} ${headers})
        target_link_libraries(BenchTGeoConeSeg_v ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom -Wl,-rpath ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libVecGeom.so)
    endif()
    add_dependencies(BenchTGeoConeSeg_v VectorizedGeom)


# for more complete test from CMS:
   file(GLOB benchsource2 src/BenchmarkFramework/BenchSelectedCMSShapes.cxx)
   add_executable(BenchSelectedCMSShapes ${benchsource2} ${headers})
   target_link_libraries(BenchSelectedCMSShapes ${ROOT_LIBRARIES} -lGeom ${TBB_LIBRARIES}) 
endif()


#----------------------------------------------------------------------------
 
