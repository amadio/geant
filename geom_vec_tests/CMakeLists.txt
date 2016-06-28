# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

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

# for optimization studies, we compile this directory with higher Ox when no debug build
# TODO: check compiler etc
if(CMAKE_BUILD_TYPE MATCHES "Debug")
	set( BASEFLAGS "-O0 -fpermissive" )	
else()
# set baseflags as a function of the compiler:
      message(STATUS "setting up compilation flags to compile vectorized geometry")
 if( CMAKE_COMPILER_IS_GNUCXX )
      message(STATUS "configuration compile flags for gcc")
      set( BASEFLAGS "-Ofast -ftree-vectorize -ftree-vectorizer-verbose=0 -fpermissive" )
  elseif( CMAKE_CXX_COMPILER_ID STREQUAL Intel )
      message(STATUS "configuring compile flags for Intel compiler")
      set( BASEFLAGS "-Ofast -restrict -vec-report=2" )
  endif()
endif()

if( VC_FOUND )
  add_definitions( "-DVEC_EXTENSIONS" )
  set( BASEVCFLAGS "-std=c++11 -fabi-version=6" )
endif()
if( USEAVX )
  	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} ${BASEVCFLAGS} -mavx")
else()
   	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} ${BASEVCFLAGS} -msse4.1")
endif()

# for the timer
if ( USE_RDTSC )
  add_definitions( "-DUSE_RDTSC" )
endif()


#----------------------------------------------------------------------------
# Add the library for vectorized Geometries
#
if( VC_FOUND )
add_library(VecGeom SHARED ${sources} ${headers} ${VC_LIBRARIES})
add_library(VecGeom_static STATIC ${sources} ${headers} ${VC_LIBRARIES})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom -L${VC_LIBRARY_DIR} -lVc)
target_link_libraries(VecGeom_static ${ROOT_LIBRARIES} -lEG -lGeom -L${VC_LIBRARY_DIR} -lVc)
else()
add_library(VecGeom SHARED ${sources} ${headers})
add_library(VecGeom_static STATIC ${sources} ${headers})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom)
endif()

add_custom_target(VectorizedGeom DEPENDS VecGeom VecGeom_static)


#----------------------------------------------------------------------------
# Add some benchmarks and standard tests if wanted
if( VECGEOBENCHMARKS )
    file(GLOB benchmarksources src/BenchmarkFramework/BenchTGeoBBox_v.cxx
			       src/BenchmarkFramework/BenchTGeoTube_v.cxx
			       src/BenchmarkFramework/BenchTGeoTubeSeg_v.cxx
			       src/BenchmarkFramework/BenchTGeoPcon_v.cxx
			       src/BenchmarkFramework/BenchTGeoPconSeg_v.cxx
			       src/BenchmarkFramework/BenchTGeoCone_v.cxx
			       src/BenchmarkFramework/BenchTGeoConeSeg_v.cxx)
    if( VC_FOUND ) # if we have VC we build some additional stuff
	 list(APPEND benchmarksources src/BenchmarkFramework/BenchTGeoBBox_v_masked.cxx)
    endif( VC_FOUND )

    set_source_files_properties(${benchmarksources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} ${BASEVCFLAGS}")

    foreach( benchmarksrc  ${benchmarksources} )
        get_filename_component( BENCHEXEC  ${benchmarksrc} NAME_WE)
	message(STATUS "building benchmark ${BENCHEXEC}")
	if( VC_FOUND )
	add_executable( ${BENCHEXEC} ${benchmarksrc} ${headers} ${VC_LIBRARIES} ${VectorizedGeom})
	target_link_libraries( ${BENCHEXEC} ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} VecGeom_static  ${VC_LIBRARIES} ${TBB_LIBRARIES} -lrt)
	else()
	add_executable( ${BENCHEXEC} ${benchmarksrc} ${headers} )
	target_link_libraries( ${BENCHEXEC} ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} VecGeom_static ${TBB_LIBRARIES} -lrt)
	endif( VC_FOUND )
	add_dependencies( ${BENCHEXEC} VectorizedGeom )
    endforeach( benchmarksrc )



    # COPY THE CMS DETECTORFILE TO BINARY BUILD DIRECTORY (FOR CONVENIENCE)
    #---------------------------------------------------------------------------------------------
    configure_file(src/BenchmarkFramework/cms.root ${CMAKE_BINARY_DIR}/bin/cms.root COPYONLY)
endif()


#----------------------------------------------------------------------------

 
