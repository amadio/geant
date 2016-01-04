cmake_minimum_required(VERSION 2.8.8)
project(Hello)

#Added to detect compiler support for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


#if (USE_ANANYA)
 # add_definitions(-DANANYA)
#endif (USE_ANANYA)

#Copies cms2015.txt to build/VecMagFieldRoutine at runtime 
configure_file(${CMAKE_SOURCE_DIR}/VecMagFieldRoutine/cms2015.txt cms2015.txt COPYONLY)

#Added to include root. Used for graph plotting.
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

include_directories(include)
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(magfieldtest ${CMAKE_CURRENT_SOURCE_DIR}/test/magfieldtest.cpp ${sources})
target_link_libraries(magfieldtest ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(testVector3D ${CMAKE_CURRENT_SOURCE_DIR}/test/testVector3D.cpp)
target_link_libraries(testVector3D ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(MainTestFunc ${CMAKE_CURRENT_SOURCE_DIR}/test/MainMagFieldTest.cpp ${sources})
target_link_libraries(MainTestFunc ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(BoundariesGraph ${CMAKE_CURRENT_SOURCE_DIR}/test/visualizations/GraphsForBoundaries.cpp ${sources})
target_link_libraries(BoundariesGraph ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(TimingTest ${CMAKE_CURRENT_SOURCE_DIR}/test/BenchmarkTiming.cpp ${sources})
target_link_libraries(TimingTest ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(TestFunc ${CMAKE_CURRENT_SOURCE_DIR}/test/MainMagFieldTest2.cpp ${sources})
target_link_libraries(TestFunc ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})
