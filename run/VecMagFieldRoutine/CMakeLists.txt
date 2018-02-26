cmake_minimum_required(VERSION 2.8.8)
project(TestVecMagField)

#Added to include root. Used for graph plotting.
find_package(ROOT)
include_directories(${VECGEOM_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

include_directories(${CMAKE_SOURCE_DIR}/run/magneticfield/inc)

#Copies cms2015.txt to build/VecMagFieldRoutine at runtime 
configure_file(${CMAKE_SOURCE_DIR}/run/VecMagFieldRoutine/cmsmagfield2015.txt cmsmagfield2015.txt COPYONLY)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmsmagfield2015.txt DESTINATION VecMagFieldRoutine)

include_directories(inc)
include_directories(${CMAKE_SOURCE_DIR}/run/magneticfield/inc)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

set( ALL_LIBRARIES Vmagfield ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(Vecmagfieldtest ${CMAKE_CURRENT_SOURCE_DIR}/test/magfieldtest.cpp ${sources})
#target_link_libraries(Vecmagfieldtest ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(smallTest ${CMAKE_CURRENT_SOURCE_DIR}/test/smallTest.cpp ${sources})
#target_link_libraries(smallTest ${ALL_LIBRARIES}) # ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(VecTimingTest ${CMAKE_CURRENT_SOURCE_DIR}/test/BenchmarkTiming.cpp ${sources})
target_link_libraries(VecTimingTest ${ALL_LIBRARIES}) # ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})
install(TARGETS VecTimingTest DESTINATION bin)

add_executable(VectorAgainstScalar ${CMAKE_CURRENT_SOURCE_DIR}/test/CheckVectorAgainstScalar.cpp ${sources})
target_link_libraries(VectorAgainstScalar  ${ALL_LIBRARIES}) # ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})
install(TARGETS VectorAgainstScalar DESTINATION bin)

add_executable(VecTestFunc ${CMAKE_CURRENT_SOURCE_DIR}/test/MainMagFieldTest.cpp ${sources})
target_link_libraries(VecTestFunc  ${ALL_LIBRARIES}) #  ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(VecFloatDoubleDiff ${CMAKE_CURRENT_SOURCE_DIR}/test/FloatDoubleDiff.cpp ${sources})
#target_link_libraries(VecFloatDoubleDiff ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(VecTestGatherAlias ${CMAKE_CURRENT_SOURCE_DIR}/test/TestGatherAlias.cpp ${sources})
#target_link_libraries(VecTestGatherAlias ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

add_executable(TestAOS3D ${CMAKE_CURRENT_SOURCE_DIR}/test/TestAOS3D.cpp ${sources})
target_link_libraries(TestAOS3D  ${ALL_LIBRARIES}) # ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(VecTestVector3D ${CMAKE_CURRENT_SOURCE_DIR}/test/testVector3D.cpp)
#target_link_libraries(VecTestVector3D ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})

#add_executable(VecBoundariesGraph ${CMAKE_CURRENT_SOURCE_DIR}/test/visualizations/GraphsForBoundaries.cpp ${sources})
#target_link_libraries(VecBoundariesGraph ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})


# add_executable(smallTest ${CMAKE_CURRENT_SOURCE_DIR}/test/smallTest.cpp ${sources})
# target_link_libraries(smallTest ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL} ${USERKERNELLIB})
