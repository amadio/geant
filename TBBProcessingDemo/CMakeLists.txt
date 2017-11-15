# See
# http://stackoverflow.com/questions/1027247/best-way-to-specify-sourcefiles-in-cmake
# including the discussion.
include_directories(${TBB_INCLUDE_DIRS})

#----------------------------------------------------------------------------------------------
# Find GeantV headers library
#
#find_package(GeantV REQUIRED)
if(GeantV_DIR)
  message(STATUS "GeantV includes: ${GeantV_DIR}/inc")
  include_directories(AFTER SYSTEM ${GeantV_DIR}/inc)
endif()

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Geom Graf3d RIO Core Physics EG)
if(ROOT_FOUND)
  message(STATUS "ROOT found: ${ROOT_FOUND}")
  message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
endif()

#----------------------------------------------------------------------------------------------
# Find VecGeom geometry headers library
#

message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
find_package(VecGeom REQUIRED)
if( VECGEOM_FOUND )
   find_package(Vc REQUIRED)
   if ( Vc_FOUND )
      include_directories( AFTER SYSTEM ${Vc_INCLUDE_DIR})
      set(VC_LIBRARIES ${Vc_LIBRARIES})
   endif()
   if( USE_VECGEOM_NAVIGATOR )
     add_definitions( "-DUSE_VECGEOM_NAVIGATOR" )
   endif()
   include_directories(${VECGEOM_INCLUDE_DIR})
else()
   message(FATAL_ERROR "Required VecGeom library not found")
endif()

if(ROOT_FOUND)
  include_directories(AFTER SYSTEM ${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})
  set(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} ${ROOT_LIBRARIES})
endif()

find_package(HepMC)
if(HepMC_FOUND)
  message(STATUS "HepMC_DIR: ${HepMC_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
  set(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} ${HEPMC_LIBRARIES})
endif()

include_directories("/Users/lima/work/geantv/base/inc")
include_directories("/Users/lima/work/geantv/xsec/inc")
include_directories("/Users/lima/work/geantv/vecprot_v2/inc")
include_directories("/Users/lima/work/geantv/vecprot_v2_tbb/inc")
include_directories("/Users/lima/work/geantv/examples/inc")
set(GEANTV_LIBRARIES ${GEANTV_LIBRARIES} ${GeantV_DIR}/lib/libGeant_v.so ${GeantV_DIR}/lib/libGeant_tbb.so ${GeantV_DIR}/lib/libGeantExamples.so ${GeantV_DIR}/lib/libXsec.so)

include_directories("/opt/local/vecgeom/install/include")

set(CORE_FILES
    TBBFrameworkCore/Event.cpp
    TBBFrameworkCore/EventProcessor.cpp
    TBBFrameworkCore/FactoryManagers.cpp
    TBBFrameworkCore/Filter.cpp
    TBBFrameworkCore/FilterOnPathWrapper.cpp
    TBBFrameworkCore/FilterWrapper.cpp
    TBBFrameworkCore/Getter.cpp
    TBBFrameworkCore/Module.cpp
    TBBFrameworkCore/ModuleThreadStack.cpp
    TBBFrameworkCore/ModuleWrapper.cpp
    TBBFrameworkCore/Path.cpp
    TBBFrameworkCore/PrefetchAndWorkWrapper.cpp
    TBBFrameworkCore/Producer.cpp
    TBBFrameworkCore/ProducerWrapper.cpp
    TBBFrameworkCore/Queues.cpp
    TBBFrameworkCore/Schedule.cpp
    TBBFrameworkCore/SerialTaskQueue.cpp
    TBBFrameworkCore/SimpleSource.cpp
    TBBFrameworkCore/WaitingList.cpp
    TBBFrameworkCore/WaitingTaskList.cpp)

include_directories(BEFORE
    ${CMAKE_CURRENT_SOURCE_DIR}/TBBFrameworkCore)
add_library(TBBFrameworkCore SHARED ${CORE_FILES})
target_link_libraries(TBBFrameworkCore ${TBB_LIBRARIES} ${GEANTV_LIBRARIES} ${VECGEOM_LIBRARIES})

set(TEST_FILES
  TBBTestModules/BusyWaitPassFilter.cpp
  TBBTestModules/BusyWaitProducer.cpp
  TBBTestModules/EventGeneratorProducer.cpp
  TBBTestModules/EventTimesBusyWaitPassFilter.cpp
  TBBTestModules/EventTimesBusyWaitProducer.cpp
  TBBTestModules/EventTimesPassFilterBase.cpp
  TBBTestModules/EventTimesProducerBase.cpp
  TBBTestModules/EventTimesSleepingPassFilter.cpp
  TBBTestModules/EventTimesSleepingProducer.cpp
  TBBTestModules/EventTimesSyncGetPassFilterBase.cpp
  TBBTestModules/EventTimesSyncGetProducerBase.cpp
  TBBTestModules/EventTimesSyncGetSleepingPassFilter.cpp
  TBBTestModules/EventTimesSyncGetSleepingProducer.cpp
  TBBTestModules/FastPassFilter.cpp
  TBBTestModules/FastProducer.cpp
  TBBTestModules/GeantVProducer.cpp
  TBBTestModules/CMSApplicationTBB.cpp
  TBBTestModules/PassFilter.cpp
  TBBTestModules/SleepingPassFilter.cpp
  TBBTestModules/SleepingProducer.cpp
  TBBTestModules/ThreadSafeProducer.cpp
  TBBTestModules/ThreadSafeWhichReadsProducer.cpp
  TBBTestModules/ThreadSaferGetterFilter.cpp
  TBBTestModules/ThrowingFilter.cpp
  TBBTestModules/ThrowingProducer.cpp
  TBBTestModules/Waiter.cpp
  TBBTestModules/thread_type_from_config.cpp
  ${CMAKE_SOURCE_DIR}/BusyWaitCalibration/busyWait.cpp
  ${CMAKE_SOURCE_DIR}/BusyWaitCalibration/busy_wait_scale_factor.cpp)

include_directories(BEFORE
    ${CMAKE_SOURCE_DIR}/BusyWaitCalibration)
add_library(TBBTestModules SHARED ${TEST_FILES})
target_link_libraries(TBBTestModules TBBFrameworkCore)

if(CPPUNIT_FOUND)
  include_directories(${CPPUNIT_INCLUDE_DIR})
  set(UNIT_TEST_FILES
      TBBUnitTests/Event_test.cpp
      TBBUnitTests/SerialTaskQueue_test.cpp
      TBBUnitTests/WaitingTaskList_test.cpp
      TBBUnitTests/main.cpp)
  add_executable(TBBRunUnitTests ${UNIT_TEST_FILES})
  target_link_libraries(TBBRunUnitTests TBBFrameworkCore ${TBB_LIBRARIES} ${CPPUNIT_LIBRARY} pthread)
endif()

add_executable(TBBDemo main.cpp)
target_link_libraries(TBBDemo TBBFrameworkCore TBBTestModules ${TBB_LIBRARIES} ${VECGEOM_LIBRARIES})
