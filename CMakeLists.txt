 #----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GEANTV)

# Fix bug on APPLE, this is the default everywhere
if(APPLE)
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif()

#----------------------------------------------------------------------------
# Additional options

option(CUDA "Enable compilation for CUDA." OFF)

option(DATA_DOWNLOAD "Enable data download" OFF)

option(CTEST "Enable CTest when building." OFF)
if (CTEST)
  set (DATA_DOWNLOAD ON)
endif()

option(USE_VECPHYS "Enable activation of vector physics." OFF)
if(USE_VECPHYS)
  add_definitions(-DUSE_VECPHYS)
  message(STATUS "USE_VECPHYS ON  : Vector physics is active !")
else()
  message(STATUS "USE_VECPHYS OFF : Vector physics is NOT active !")
endif()


#---Set pathes where to put the libraries, executables and headers------------------------------
if( GEANT_OUTPUT_DIR )
    set(OUTDIR ${GEANT_OUTPUT_DIR})
else()
    set(OUTDIR $ENV{GEANT_OUTPUT_DIR})
endif()
if(NOT OUTDIR)
  set(OUTDIR ${PROJECT_SOURCE_DIR})
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
#----------------------------------------------------------------------------
# Find VecGeom geometry headers library 
#
find_package(VecGeom REQUIRED)
if( VECGEOM_FOUND )
    find_package(Vc)
    if ( VC_FOUND )
        include_directories( AFTER SYSTEM ${VC_INCLUDE_DIR})
    endif()
    if( USE_VECGEOM_NAVIGATOR )
     add_definitions( "-DUSE_VECGEOM_NAVIGATOR" )
     endif() 
   include_directories(${VECGEOM_INCLUDE_DIR})
else()
   message(FATAL_ERROR "Required VecGeom library not found")
endif()

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
find_package(ROOT REQUIRED)

include_directories(AFTER SYSTEM ${ROOT_INCLUDE_DIR})

if( USE_ROOT )
   add_definitions("-DUSE_ROOT")
endif()

#----------------------------------------------------------------------------
# Load macros
#
include(cmake/modules/GeantMacros.cmake)

#---------------------------------------------------------------------------
# Include creation and submission to dashboard cdash.cern.ch

include(CTest)

if (CTEST)
  message(STATUS "Testing with CTest enabled.")
  enable_testing()
endif()

#---------------------------------------------------------------------------
# Set data to be dowloaded

 function(DATA_DOWNLOAD)
    set(DATA_FILE data.tar.gz)
    if(APPLE)
      execute_process(COMMAND curl -o ${DATA_FILE} ${DATA_URL})      
    else()  
      execute_process(COMMAND wget -q ${DATA_URL} -O ${DATA_FILE})
    endif()
 endfunction(DATA_DOWNLOAD)

message(STATUS "Downloading data files")
if(DATA_DOWNLOAD)
  set(DATA_URL "http://mnovak.web.cern.ch/mnovak/data/data.tar.gz")
  set(EXPECTED_MD5SUM "8de0246795899f373c08c721161b8d42")
message(STATUS "Searching for already downloaded data file(s)")
find_file(DATA_FILE data.tar.gz ${CMAKE_BINARY_DIR})
  if(DATA_FILE STREQUAL "DATA_FILE-NOTFOUND")
    message("File not found: data.tar.gz")
    DATA_DOWNLOAD()
  else()
    message("File found:" ${DATA_FILE})
    if(APPLE)
      execute_process(COMMAND md5sum ${DATA_FILE} OUTPUT_VARIABLE MD5SUM)
      string(LENGTH ${MD5SUM} MD5LENGTH)
      MATH(EXPR START "${MD5LENGTH} - 33")
      string(SUBSTRING ${MD5SUM} ${START} 32 MD5SUM)         
    else()
      execute_process(COMMAND md5sum ${DATA_FILE} OUTPUT_VARIABLE MD5SUM)
      string(SUBSTRING ${MD5SUM} 0 32 MD5SUM)
    endif()
    if(MD5SUM STREQUAL EXPECTED_MD5SUM)
      message("MD5sum of file is matching")
    else()
      message("MD5sum of file is not matching " ${EXPECTED_MD5SUM})
      DATA_DOWNLOAD()
    endif()
  endif()

message("Untar data.tar.gz")
execute_process(COMMAND tar xzf ${CMAKE_BINARY_DIR}/data.tar.gz -C ${CMAKE_SOURCE_DIR}/vecprot_v2/tests/)
#execute_process(COMMAND rm -rf ${CMAKE_BINARY_DIR}/data.tar.gz)
endif()

################################################################################

#----------------------------------------------------------------------------
# Define compiler flags
#

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++ -ggdb")
  if(APPLE) 
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunknown-pragmas")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

get_property(GEANT5_COMPILE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
get_property(GEANT5_COMPILE_FLAGS_NOW GLOBAL PROPERTY COMPILE_FLAGS)
set(GEANT5_CXX_NOW_FLAGS                "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG               "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-2 -std=c++11")
set(GEANT5_CXX_NOW_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
set(GEANT5_CXX_NOW_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} -std=c++11")
set(GEANT5_CXX_NOW_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
set(GEANT5_CXX_NOW_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -std=c++11")
set(GEANT5_EXE_NOW_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} -std=c++11")
if (CUDA)
   # Currently we are stuck with gcc 4.8 and a broken Weffc++ so let's not us it there
   set(GEANT5_WARNING_FLAGS "-W -Wall -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")
else()
   set(GEANT5_WARNING_FLAGS "-W -Wall -Weffc++ -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")
endif()

set(GEANT5_CXX_FLAGS                "${CMAKE_CXX_FLAGS} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ${GEANT5_WARNING_FLAGS}")
set(GEANT5_COMPILE_FLAGS            "${GEANT5_COMPILE_FLAGS_NOW} ${GEANT5_WARNING_FLAGS}")

#---------------------------------------------------------------------------
if(CUSTOM_ALLOCATOR)
  find_package(CAllocator)
endif()
if (CUDA)
  find_package(CUDA REQUIRED)
  include(FindCUDA)

  if (NOT CUDA_ARCH)
    if (CUDA_SEPARABLE_COMPILATION)
      # Separable compilation will 'succeed' with sm50 and sm21 but sm50 is
      # not supported by the K20 and sm21 is too old for our code ....
      # The related bug is fixed in nvcc.  This fix will be part of nvcc
      # released after June 2014 (so likely 6.0.2)
      set(CUDA_ARCH 50)
    else()
      set(CUDA_ARCH 35)
    endif()
  endif()
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")

  #set(CUDA_NVCC_FLAGS --ptxas-options=-v)
  set(
    CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}
   -std=c++11
    -O2 --use_fast_math -Xcudafe "--diag_suppress=code_is_unreachable"
    -Xcudafe "--diag_suppress=initialization_not_reachable"
  )
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -g -G)
  endif()
  # set(CUDA_SEPARABLE_COMPILATION ON)
  # set(CUDA_VERBOSE_BUILD ON)

  if (VECGEOM_CUDA_LIBRARY)
     SET(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} ${CUDA_LIBRARIES} )
     message(STATUS "Updated VecGeom libraries to ${VECGEOM_LIBRARIES}")
  endif()

  include_directories(AFTER SYSTEM ${CUDA_INCLUDE_DIRS})
  message(STATUS "Adde to include directorys:  ${CUDA_INCLUDE_DIRS}")

  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")

endif()

#---------------------------------------------------------------------------
# Include directories
include_directories("${PROJECT_SOURCE_DIR}/vecprot_v2/inc")
include_directories("${PROJECT_SOURCE_DIR}/base/inc")
include_directories("${PROJECT_SOURCE_DIR}/xsec/inc")
include_directories("${PROJECT_SOURCE_DIR}/examples/inc")


# Work around VecCore not being its own project yet.
include_directories("${VECGEOM_INSTALL_DIR}/include/VecCore")
include_directories("${VECGEOM_INCLUDE_DIR}/include/VecCore")

message(STATUS "VECGEOM_INSTALL_DIR is -${VECGEOM_INSTALL_DIR}-")
message(STATUS "VECGEOM_INCLUDE_DIR is -${VECGEOM_INCLUDE_DIR}-")

#----------------------------------------------------------------------------
# Add subdirectories
#
add_custom_target(_headers ALL)

add_subdirectory (base)
add_subdirectory (xsec)
add_subdirectory (vecprot_v2)
add_subdirectory (examples)
add_subdirectory (tabxsec)
add_subdirectory (micromc)
#add_subdirectory (test-small)
add_subdirectory (cmstrack)
add_subdirectory (test-complex)

if (CUDA)
  add_subdirectory (vecprot_v2/cuda)
endif()

if(USE_VECPHYS)
 add_subdirectory (vectphysproc)
endif()

if(CTEST)
  add_subdirectory(vecprot_v2/tests)
endif()
