include(cmake/modules/RecordCmdLine.cmake)

# uninstall target
configure_file(
    "cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

include(cmake/modules/DoxygenSetup.cmake)

#----------------------------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(GEANTV)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ ISO Standard")
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(APPLE)
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
# Fix bug on APPLE, this is the default everywhere
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
  set(CMAKE_MACOSX_RPATH ON)
endif()

#----------------------------------------------------------------------------------------------
# Additional options
#
option(KNC "Enable compilation for KNC." OFF)

option(KNL "Enable compilation for KNL." OFF)

option(CUDA "Enable compilation for CUDA." OFF)

option(DATA_DOWNLOAD "Enable data download" OFF)

option(MAGIC_DEBUG "Debug xsec serialisation with magic number" ON)
if(MAGIC_DEBUG)
  add_definitions(-DMAGIC_DEBUG)
endif()


set(CMAKE_MACOSX_RPATH ON)

option(GEANT_INSTALL_RPATH "Enable RPATH injection for installed binaries" ON)

option(WITH_GEANT4 "Enable libraries and test using the Geant4 libraries" ON)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

#----------------------------------------------------------------------------------------------
# Use absolute RPATH when linking to avoid needing to export LD_LIBRARY_PATH
#
if (GEANT_INSTALL_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")
endif()

set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

#----------------------------------------------------------------------------------------------
# Process various options
#
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

option(VECTORIZED_GEOMETRY "Enable vectorized geometry" OFF)
if(VECTORIZED_GEOMETRY)
  add_definitions(-DVECTORIZED_GEOMETRY)
  message(STATUS "VECTORIZED_GEOMETRY ON  : Vectorized geometry ON !")
else()
  message(STATUS "VECTORIZED_GEOMETRY OFF : Vectorized geometry OFF !")
endif()

option(USE_CALLGRIND_CONTROL "Enable callgrind profiling" OFF)
if(USE_CALLGRIND_CONTROL)
  add_definitions(-DUSE_CALLGRIND_CONTROL)
  message(STATUS "USE_CALLGRIND_CONTROL ON  : Callgrind profiling enabled")
else()
  message(STATUS "USE_CALLGRIND_CONTROL OFF : Callgrind profiling disabled")
endif()

option(BUG_HUNT "Enable crosschecks for finding bugs." OFF)
option(USE_VECGEOM_NAVIGATOR "Use the navigator from VecGeom" ON)
option(USE_ROOT "Use ROOT" ON)
option(USE_NUMA "Use NUMA" ON)
option(USE_TBB "Use threading building blocks" OFF)

option(USE_REAL_PHYSICS "Use the real physics interface in the simulation" OFF)
if(USE_REAL_PHYSICS)
  message(STATUS "USE_REAL_PHYSICS ON  : Simulation with the real physics interface")
  if (NOT (USE_VECGEOM_NAVIGATOR))
    message (FATAL_ERROR "Real physics depends on VecGeom but USE_VECGEOM_NAVIGATOR=OFF")
  endif()
  add_definitions(-DUSE_REAL_PHYSICS)
else()
  message(STATUS "USE_REAL_PHYSICS OFF : Simulation with the old physics interface")
endif()

option(USE_VECGEOM_CONTAINERS "Replace std::vector, std::map with vecgeom::Vector, vecgeom::map in the material lib." ON)
if (USE_VECGEOM_NAVIGATOR AND USE_VECGEOM_CONTAINERS)
  message(STATUS "USE_VECGEOM_CONTAINERS ON : Material library will be built with vecgeom::Vector and vecgeom::map.")
  add_definitions(-DUSE_VECGEOM_CONTAINERS)
else()
  message(STATUS "USE_VECGEOM_CONTAINERS OFF : Material library will be built with std::vector and std::map.")
endif()


#---------------------------------------------------------------------------
# Debuging options
#
if(BUG_HUNT)
  add_definitions(-DBUG_HUNT)
  message(STATUS "BUG_HUNT ON  : Bug hunting mode !")
else()
  message(STATUS "BUG_HUNT OFF  : No bug hunting...")
endif()

#----------------------------------------------------------------------------------------------
# Find packages
#
find_package(HepMC REQUIRED)
if(HepMC_FOUND)
  message(STATUS "HepMC includes: ${HEPMC_INCLUDE_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
endif()

if (KNC)
  if (USE_ROOT)
    message (WARNING "Root has to be disabled in order to compile for KNC, setting USE_ROOT=OFF")
    set (USE_ROOT OFF CACHE BOOL "Use ROOT" FORCE)
  endif()
  if (WITH_GEANT4)
    message (WARNING "Disabling Geant4 in order to compile for KNC, setting WITH_GEANT4=OFF")
    set (WITH_GEANT4 OFF CACHE BOOL "Enable libraries and test using the Geant4 libraries" FORCE)
  endif()
  if (CUDA)
    message (FATAL "CUDA and KNC are incompatible flags")
  endif()
endif()

#----------------------------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
if(WITH_GEANT4)
  if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
  else()
    find_package(Geant4 REQUIRED)
  endif()
endif()

#----------------------------------------------------------------------------------------------
# Set pathes where to put the libraries, executables and headers
#
if( GEANT_OUTPUT_DIR )
    set(OUTDIR ${GEANT_OUTPUT_DIR})
else()
    set(OUTDIR $ENV{GEANT_OUTPUT_DIR})
endif()
if(NOT OUTDIR)
  set(OUTDIR ${PROJECT_BINARY_DIR})
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

#----------------------------------------------------------------------------------------------
# Find VecGeom geometry headers library
#
set(VECGEOM_DIR "" CACHE DIRECTORY "Location of VecGeom")
if(VECGEOM_DIR)
  list(APPEND CMAKE_PREFIX_PATH ${VECGEOM_DIR})
endif()
if($ENV{VECGEOM_DIR})
  list(APPEND CMAKE_PREFIX_PATH $ENV{VECGEOM_DIR})
endif()

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

#----------------------------------------------------------------------------------------------
# Find ROOT and setup the libraries
#

if(USE_ROOT)
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
   list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
   find_package(ROOT REQUIRED COMPONENTS EG Gui Geom TreePlayer)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
   include(${ROOT_USE_FILE})
   include_directories(AFTER SYSTEM ${ROOT_INCLUDE_DIRS})

   add_definitions("-DUSE_ROOT")
else()
#----------------------------------------------------------------------------------------------
# ROOT adds some of these flags. If we do not use ROOT, we have to add them by hand
#
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -pipe -std=c++11")
endif()

#----------------------------------------------------------------------------
#Find TBB and setup the libraries
#
set (TBB_LIBRARIES)
if(USE_TBB)
  find_package(TBB)
  if( TBB_FOUND )
    add_definitions("-DGEANT_TBB")
    message(STATUS "USE_TBB ON  : TBB enabled")
  else()
    message(STATUS "TBB requested but not found")
    set(USE_TBB OFF CACHE BOOL "" FORCE)
  endif()
else()
  message(STATUS "USE_TBB OFF  : TBB disabled")
endif()

#----------------------------------------------------------------------------
# Load macros
#
include(cmake/modules/GeantMacros.cmake)

#----------------------------------------------------------------------------------------------
# see if we have CCACHE ( to avoid recompilation on branch switches )
#
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message(STATUS "found ccache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)


#----------------------------------------------------------------------------------------------
# Include creation and submission to dashboard cdash.cern.ch
#
include(CTest)

if (CTEST)
  set(DATA_DOWNLOAD ON)
  message(STATUS "Testing with CTest enabled.")
  enable_testing()
endif()

#----------------------------------------------------------------------------------------------
# Set data to be dowloaded
#
option(DATA_DOWNLOAD "Set of data files for examples" OFF)

if(DATA_DOWNLOAD)
  message(STATUS "Downloading data files")

  file(MAKE_DIRECTORY data)
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/ExN03.root"
    "ExN03.root" "${CMAKE_SOURCE_DIR}/data/" "b6b0cfdd5035117171bfe1b7f1f40c3f")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/cms2015.root"
    "cms2015.root" "${CMAKE_SOURCE_DIR}/data/" "fb3ede867532b372c5e6f7138d00c07e")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/fstate_FTFP_BERT_G496p02_1mev.root"
    "fstate_FTFP_BERT_G496p02_1mev.root" "${CMAKE_SOURCE_DIR}/data/" "d7d3fc9fe2a835f9021bd3747cabf582")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/fstate_FTFP_BERT.root"
    "fstate_FTFP_BERT.root" "${CMAKE_SOURCE_DIR}/data/" "4eea5710c691749980f5f8ad02a42231")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/xsec_FTFP_BERT_G496p02_1mev.root"
    "xsec_FTFP_BERT_G496p02_1mev.root" "${CMAKE_SOURCE_DIR}/data/" "648b4f4136d9a4cb6b590d98dc9182d8")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/xsec_FTFP_BERT.root"
    "xsec_FTFP_BERT.root" "${CMAKE_SOURCE_DIR}/data/" "90c8ef05770cea48c7807e0241921cab")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/pp14TeVminbias.root"
    "pp14TeVminbias.root" "${CMAKE_SOURCE_DIR}/data/" "109b827c595ae4a62b74f08f860abbb7")
  DOWNLOAD_IF_NOT_INSTALLED("http://geant-server.web.cern.ch/geant-server/ATLASbar.root"
    "ATLASbar.root" "${CMAKE_SOURCE_DIR}/data/" "b0a579a0158b2bf704ace490744d858d")
endif()

#----------------------------------------------------------------------------------------------
# Handle compiler flags
#
if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)
  string(REPLACE "." ";" GCC_VERSION_LIST ${GCC_VERSION})
  set(GCC_VERSION_MAJOR 0)
  set(GCC_VERSION_MINOR 0)
  set(GCC_VERSION_PATCH 0)
  if(GCC_VERSION_LENGTH GREATER 0)
    list(GET GCC_VERSION_LIST 0 GCC_VERSION_MAJOR)
  endif()
  if(GCC_VERSION_LENGTH GREATER 1)
    list(GET GCC_VERSION_LIST 1 GCC_VERSION_MINOR)
  endif()
  if(GCC_VERSION_LENGTH GREATER 2)
    list(GET GCC_VERSION_LIST 2 GCC_VERSION_PATCH)
  endif()

  if(GCC_VERSION_MAJOR LESS 5 AND GCC_VERSION_MINOR LESS 9)
    set(DISABLE_WEFFC ON)
  endif()
  if(GCC_VERSION_MAJOR LESS 5)
    if ( Vc_FOUND )  # Actually this should really be if VecGeomConfigureToUseVC
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=6")
    endif()
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  if("${CMAKE_VERSION}" VERSION_LESS "3.6")
    # add flag by hand if CMake does not support ICC
    if("${CMAKE_CXX_STANDARD}" MATCHES "11|14")
      add_compile_options(-std=c++${CMAKE_CXX_STANDARD})
    else()
      message(FATAL_ERROR "Unsupported C++ standard requested")
    endif()
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

if (NOT DISABLE_WEFFC)
  set(WEFFC_FLAG -Weffc++)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WEFFC_FLAG} -ggdb")
  if(APPLE)
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunknown-pragmas")
   else()
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
   endif()
 else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
 endif()

 # set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

 if (KNC)
   if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
   else()
     message(FATAL_ERROR "KNC native or offload compilation requires Intel Compiler.")
   endif()
 endif()


 get_property(GEANTV_COMPILE_DEFINITIONS GLOBAL PROPERTY COMPILE_DEFINITIONS)
 get_property(GEANTV_COMPILE_FLAGS_NOW GLOBAL PROPERTY COMPILE_FLAGS)
 set(GEANTV_CXX_NOW_FLAGS                "${CMAKE_CXX_FLAGS}")
 set(CMAKE_CXX_FLAGS_DEBUG               "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-2 ")
 set(GEANTV_CXX_NOW_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ")
 set(GEANTV_CXX_NOW_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ")
 set(GEANTV_CXX_NOW_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ")
 set(GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ")
 set(GEANTV_EXE_NOW_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ")
 if (CUDA)
   # Currently we are stuck with gcc 4.8 and a broken Weffc++ so let's not us it there
   set(GEANTV_WARNING_FLAGS "-W -Wall -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")
 else()
   set(GEANTV_WARNING_FLAGS "-W -Wall ${WEFFC_FLAG} -Wno-narrowing -Wwrite-strings -Wcast-qual -Wshadow -Wno-long-long -pedantic")
 endif()

 set(GEANTV_CXX_FLAGS                "${CMAKE_CXX_FLAGS} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ${GEANTV_WARNING_FLAGS}")
 set(GEANTV_COMPILE_FLAGS            "${GEANTV_COMPILE_FLAGS_NOW} ${GEANTV_WARNING_FLAGS}")

#----------------------------------------------------------------------------------------------
# Handle custom alloactor
#
if(CUSTOM_ALLOCATOR)
  find_package(CAllocator)
endif()

#----------------------------------------------------------------------------------------------
# CUDA section
#
if (CUDA)
  option(GEANT_CUDA_SINGLE_SOURCE "Compile all the GeantV CUDA code in a single compilation unit." ON)

  # We must have those over-ride placed here as they need to be done before loading FindCUDA.cmake
  set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "Use the static version of the CUDA runtime library if available")
  set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "Compile CUDA objects with separable compilation enabled.")

  find_package(CUDA REQUIRED)
  include(FindCUDA)

  if (NOT CUDA_ARCH)
    if (CUDA_SEPARABLE_COMPILATION)
      set(CUDA_ARCH 30)
    else()
      message(FATAL_ERROR "VecGeom requires using separable compilation")
    endif()
    if (NOT VECGEOM_CUDA_ARCH)
      message(FATAL_ERROR "VecGeom not configured with CUDA/NVCC")
    endif()
    string(REGEX REPLACE "-arch=sm_([0-9]+)" "\\1" VECGEOM_CUDA_ARCH_NUM "${VECGEOM_CUDA_ARCH}")
    if (${VECGEOM_CUDA_ARCH_NUM} LESS ${CUDA_ARCH})
      message(FATAL_ERROR "GeantV must be configured with a CUDA SM (${CUDA_ARCH}) less than VecGeom's (${VECGEOM_CUDA_ARCH_NUM})" )
    elseif(${CUDA_ARCH} LESS ${VECGEOM_CUDA_ARCH_NUM})
      # all good
    elseif(${CUDA_ARCH} EQUAL ${VECGEOM_CUDA_ARCH_NUM})
      # all good
    else()
      message(WARNING "GeantV must be configured with a CUDA SM (${CUDA_ARCH}) less than VecGeom's (${VECGEOM_CUDA_ARCH_NUM}) but one of the two is invalid." )
    endif()
  endif()
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")

  set(
    CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}
    -Xcompiler ,\"-Wno-unused-function\"
    -Xcompiler ,\"-Wno-unused-parameter\"
    -O2 --use_fast_math -Xcudafe "--diag_suppress=code_is_unreachable"
    -Xcudafe "--diag_suppress=initialization_not_reachable"
    -D_MWAITXINTRIN_H_INCLUDED
    -D_FORCE_INLINES
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
  message(STATUS "Added to include directories:  ${CUDA_INCLUDE_DIRS}")

endif()

#----------------------------------------------------------------------------------------------
# Include directories
#
include_directories("${PROJECT_SOURCE_DIR}/vecprot_v2/inc")
include_directories("${PROJECT_SOURCE_DIR}/base/inc")
if(USE_REAL_PHYSICS)
  include_directories("${PROJECT_SOURCE_DIR}/examplesRP/inc")
else()
  include_directories("${PROJECT_SOURCE_DIR}/xsec/inc")
  include_directories("${PROJECT_SOURCE_DIR}/fastsim/inc")
  include_directories("${PROJECT_SOURCE_DIR}/examples/inc")
endif()
include_directories("${PROJECT_SOURCE_DIR}/vecprot_v2_tbb/inc")

string(REPLACE " -DVECGEOM_ROOT" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------------------------
# Add subdirectories
#
add_custom_target(_headers ALL)

message(STATUS "COMPILATION FLAGS ARE -${GEANTV_CXX_FLAGS} - ${GEANTV_COMPILATION_FLAGS}")
add_subdirectory (base)
add_subdirectory (numa)
if(NOT USE_REAL_PHYSICS)
  add_subdirectory (xsec)
endif()
if(USE_TBB)
  add_subdirectory (vecprot_v2_tbb)
endif()
add_subdirectory (vecprot_v2)
if(USE_REAL_PHYSICS)
  add_subdirectory (examplesRP)
else()
  add_subdirectory (examples)
endif()
if(WITH_GEANT4 AND NOT (CUDA) AND USE_ROOT AND NOT (USE_REAL_PHYSICS))
  add_subdirectory (tabxsec)
  add_subdirectory (test-small)
  add_subdirectory (test-complex)
  add_subdirectory (test-LHCb)
  add_subdirectory (cmstrack)
endif()
add_subdirectory (magneticfield)
if(USE_ROOT AND NOT (USE_REAL_PHYSICS))
  add_subdirectory (fastsim)
endif()
if (ROOT_mathmore_FOUND)
  add_subdirectory (Nudy)
endif()

if(USE_VECPHYS)
  add_subdirectory(vecphys)
endif()

if (NOT (USE_REAL_PHYSICS))
  add_subdirectory(physics)
endif()

if (USE_REAL_PHYSICS)
 add_subdirectory(realphysics)
endif()

if (USE_VECGEOM_NAVIGATOR AND NOT (CUDA))
  add_subdirectory(material)
endif()




#----------------------------------------------------------------------------------------------
# Process Doxygen file
#
if(DOXYGEN_FOUND)
  message(${DOXYFILE_SOURCE_DIRS})
  configure_file("${DOXYFILE_IN}" "${DOXYFILE}" @ONLY)
endif()
