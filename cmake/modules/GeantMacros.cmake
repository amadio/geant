 # define a function downloading from a URL into a local file LOCALFILE
 function(FILE_DOWNLOAD FILE_URL LOCALFILE )
  if(APPLE)
      #execute_process(COMMAND curl ${FILE_URL} -o  ${LOCALFILE})
      file(DOWNLOAD ${FILE_URL} ${LOCALFILE})
  else()
     execute_process(COMMAND wget -q ${FILE_URL} -O ${LOCALFILE})
     #file(DOWNLOAD ${FILE_URL} ${FILE_URL})
  endif()
 endfunction(FILE_DOWNLOAD)
 # end of function FILE DOWNLOAD

 # define a function checking md5 hashes
 # result is stored in MD5MATCHES ( 1 == true, 0 == false )
 function(CHECKMD5 FILETOCHECK EXPECTEDMD5HASH MD5MATCHES)
     if(APPLE)
         execute_process(COMMAND md5 ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
         string(LENGTH ${MD5SUM} MD5LENGTH)
         MATH(EXPR START "${MD5LENGTH} - 33")
         string(SUBSTRING ${MD5SUM} ${START} 32 MD5SUM)
     else()
         execute_process(COMMAND md5sum ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
         string(SUBSTRING ${MD5SUM} 0 32 MD5SUM)
     endif()
     if(MD5SUM STREQUAL EXPECTEDMD5HASH)
       set(${MD5MATCHES} 1 PARENT_SCOPE)
     else()
       set(${MD5MATCHES} 0 PARENT_SCOPE)
     endif()
 endfunction(CHECKMD5)

 # actual function for managing the download
 function(DOWNLOAD_IF_NOT_INSTALLED FILE_URL LOCALFILE TARGETPATH MD5HASH )
   find_file(FOUNDFILE NAMES ${LOCALFILE} PATHS ${TARGETPATH})
   if(FOUNDFILE STREQUAL "FOUNDFILE-NOTFOUND")
       # set need download
       message(STATUS "need download of ${LOCALFILE} since not found")
       set( NEEDTODOWNLOAD 1 )
   else()
       # check md5
       message(STATUS "found existing file ${LOCALFILE}")
       CHECKMD5( ${FOUNDFILE} ${MD5HASH} MD5CORRECT )
       if( ${MD5CORRECT} STREQUAL "1" )
           # do not set download flag
           set( NEEDTODOWNLOAD 0 )
       else( )
           # set need download
           message(STATUS "hash ${MD5HASH} not correct for file ${FOUNDFILE} ${MD5CORRECT}" )
           set( NEEDTODOWNLOAD 1 )
       endif( )
   endif()

   if( ${NEEDTODOWNLOAD} STREQUAL 1 )
       message(STATUS " downloading ... ")
       set(DOWNLOADLOCATION "${TARGETPATH}/${LOCALFILE}")
       FILE_DOWNLOAD( ${FILE_URL} ${DOWNLOADLOCATION} )
   else()
       message(STATUS " doing nothing ... ")
   endif()
   # in principle have to check now if download succeeded and has right MD5
   # TOBEDONE

   # this is annoying but we have to clear FOUNDFILE SINCE THIS IS TREATED LIKE A STATIC VARIABLE
   unset(FOUNDFILE CACHE)
 endfunction(DOWNLOAD_IF_NOT_INSTALLED)

# function to calculate the path of the current source directory related to the project source directory
# The result is stored in ${RelativeCurrentSourceDir}
function(GeantRelPath)
  get_filename_component(RealSourceDir ${CMAKE_SOURCE_DIR} REALPATH)
  get_filename_component(RealCurDir ${CMAKE_CURRENT_SOURCE_DIR} REALPATH)
  string(REGEX REPLACE "${RealSourceDir}[/\\]" ""  RelativeCurrentSourceDir ${RealCurDir})
  set(RelativeCurrentSourceDir ${RelativeCurrentSourceDir} PARENT_SCOPE)
  get_filename_component(RelativeCurrentSourceParent ${RelativeCurrentSourceDir} DIRECTORY)
  set(RelativeCurrentSourceParent  ${RelativeCurrentSourceParent} PARENT_SCOPE)
endfunction(GeantRelPath)

# Generic Physics Test CMakeLists.txt
# If the source file name is not specified, assumes that it is the same name as the directory name.
function(GeantPhysicsTest)

CMAKE_PARSE_ARGUMENTS(GeantPhysicsTest "" "OPTIONAL MAIN" "SCRIPTS;INCDIRS" ${ARGN})

GeantRelPath()

get_filename_component(Category ${RelativeCurrentSourceDir} NAME)

set(IsGeantV "${Category}" STREQUAL "GeantV")

if (NOT DEFINED GeantPhysicsTest_MAIN)
  get_filename_component(GeantPhysicsTest_MAIN ${RelativeCurrentSourceParent} NAME)
  get_filename_component(Category ${RelativeCurrentSourceDir} NAME)
  if( (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GeantPhysicsTest_MAIN}.cc) AND ${IsGeantV})
    set(GeantPhysicsTest_MAIN "${GeantPhysicsTest_MAIN}_GV")
  endif()
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GeantPhysicsTest_MAIN}.cc)
     message(FATAL_ERROR "Test filename can not be guessed from the directory.
     Either create the file ${RelativeCurrentSourceDir}/${GeantPhysicsTest_MAIN}.cc
     or pass the name via the MAIN parameter to GeantPhysicsTest (without the extension)" )
  endif()
endif()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/README)
  set(GeantPhysicsTest_SCRIPTS "README;${GeantPhysicsTest_SCRIPTS}")
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cc)
file(GLOB headers inc/*.h)

include_directories( inc )

if (${IsGeantV})
  include_directories(
    ${VECGEOM_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/core/numa/inc
    ${CMAKE_SOURCE_DIR}/core/interfaces/inc
    ${CMAKE_SOURCE_DIR}/core/concurrency/inc
    ${CMAKE_SOURCE_DIR}/run/scheduler/inc
    ${CMAKE_SOURCE_DIR}/core/base/inc
    ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
    ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
    ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
    ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
    ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
    ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
    ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
  )
else()
  # Geant4
  include_directories(${Geant4_INCLUDE_DIR})
endif()

foreach(_incdir ${GeantPhysicsTest_INCDIRS})
  include_directories(${CMAKE_SOURCE_DIR}/${_incdir})
endforeach()

#----------------------------------------------------------------------------------------------
# Executable
#

string(REGEX REPLACE
       "([^/]*)/tests"
       "bin/tests/\\1"
       TestOutputDir
       ${RelativeCurrentSourceDir})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TestOutputDir})
add_executable(${GeantPhysicsTest_MAIN} ${GeantPhysicsTest_MAIN}.cc ${sources})
target_link_libraries(${GeantPhysicsTest_MAIN} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material RealPhysics ${VECGEOM_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build/install directory.
#

foreach(_script ${GeantPhysicsTest_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${GeantPhysicsTest_MAIN} DESTINATION ${TestOutputDir})
install(FILES ${GeantPhysicsTest_SCRIPTS} DESTINATION ${TestOutputDir})

endfunction(GeantPhysicsTest)