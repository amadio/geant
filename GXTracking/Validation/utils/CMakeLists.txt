################
# cmake options
################
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT(utils)

##.. inherits CMAKE_MODULE_PATH and ROOT_* variables from top CMakeLists
#message(STATUS "From utils: CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
#message(STATUS "From utils: ROOT_INCLUDE_DIR = " ${ROOT_INCLUDE_DIR})
#message(STATUS "From utils: GPUNONRANDOM = " ${GPUNONRANDOM})
#message(STATUS "From utils: GPUDEBUG = " ${GPUDEBUG})
#message(STATUS "From utils: GPUPLOTS = " ${GPUPLOTS})

## Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_CMAKE_DIR   cmake   CACHE PATH "Installation directory for CMake files")
#
##.. Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()


##.. these are interesting to automatically build a setup script using CMake
# configure a header file to pass some of the CMake settings to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/utilsConfig.hh.in"
#  "${PROJECT_BINARY_DIR}/utilsConfig.hh" @ONLY
#  )

#.. Add the binary tree to the search path for include files, so that we will find ValidationConfig.h
include_directories(
   ${PROJECT_SOURCE_DIR})   # to find GPHistoManager.hh
#   ${PROJECT_BINARY_DIR})  # to find ValidationConfig.hh

################
### Adding ROOT libraries - see http://www.cmake.org/cmake/hep/cmake_tutorial.html#s2
################

include_directories ("${ROOT_INCLUDE_DIR}")

add_library(utils src/GPHistoManager.cc)

##set_target_properties(Validation PROPERTIES PUBLIC_HEADER "include/GPHistoManager.hh;${PROJECT_BINARY_DIR}/ValidationConfig.hh")
#set_target_properties(utils PROPERTIES PUBLIC_HEADER "include/GPHistoManager.hh")

install(TARGETS utils
  # IMPORTANT: Add the utils library to the "export-set"
  EXPORT utilsTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin   # for binaries
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # for shared lib
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib   # for static lib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/utils" COMPONENT include)
#  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/utils")

############
# Automation: taken from http://www.cmake.org/Wiki/CMake/Tutorials/How_to_create_a_ProjectConfig.cmake_file
############

# Add all targets to the build-tree export set
#export(TARGETS Validation  FILE "${PROJECT_BINARY_DIR}/ValidationTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE utils)
 
# Create the utilsConfig.cmake and utilsConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")

#configure_file(UtilsConfig.cmake.in "${PROJECT_BINARY_DIR}/UtilsConfig.cmake" @ONLY)

# ... for the install tree
set(CONF_INCLUDE_DIRS "\${UTILS_CMAKE_DIR}/${REL_INCLUDE_DIR}")
#configure_file(UtilsConfig.cmake.in
#  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UtilsConfig.cmake" @ONLY)

# ... for both
#configure_file(UtilsConfigVersion.cmake.in
#  "${PROJECT_BINARY_DIR}/UtilsConfigVersion.cmake" @ONLY)
 
# Install the UtilsConfig.cmake and UtilsConfigVersion.cmake
#install(FILES
#  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UtilsConfig.cmake"
#  "${PROJECT_BINARY_DIR}/UtilsConfigVersion.cmake"
#  DESTINATION "${INSTALL_CMAKE_DIR}")
 
# Install the export set for use with the install-tree
install(EXPORT utilsTargets DESTINATION "${INSTALL_CMAKE_DIR}")
