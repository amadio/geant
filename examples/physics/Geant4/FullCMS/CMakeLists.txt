#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(full_cms)

#----------------------------------------------------------------------------
# Find Geant4 package, batch mode only executable (i.e. no need ui and vis).
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(OUTPUT bin/examples/FullCMS/Geant4)

add_executable(full_cms full_cms.cc ${sources} ${headers})
target_link_libraries(full_cms ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build/OUTPUT directory. This is so that, after
# install, we can run the executable directly because it relies on these
# scripts being in the current working directory.
#
set(FULLCMS_SCRIPTS
  bench.g4
  cms.gdml
  README
  )

foreach(_script ${FULLCMS_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${OUTPUT}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin/examples/FullCMS/Geant4' directory under the
# CMAKE_INSTALL_PREFIX
#
install(TARGETS full_cms DESTINATION ${OUTPUT})
install(FILES ${FULLCMS_SCRIPTS} DESTINATION ${OUTPUT})
