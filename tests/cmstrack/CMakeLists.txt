#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include(${Geant4_USE_FILE})

if( SNEXTG4 )
  add_definitions("-DSNEXTG4")
endif()

if ( NAVTEST ) 
  add_definitions("-DNAVTEST")
endif()

find_package(Pythia8 REQUIRED)

include_directories(inc)
include_directories(AFTER SYSTEM ${Geant4_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR})
file(GLOB sources src/*.cc)
file(GLOB headers inc/*.hh)

#----------------------------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(load_gdml load_gdml.cc ${sources} ${headers})
target_link_libraries(load_gdml -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} 
  ${PYTHIA8_LIBRARIES})

#----------------------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G01_SCRIPTS
  g01.in event.in vis.mac cms.gdml gungeantino.in materials.xml scanevent.C vis.mac)

foreach(_script ${G01_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_CURRENT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------------------------
# Add program to the project targets
#
add_custom_target(G01 DEPENDS load_gdml)

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS load_gdml DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

