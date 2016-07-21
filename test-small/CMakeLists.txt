#----------------------------------------------------------------------------------------------
# Define options
#
option(MAKESTAT "Collect detailed statistics during exampleN03 run!" OFF)
if(MAKESTAT)
  add_definitions(-DMAKESTAT)
  message(STATUS "MAKESTAT ON  : Detailed statistics IS collected during exampleN03 run!")
else()
  message(STATUS "MAKESTAT OFF : Detailed statistics IS NOT collected during exampleN03 run!")
endif()

#----------------------------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
message(STATUS "Using Geant4 file ${Geant4_USE_FILE}")
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(include 
  ${Geant4_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
  ../xsec/inc 
  ../physics/particle/inc
  )

file(GLOB sources src/*.cc)
file(GLOB headers include/*.hh)

#----------------------------------------------------------------------------------------------
# Deal with compilation flags, since they are superseded by Geant4
#
if(CMAKE_BUILD_TYPE MATCHES None)
  set(_cf ${CMAKE_CXX_FLAGS})
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
elseif(CMAKE_BUILD_TYPE MATCHES MiniSizeRel)
  set(_cf ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINISIZEREL})
endif()

separate_arguments(_cf)
set_directory_properties(PROPERTIES COMPILE_OPTIONS "${_cf}")

set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO}")

#----------------------------------------------------------------------------------------------
# Add the executable
#
add_executable(exampleN03 exampleN03.cc ${sources} ${headers})
set_target_properties(tabxsec PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(exampleN03 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} 
  ${XSEC_LINK_LIBS} Geant_v)
if(USE_VECGEOM_NAVIGATOR)
   target_link_libraries(exampleN03 -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
endif()
add_dependencies(exampleN03 xseclib Geantlib_v)
set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")

#----------------------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(N03_SCRIPTS
    exampleN03.in exampleN03.large_N.in exampleN03.out exo1.mac exo2.mac gui.mac newgeom.mac rndmSeed.mac run1.mac run2.mac vis.mac icons.mac
    visTutor/exN03Tree0.mac visTutor/exN03Vis12.loop visTutor/exN03Vis12.mac visTutor/exN03Vis13.loop visTutor/exN03Vis13.mac visTutor/exN03Vis5.mac visTutor/g4Xt.xrm visTutor/heprep.mac
  )

foreach(_script ${N03_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_CURRENT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------------------------
# Add program to the project targets
#
add_custom_target(N03 DEPENDS exampleN03)

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS exampleN03 DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()



