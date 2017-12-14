#----------------------------------------------------------------------------------------------
# Add source files & include directories
#

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

include_directories(include
  ${Geant4_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/physics/particle/inc)

if (USE_VECGEOM_NAVIGATOR)
 include_directories(${CMAKE_SOURCE_DIR}/material/inc)
endif()

if(USE_NUMA)
  include_directories(${CMAKE_SOURCE_DIR}/numa/inc)
endif()

file(GLOB sources src/*.cc)
file(GLOB headers include/*.hh)

#----------------------------------------------------------------------------------------------
# Redefine compiler flags since Geant4 is overriding them
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
add_executable(testcomplex testcomplex.cc ${sources} ${headers})
set_target_properties(testcomplex PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(testcomplex -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}
  ${XSEC_LINK_LIBS} Geant_v Xsec)

set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")

add_dependencies(testcomplex Geant_v Xsec)

#----------------------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TESTCOMPLEX. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(TESTCOMPLEX_SCRIPTS
  g4macro_TPHYS.mac g4macro_G4.mac testcomplex.sh
  )

foreach(_script ${TESTCOMPLEX_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
#
add_custom_target(TESTCOMPLEX_T DEPENDS testcomplex)

#----------------------------------------------------------------------------
# Setting Geant4 version compiler definition
#
if (Geant4_VERSION_MAJOR LESS 10)
  set_target_properties(testcomplex PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=1")
else()
  set_target_properties(testcomplex PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=0")
endif()


#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testcomplex DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
