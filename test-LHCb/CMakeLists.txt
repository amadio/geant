# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(include 
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Define compiler flags : already done in the main CMakeLists
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

file(GLOB sources src/*.cc)
file(GLOB headers include/*.hh)

set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(testlhcb testLHCb.cc ${sources} ${headers})
set_target_properties(tabxsec PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(testlhcb ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}
  Geant_v Xsec)

set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")

add_dependencies(testlhcb Geant_v Xsec)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TESTLHCB. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(TESTLHCB_SCRIPTS
    g4macro_TPHYS.mac g4macro_G4.mac testLHCb.sh 
   )

foreach(_script ${TESTLHCB_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
#
add_custom_target(TESTLHCB_T DEPENDS testlhcb)

#----------------------------------------------------------------------------
# Setting Geant4 version compiler definition  
#
if (Geant4_VERSION_MAJOR LESS 10)
  set_target_properties(testlhcb PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=1") 
else()
  set_target_properties(testlhcb PROPERTIES COMPILE_DEFINITIONS "GEANT4_LESS_10=0") 
endif()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS testlhcb DESTINATION bin)

