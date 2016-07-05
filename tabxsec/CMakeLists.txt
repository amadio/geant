# Add this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
message(STATUS "Using Geant4 file ${Geant4_USE_FILE}")
include(${Geant4_USE_FILE})

include_directories(inc inc)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

if(CMAKE_BUILD_TYPE MATCHES None)
  set(_cf "${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
   set(_cf "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
   set(_cf "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
   set(_cf "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
elseif(CMAKE_BUILD_TYPE MATCHES MiniSizeRel)
   set(_cf "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINISIZEREL}")
endif()

file(GLOB c_sources src/*.c)

file(GLOB sources tabxsec.cc src/*.cc src/*.cxx)
file(GLOB headers inc/*.hh src/*.h)

set_property(SOURCE ${sources} PROPERTY COMPILE_FLAGS ${_cf})

set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO}")


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(tabxsec tabxsec.cc ${sources} ${c_sources} ${headers})
set_target_properties(tabxsec PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(tabxsec ${CALLOC_LIBRARY} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} 
                              ${VECGEOM_LIBRARIES} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${XSEC_LINK_LIBS} -lGeant_v)
add_dependencies(tabxsec xseclib)

set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build tabxsec. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(_macros 
macros)

set(TABXSEC_SCRIPTS
  ${_macros}/tabxsec.in
  ${_macros}/tabxsec.out
  ${_macros}/init.mac
  ${_macros}/init_vis.mac
  ${_macros}/run1.mac
  ${_macros}/run2.mac
  ${_macros}/vis.mac
  )

foreach(_script ${TABXSEC_SCRIPTS})
  configure_file(
    ${_script}
    ${CMAKE_CURRENT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ExtracXS DEPENDS tabxsec)
install(TARGETS tabxsec DESTINATION bin)

