#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include(${Geant4_USE_FILE})

file(GLOB c_sources src/*.c)

file(GLOB sources tabxsec.cc src/*.cc src/*.cxx)
file(GLOB headers inc/*.hh src/*.h)
include_directories(inc ${CMAKE_SOURCE_DIR}/physics/particle/inc)

if (USE_VECGEOM_NAVIGATOR)
 include_directories(${CMAKE_SOURCE_DIR}/material/inc)
endif()

if(USE_NUMA)
  include_directories(${CMAKE_SOURCE_DIR}/numa/inc)
endif()

#----------------------------------------------------------------------------------------------
# Redefine compilation flags since Geant4 is overriding them
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

set_property(SOURCE ${sources} PROPERTY COMPILE_FLAGS ${_cf})

set(CMAKE_CXX_FLAGS                "${GEANTV_CXX_NOW_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${GEANTV_CXX_NOW_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${GEANTV_CXX_NOW_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${GEANTV_CXX_NOW_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${GEANTV_CXX_NOW_FLAGS_RELWITHDEBINFO}")

#----------------------------------------------------------------------------------------------
# Add the executable
#
add_executable(tabxsec tabxsec.cc ${sources} ${c_sources} ${headers})
set_target_properties(tabxsec PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(tabxsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CALLOC_LIBRARY} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}
                              ${VECGEOM_LIBRARIES} ${XSEC_LINK_LIBS} Geant_v NewPhysics)
add_dependencies(tabxsec xseclib NewPhysics)
set(CMAKE_EXE_LINKER_FLAGS         "${GEANTV_EXE_NOW_LINKER_FLAGS}")
add_custom_target(ExtracXS DEPENDS tabxsec)

#----------------------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build tabxsec. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(_macros macros)

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

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS tabxsec DESTINATION bin)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
