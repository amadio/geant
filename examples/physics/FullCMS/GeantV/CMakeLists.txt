#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/physics/generators/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
)

include_directories(
  inc
)

if (USE_TBB)
  include_directories(${TBB_INCLUDE_DIR})
endif()

#----------------------------------------------------------------------------------------------
# Build executables
#

set(NOTEST_EXECUTABLES
   FullCMS.cc
)

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    target_link_libraries(${TARGET_NAME} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics GeantExamplesRP ${ROOT_LIBRARIES} -lpthread ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
    add_dependencies(${TARGET_NAME} Geant_v RealPhysics GeantExamplesRP GeantNuma)
    set(GEANTV_EXECUTABLES ${EXECUTABLES} ${TARGET_NAME})
  endforeach()
endmacro()

function(add_to_ctest EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    if (CTEST)
      get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
      if(TARGET_NAME STREQUAL "TestEm3_GV")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      elseif(TARGET_NAME STREQUAL "TestEm5_GV")
          add_test(NAME ${TARGET_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND  $<TARGET_FILE:${TARGET_NAME}>)
      else()
        add_test(NAME ${TARGET_NAME} ${TARGET_NAME})
      endif()
      set_tests_properties(${TARGET_NAME} PROPERTIES ENVIRONMENT GEANT_PHYSICS_DATA=${CMAKE_SOURCE_DIR}/physics/data)
    endif()
  endforeach()
endfunction()

link_directories(${PROJECT_BINARY_DIR}/lib)

build_executables("${TEST_EXECUTABLES}")
build_executables("${NOTEST_EXECUTABLES}")
if(CTEST)
  add_to_ctest("${TEST_EXECUTABLES}")
endif()


#----------------------------------------------------------------------------
# Copy all scripts and install the executable(s)
#
### FullCMS ###
set(FullCMS_OUTPUT bin/examples/FullCMS/GeantV)
set(FullCMS_SCRIPTS FullCMS.mac README ${CMAKE_SOURCE_DIR}/examples/physics/FullCMS/Geant4/cms.gdml)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${FullCMS_OUTPUT})
foreach(_script ${FullCMS_SCRIPTS})
  get_filename_component(script_name ${_script} NAME)
  # message(WARNING "conf ${_script} to ${CMAKE_BINARY_DIR}/${FullCMS_OUTPUT}/${script_name} ")
  configure_file(
    ${_script}
    ${CMAKE_BINARY_DIR}/${FullCMS_OUTPUT}/${script_name}
    COPYONLY
  )
endforeach()
install(TARGETS FullCMS DESTINATION ${FullCMS_OUTPUT})
install(FILES ${FullCMS_SCRIPTS} DESTINATION ${FullCMS_OUTPUT})

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
