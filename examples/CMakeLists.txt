#----------------------------------------------------------------------------------------------
# Template/Macro to use in examples' CMakeLists.txt
#
macro(GeantExampleDef)
   GeantRelPath()

   include_directories( inc )

   if (USE_TBB)
   include_directories(${TBB_INCLUDE_DIR})
   endif()

   # This assumes that the files are examples/ExampleName/Geant[4V]

   get_filename_component(Category ${RelativeCurrentSourceDir} NAME)
   set(IsGeantV "${Category}" STREQUAL "GeantV")

   if (NOT DEFINED ExampleName) 
   if (DEFINED ExampleMainSource)
      get_filename_component(ExampleName ${ExampleMainSource} NAME_WE)
   else()
      get_filename_component(ExampleName ${RelativeCurrentSourceParent} NAME)
   endif()
   endif()
   if (NOT DEFINED ExampleMainSource)
   set(ExampleMainSource ${ExampleName}.cc)
   endif()

   #----------------------------------------------------------------------------------------------
   # Build executables
   #
   if (NOT DEFINED ExampleOutputDir) 
   string(REGEX REPLACE
         "^examples/([^/]*)"
         "bin/examples/\\1"
         ExampleOutputDir
         ${RelativeCurrentSourceDir})
   endif()

   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${ExampleOutputDir})


   add_executable(${ExampleName} ${ExampleMainSource})
   target_link_libraries(${ExampleName} -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v RealPhysics GeantExamplesRP ${ROOT_LIBRARIES} -lpthread ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${NUMA_ALL_LIBRARIES})
   add_dependencies(${ExampleName} Geant_v RealPhysics GeantExamplesRP GeantNuma)


   if (CTEST AND ${ExampleCTestLevel} LESS_EQUAL ${CTEST_LEVEL})
   add_test(NAME ${ExampleName}
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} 
            COMMAND $<TARGET_FILE:${ExampleName}>)
   set_tests_properties(${ExampleName} PROPERTIES ENVIRONMENT GEANT_PHYSICS_DATA=${CMAKE_SOURCE_DIR}/physics/data)
   endif()

   link_directories(${PROJECT_BINARY_DIR}/lib)


   #----------------------------------------------------------------------------
   # Copy all scripts and install the executable(s)
   #
   foreach(_script ${ExampleScripts})
   get_filename_component(script_name ${_script} NAME)
   # message(WARNING "conf ${_script} to ${CMAKE_BINARY_DIR}/${FullCMS_OUTPUT}/${script_name} ")
   configure_file(
      ${_script}
      ${CMAKE_BINARY_DIR}/${ExampleOutputDir}/${script_name}
      COPYONLY
   )
   endforeach()
   install(TARGETS ${ExampleName} DESTINATION ${ExampleOutputDir})
   install(FILES ${ExampleScripts} DESTINATION ${ExampleOutputDir})

   #----------------------------------------------------------------------------------------------
   # Add this directory to Doxygen
   #
   if(DOXYGEN_FOUND)
   set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
   endif()

endmacro(GeantExampleDef)

#----------------------------------------------------------------------------------------------
# Add subdirectories
#

add_subdirectory (physics)
add_subdirectory (magneticfield/simplifiedCMS)
