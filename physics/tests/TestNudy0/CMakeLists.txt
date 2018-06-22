#----------------------------------------------------------------------------------------------
# The GeantV version of the physics test
#
add_subdirectory(GeantV)

#----------------------------------------------------------------------------------------------
# The corresponding Geant4 version of the test (only if build with Geant4 was requested)
#
if (WITH_GEANT4 AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Geant4)
  add_subdirectory(Geant4)
endif()
