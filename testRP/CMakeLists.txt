#----------------------------------------------------------------------------------------------
# The GeantV version of the real-physics tests
#
add_subdirectory(GeantV/TestEm0)
add_subdirectory(GeantV/bremTest)
add_subdirectory(GeantV/comptonTest)
add_subdirectory(GeantV/xsectionsTest)
add_subdirectory(GeantV/conversionTest)

#----------------------------------------------------------------------------------------------
# The corresponding Geant4 version of the tests (only if build with Geant4 was requested)
#
if (WITH_GEANT4)
  add_subdirectory(Geant4/TestEm0)
  add_subdirectory(Geant4/bremTest)
  add_subdirectory(Geant4/comptonTest)
  add_subdirectory(Geant4/conversionTest)
endif()
