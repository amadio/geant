#----------------------------------------------------------------------------------------------
# The GeantV version of the real-physics tests
#
add_subdirectory(GeantV/TestEm0)
add_subdirectory(GeantV/bremTest)
add_subdirectory(GeantV/comptonTest)
add_subdirectory(GeantV/xsectionsTest)

#----------------------------------------------------------------------------------------------
# The corresponding Geant4 version of the tests
#
add_subdirectory(Geant4/TestEm0)
add_subdirectory(Geant4/bremTest)
add_subdirectory(Geant4/comptonTest)
