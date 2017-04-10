#----------------------------------------------------------------------------------------------
# The GeantV version of the real-physics tests
#
add_subdirectory(GeantV/TestEm0)
add_subdirectory(GeantV/bremTest)


#----------------------------------------------------------------------------------------------
# The corresponding Geant4 version of the tests
#
add_subdirectory(Geant4/TestEm0)
add_subdirectory(Geant4/bremTest)
