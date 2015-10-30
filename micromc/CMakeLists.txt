#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/*.h)



#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
include_directories(inc ${PROJECT_SOURCE_DIR}/vecprot_v2/inc)
include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

add_executable(micromc src/micromc.cxx ${sources} ${headers})
target_link_libraries(micromc ${ROOT_LIBRARIES} -lGeom -lEG -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${XSEC_LINK_LIBS} -lGeant_v)
target_link_libraries(micromc ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})

add_dependencies(micromc xseclib Geantlib_v)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(MuMC DEPENDS micromc)
install(TARGETS micromc DESTINATION bin)
