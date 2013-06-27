#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
include_directories(inc ${PROJECT_SOURCE_DIR}/vecprot/inc)
include_directories(${PROJECT_SOURCE_DIR}/xsec/inc)

add_executable(micromc src/micromc.cxx ${sources} ${headers})
target_link_libraries(micromc ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lGeom -lEG -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec -lGeant)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(MuMC DEPENDS micromc)
