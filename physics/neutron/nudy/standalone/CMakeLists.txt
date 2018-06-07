set(NUDYSRC "${CMAKE_SOURCE_DIR}/physics/neutron/nudy")

file (GLOB sources 
 ${NUDYSRC}/src/*.cc
 ${NUDYSRC}/standalone/NudyMain.cc
)

file (GLOB headers ${NUDYSRC}/inc/Geant/*.h)

include_directories( ${NUDYSRC}/inc)

#---------------------------------
#  Executable
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/nudy/bin)


add_executable(NudyMain ${NUDYSRC}/standalone/NudyMain.cc ${sources})
target_link_libraries(NudyMain -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

#install(TARGETS NudyMain DESTINATION bin)