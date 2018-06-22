set(NUDYSRC "${CMAKE_SOURCE_DIR}/physics/neutron/nudy")
set(NUDYDATA "${CMAKE_SOURCE_DIR}/physics/data/neutron/nudy")

file (GLOB sources 
 ${NUDYSRC}/EndfToRoot/*.cc
)

file (GLOB headers 
  ${NUDYSRC}/inc/Geant/*.h
)

include_directories( 
  ${CMAKE_SOURCE_DIR}/physics/neutron/nudy/inc
)

#---------------------------------
#  Executable
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/nudy/bin)

add_executable(EndfToRoot ${NUDYSRC}/EndfToRoot/EndfToRoot.cc)
target_link_libraries(EndfToRoot -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

add_executable(RootToSigma ${NUDYSRC}/EndfToRoot/RootToSigma.cc)
target_link_libraries(RootToSigma -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

add_executable(EndfToPointRoot ${NUDYSRC}/EndfToRoot/EndfToPointRoot.cc)
target_link_libraries(EndfToPointRoot -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

add_executable(RootToRetrive ${NUDYSRC}/EndfToRoot/RootToRetrive.cc)
target_link_libraries(RootToRetrive -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

add_executable(Xsec ${NUDYSRC}/EndfToRoot/Xsec_plot.cc)
target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

#add_executable(runjendl ${NUDYSRC}/EndfToRoot/main2.cc)
#target_link_libraries(runjendl -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy )

#install(TARGETS EndfToRoot DESTINATION ${NUDYDATA})

#install(TARGETS RootToSigma DESTINATION ${NUDYDATA})

#install(TARGETS EndfToPointRoot DESTINATION ${NUDYDATA})

#install(TARGETS RootToRetrive DESTINATION ${NUDYDATA})

#install(TARGETS Xsec DESTINATION ${NUDYDATA})

if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
