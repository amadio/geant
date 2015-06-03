#----------------------------------------------------------------------------
# Enable the Language Fortran

enable_language(Fortran)


#----------------------------------------------------------------------------
# Locate Sources and Headers for this project
# NB : Headers are also included

file(GLOB sources src/*.cxx src/*.f)
file(GLOB headers inc/T*.h)

include_directories(inc)


#----------------------------------------------------------------------------
# Generate the Dictionary

ROOT_GENERATE_DICTIONARY(NudyDict ${headers} MODULE libNudy LINKDEF inc/LinkDef.h)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Nudy SHARED ${sources} NudyDict.cxx)

target_link_libraries(Nudy ${ROOT_LIBRARIES} Geom EG MathMore Graf)

# for some reason the installation process alters the libraries on Mac, and the program fails
if(APPLE) 
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy.dylib DESTINATION lib)
else()
  install(TARGETS Nudy DESTINATION lib)
endif()

install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libNudy_rdict.pcm DESTINATION lib)

