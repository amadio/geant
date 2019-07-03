#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(
  ${CMAKE_SOURCE_DIR}/core/numa/inc
  ${CMAKE_SOURCE_DIR}/core/concurrency/inc
  ${CMAKE_SOURCE_DIR}/core/interfaces/inc
  ${CMAKE_SOURCE_DIR}/run/userapp/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/physics/generators/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
  ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
)

set(ExampleScripts FullCMS.mac README ${CMAKE_SOURCE_DIR}/data/cms2018.gdml
                   ${CMAKE_SOURCE_DIR}/data/cmsmagfield2015.txt )

# Set to 2 or more to exclude from the default set
set(ExampleCTestLevel 2)

## The additional source/header files are aglomerated into a single library in
#     examples/physics/CMakeLists.txt
#  To add files to a GeantV example, currently you must edit that file.
#
# file(GLOB ExampleHeaders inc/*.h)
# file(GLOB ExampleSources src/*.cxx)

# set(ExampleHeaders CMSData.h CMSDetectorConstruction.h CMSFieldConstruction.h CMSFullApp.h CMSParticleGun.h CMSPhysicsList.h CMSmagField.h)
# set(ExampleSources CMSData.cxx CMSDetectorConstruction.cxx CMSFieldConstruction.cxx CMSFullApp.cxx CMSParticleGun.cxx CMSPhysicsList.cxx)

# The following can usually be guessed from the directory name.
set(ExampleName FullCMS)
set(ExampleMainSource ${ExampleName}.cc)
set(ExampleOutputDir bin/examples/FullCMS/GeantV)

# Generic Part that can be used for other examples.
# The macro is declared in ./examples/CMakeLists.txt

GeantExampleDef()
