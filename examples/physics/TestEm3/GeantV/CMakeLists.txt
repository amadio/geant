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
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
)

set(ExampleScripts TestEm3_GV_ATLASbar.mac TestEm3_GV_magField.mac README)
# Set to 2 or more to exclude from the default set
set(ExampleCTestLevel 1)

file(GLOB ExampleSources src/*.cc)
file(GLOB ExampleHeaders inc/*.h)

# The following can usually be guessed from the directory name.
set(ExempleName TestEm3_GV)
set(ExampleMainSource ${ExempleName}.cc)
set(ExampleOutputDir bin/examples/TestEm3/GeantV)

# Generic Part that can be used for other examples.
# The macro is declared in ./examples/CMakeLists.txt

GeantExampleDef()
