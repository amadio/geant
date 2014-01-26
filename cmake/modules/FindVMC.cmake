# - Finds VMC installation (Virtual Monte Carlo)
# This module sets up VMC information 
# It defines:
# VMC_FOUND          If the ROOT is found
# VMC_INCLUDE_DIR    PATH to the include directory
# VMC_LIBRARIES      Most common libraries
# VMC_LIBRARY_DIR    PATH to the library directory 

# look if an environment variable VMCROOT exists
set(VMCROOT $ENV{VMCROOT})

set(PLATFORM $ENV{PLATFORM})
if(NOT PLATFORM)
execute_process(COMMAND root-config --arch
                OUTPUT_VARIABLE PLATFORMCR)
string(REPLACE "\n" "" PLATFORM ${PLATFORMCR})
endif()

find_library(VMC_LIBRARIES libgeant4vmc.so PATHS /usr/local/geant4_vmc/lib/tgt_${PLATFORM} ${VMCROOT}/lib/tgt_${PLATFORM})

if (VMC_LIBRARIES) 
   set(VMC_FOUND TRUE)	
   string(REPLACE "/lib/tgt_${PLATFORM}/libgeant4vmc.so" "" VMCROOT  ${VMC_LIBRARIES})
   set(VMC_INCLUDE_DIR ${VMCROOT}/include)
   set(VMC_LIBRARY_DIR ${VMCROOT}/lib/tgt_${PLATFORM})
   message(STATUS "Found VMC library in ${VMC_LIBRARIES}")		
else()
   message(STATUS "VMC library not found; try to set a VMCROOT environment variable to the base installation path or add -DVMCROOT= to the cmake command")
endif()



