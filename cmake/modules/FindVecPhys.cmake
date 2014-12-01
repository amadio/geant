# !!! This should be changed after we have a proper install for vecphys !!!
# Finds VecPhys ( the vectorized physics library )
# This module sets up VecPhys information 
# It defines:
# VECPHYS_FOUND          If the library is found
# VECPHYS_INCLUDE_DIR    PATH to the include directory
# VECPHYS_LIBRARIES      Most common libraries
# VECPHYS_LIBRARY_DIR    PATH to the library directory 

# !! this should be changed after vecphys will have a proper install
set(WHERE_TO_LOOK ${PROJECT_SOURCE_DIR}/GXTracking/GUTracking)

find_library(VECPHYS_LIBRARIES libvecphys.a PATHS ${WHERE_TO_LOOK}/build)
if (VECPHYS_LIBRARIES) 
   set(VECPHYS_FOUND TRUE)	
   set(VECPHYS_INCLUDE_DIR ${WHERE_TO_LOOK}/inc)
   set(VECPHYS_LIBRARY_DIR ${WHERE_TO_LOOK}/build)

   message(STATUS "Found VecPhys: ${VECPHYS_LIBRARIES}")
   message(STATUS "VECPHYS_INCLUDE_DIR = ${VECPHYS_INCLUDE_DIR}")
   message(STATUS "VECPHYS_LIBRARY_DIR = ${VECPHYS_LIBRARY_DIR}")
else()
   message(STATUS "VecPhys library not found; was looked for ${WHERE_TO_LOOK}/build")	
endif()



