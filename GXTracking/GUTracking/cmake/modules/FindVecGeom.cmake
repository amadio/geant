# Finds VecGeom installation ( the vectorized geometry library )
# This module sets up VecGeom information 
# It defines:
# VECGEOM_FOUND          If the library is found
# VECGEOM_INCLUDE_DIR    PATH to the include directory
# VECGEOM_LIBRARIES      Most common libraries
# VECGEOM_LIBRARY_DIR    PATH to the library directory 

# look if an environment variable VCROOT exists

set(VECGEOMROOT $ENV{VECGEOMROOT})

find_library(VECGEOM_LIBRARIES libvecgeom.a PATHS ${VECGEOMROOT}/lib)
# find_library(USOLIDS_LIBRARIES libusolids.a PATHS ${VECGEOMROOT}/lib)
# if (VECGEOM_LIBRARIES AND USOLIDS_LIBRARIES) 
if (VECGEOM_LIBRARIES) 
   set(VECGEOM_FOUND TRUE)	
   string(REPLACE "/lib/libvecgeom.a" "" VECGEOMROOT  ${VECGEOM_LIBRARIES})
#   set(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} ${USOLIDS_LIBRARIES})

   set(VECGEOM_INCLUDE_DIR ${VECGEOMROOT}/include/VecGeom)
   set(VECGEOM_LIBRARY_DIR ${VECGEOMROOT}/lib)

   if (CUDA)
      find_library(VECGEOM_CUDA_LIBRARY libvecgeomcuda.a PATHS ${VECGEOMROOT}/lib)
      if (VECGEOM_CUDA_LIBRARY)
         SET(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} ${VECGEOM_CUDA_LIBRARY} )
      endif()
   endif()
   message(STATUS "Found VecGeom in ${VECGEOMROOT} with ${VECGEOM_LIBRARIES}")
   message(STATUS "VECGEOM_INCLUDE_DIR = ${VECGEOM_INCLUDE_DIR}")
   message(STATUS "VECGEOM_LIBRARY_DIR = ${VECGEOM_LIBRARY_DIR}")
else()
   message(STATUS "VecGeom library not found; try to set a VECGEOMROOT environment variable to the base   installation path or add -DVECGEOMROOT = to the cmake command")	
endif()



