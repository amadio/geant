# - Finds TBB installation 
# This module sets up TBB information 
# It defines:
# TBB_FOUND          
# TBB_INCLUDE_DIR    PATH to the include directory
# TBB_LIBRARIES      Most common libraries
# TBB_LIBRARY_DIR    PATH to the library directory 

# look if an environment variable TBBROOT exists
set(TBBROOT $ENV{TBBROOT})

if( NOT TBBROOT )
       	set(TBB_FOUND FALSE)	
	message(STATUS "TBB not found; try to set the TBBROOT environment variable to the base installation path or add -DTBBROOT= to the cmake command")	
else()	
	find_library(TBB_LIBRARIES libtbb.so PATHS ${TBBROOT}/lib )
	if( TBB_LIBRARIES  )	
		set(TBB_FOUND TRUE)	
		set(TBB_INCLUDE_DIR ${TBBROOT}/include)
		message(STATUS "Found TBB library in ${TBB_LIBRARIES}")		
	else()
		message(STATUS "No TBB library found in ${TBBROOT}/lib")		
	endif()
endif()


