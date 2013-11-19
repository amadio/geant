# - Finds CUSTOM_ALLOCATOR
# This module sets up ROOT information 
# It defines:
# CALLOC_FOUND        If custom allocator is found
# CALLOC_LIBRARY      Location of the custo allocator library

find_library(CALLOC_LIBRARY ${CUSTOM_ALLOCATOR} PATHS ${CUSTOM_ALLOCATOR_DIR} /usr/local/lib)

if(NOT CALLOC_LIBRARY) 
  set(CALLOC_FOUND FALSE)
  message(STATUS "Custom Allocator ${CUSTOM_ALLOCATOR} not found")
else()    
  set(CALLOC_FOUND TRUE)
  message(STATUS "Found ${CALLOC_LIBRARY} Custom Allocator")

endif()

  
