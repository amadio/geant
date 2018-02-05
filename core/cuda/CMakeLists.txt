#----------------------------------------------------------------------------------------------
# Add source files & include directories
#

# All currently handled by run/scheduler/cuda/CMakeLists.txt

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
