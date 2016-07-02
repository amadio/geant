# Adding this directory to Doxygen
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()


#------------------------------------------------------------------------------
# - CMakeLists.txt for CUDA coprocessor broker
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.8)

# Steer the name of the library (main library replacement or plug-in).

set(GeantCudaLib Geant_v) # or GeantCuda
set(GeantCudaReplace ON)

set(CUDA_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CUDA_INCDIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(GEANT_COMPILATION_FLAGS ${GEANT_COMPILATION_FLAGS} -DGEANT_CUDA)
add_definitions(-DGEANT_CUDA)

if (GEANT_CUDA_SINGLE_SOURCE)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA_SEP_COMP)
endif()

INCLUDE_DIRECTORIES(${CUDA_INCDIR} ${CMAKE_SOURCE_DIR}/xsec/inc ${CMAKE_SOURCE_DIR}/magneticfield/inc )

set(dict_headers ${CUDA_INCDIR}/CoprocessorBroker.h)

set(cuda_cxx_sources
    ${CUDA_SRCDIR}/GeantCudaUtils.cxx
    ${CUDA_SRCDIR}/CoprocessorBroker.cxx
    G__GeantCuda.cxx
)


if (GEANT_CUDA_SINGLE_SOURCE)
set(cuda_sources
    ${CUDA_SRCDIR}/SingleSource.cu
 )
else()
message(WARNING "Use GEANT_CUDA_SINGLE_SOURCE=OFF is not working properly yet as the object fiels for the cxx files are not passed to the device link phase")
set(cuda_sources
   ${CUDA_SRCDIR}/../../src/GeantTrack.cxx
   ${CUDA_SRCDIR}/../../src/GeantTaskData.cxx
   ${CUDA_SRCDIR}/../../src/ScalarNavInterfaceVGM.cxx
   ${CUDA_SRCDIR}/GeantTrack.cu
   ${CUDA_SRCDIR}/PropagateTracks.cu
   ${CUDA_SRCDIR}/GeantCudaUtils.cu
)
endif()

if (GeantCudaReplace)
   set (dict_headers ${dict_headers}
         ${GEANTV_HEADERS}
         # ${XSEC_HEADERS}
   )
  set (cuda_cxx_sources ${cuda_cxx_sources}
         ${GEANTV_SOURCES}
         # ${XSEC_SOURCES}
   )
   # Remove source that are (indirectly) compiled by nvcc.
#   list(REMOVE_ITEM  cuda_cxx_sources ${CMAKE_HOME_DIRECTORY}/vecprot_v2/src/GeantTrack.cxx)
#   list(REMOVE_ITEM  cuda_cxx_sources ${CMAKE_HOME_DIRECTORY}/vecprot_v2/src/GeantTaskData.cxx)

   set (extra_dict_options -DGEANTCUDA_REPLACE)
endif()

ROOT_GENERATE_DICTIONARY(G__GeantCuda ${dict_headers}
                MODULE ${GeantCudaLib}
                LINKDEF ${CUDA_INCDIR}/CudaLinkDef.h
                OPTIONS -I${CUDA_INCLUDE_DIRS} -DVC_COMMON_SIMDIZE_H_=1 ${extra_dict_options} )

# add the command to generate the source code
#add_custom_command (
#  OUTPUT ${CUDA_DATA_FILE}
# COMMAND (cd data; tar xfz $(notdir $(CUDA_DATA_TARFILE)) ) && touch $@
#  DEPENDS ${CUDA_DATA_TARFILE}
#)

#add_custom_command (
# OUTPUT ${CUDA_DATA_TARFILE}
#  COMMAND mkdir -p data && wet https://oink.fnal.gov/perfanalysis/vp${OUTPUT} -O ${OUTPUT}
##  DEPENDS
#)


CUDA_ADD_LIBRARY( ${GeantCudaLib}
                  ${cuda_sources} ${cuda_cxx_sources}
                 SHARED
                 OPTIONS ${CUDA_ARCH}
)

cuda_add_library_depend(${GeantCudaLib} ${VECGEOM_CUDA_STATIC_LIBRARY})
target_link_libraries(${GeantCudaLib} ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL})

if( VECGEOM_FOUND )
  target_link_libraries(${GeantCudaLib} Vmagfield ${ROOT_LIBRARIES} -lEG -lGeom -lGui  ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
else()
  target_link_libraries(${GeantCudaLib} Vmagfield ${ROOT_LIBRARIES} -lEG -lGeom -lGui )
endif()

add_custom_target(Geantlib_v DEPENDS ${GeantCudaLib})

if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${GeantCudaLib}${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS ${GeantCudaLib} DESTINATION lib)
endif()

