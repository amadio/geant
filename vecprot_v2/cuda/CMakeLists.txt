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

set(CudaCode ${PROJECT_SOURCE_DIR})
set(CUDA_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CUDA_INCDIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(GEANT_COMPILATION_FLAGS ${GEANT_COMPILATION_FLAGS} -DGEANT_CUDA)
add_definitions(-DGEANT_CUDA)

#NVCC_FLAGS=--device-debug -m64 -Xcompiler ,\"-g\" -Xcompiler ,\"-fPIC\"  -arch=sm_20 --use_fast_math --ptxas-options=-v -DNVCC -DGX_MULTI_OBJ_FILES

if (CUDA_SEPARABLE_COMPILATION)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA_SEP_COMP -DGX_MULTI_OBJ_FILES)
endif()


CUDA_INCLUDE_DIRECTORIES(
       ${CudaCode}/GXTracking/gpu/cuda/kernel/include
       ${CudaCode}/GXTracking/gpu/cuda/kernel/src
)

#CXXFLAGS += -I$(CUDA_INC) -I$(CUDA_INCDIR)

INCLUDE_DIRECTORIES( ${CudaCode}/GXTracking/gpu/cuda/kernel/include
       ${CUDA_INCDIR} ${CudaCode}/xsec/inc ${CudaCode}/magneticfield/inc )


#set(CUDA_DATA_TARFILE ${PROJECT_SOURCE_DIR}/data/vpdata_v1.tar.gz)
#set(CUDA_DATA_FILE ${PROJECT_SOURCE_DIR}/data/Lambda.eBrem.e-.asc)

set(dict_headers ${CUDA_INCDIR}/CoprocessorBroker.h)

set(cuda_cxx_sources
    ${CUDA_SRCDIR}/GeantCudaUtils.cxx
    ${CUDA_SRCDIR}/CoprocessorBroker.cxx
    G__GeantCuda.cxx
)


if (CUDA_SEPARABLE_COMPILATION)
set(cuda_sources
   ${CUDA_SRCDIR}/../src/GeantTrack.cxx
   ${CUDA_SRCDIR}/../src/GeantTaskData.cxx
   ${CUDA_SRCDIR}/GeantTrack.cu
   ${CUDA_SRCDIR}/PropagateTracks.cu
   ${CUDA_SRCDIR}/GeantCudaUtils.cu
)
else()
set(cuda_sources
    ${CUDA_SRCDIR}/SingleSource.cu
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
                MODULE lib${GeantCudaLib}
                LINKDEF ${CUDA_INCDIR}/CudaLinkDef.h
                OPTIONS -I${CUDA_INCLUDE_DIRS} ${extra_dict_options} )

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

if( VECGEOM_FOUND )
  target_link_libraries(${GeantCudaLib} Vmagfield ${ROOT_LIBRARIES} -lEG -lGeom -lGui  ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
else()
  target_link_libraries(${GeantCudaLib} Vmagfield ${ROOT_LIBRARIES} -lEG -lGeom -lGui )
endif()

add_custom_target(Geantlib_v DEPENDS ${GeantCudaLib})

if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${GeantCudaLib}.dylib DESTINATION lib)
else()
  install(TARGETS ${GeantCudaLib} DESTINATION lib)
endif()
