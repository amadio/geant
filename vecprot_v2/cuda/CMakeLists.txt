#------------------------------------------------------------------------------
# - CMakeLists.txt for CUDA coprocessor broker
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.8)

# Steer the name of the library (main library replacement or
# plug-in).

set(GeantCudaLib Geant_v) # or GeantCuda
set(GeantCudaReplace ON)

set(CudaCode ${PROJECT_SOURCE_DIR})
set(CUDA_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CUDA_INCDIR ${CMAKE_CURRENT_SOURCE_DIR}/inc)

#NVCC_FLAGS=--device-debug -m64 -Xcompiler ,\"-g\" -Xcompiler ,\"-fPIC\"  -arch=sm_20 --use_fast_math --ptxas-options=-v -DNVCC -DGX_MULTI_OBJ_FILES

if (CUDA_SEPARABLE_COMPILATION)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA_SEP_COMP -DGX_MULTI_OBJ_FILES)
endif()


CUDA_INCLUDE_DIRECTORIES( ${CudaCode}/GXTracking/SimCore/Common/include 
       ${CudaCode}/GXTracking/SimCore/EMPhysics/include
       ${CudaCode}/GXTracking/SimCore/Geometry/include
       ${CudaCode}/GXTracking/SimCore/Transportation/include
       ${CudaCode}/GXTracking/SimCore/Manager/include
       ${CudaCode}/GXTracking/FWCore/FrameWork/include
       ${CudaCode}/GXTracking/FWCore/Application/include
       ${CudaCode}/GXTracking/GPTracking/EMPhysics/include
       ${CudaCode}/GXTracking/gpu/cuda/kernel/include
       ${CudaCode}/GXTracking/gpu/cuda/kernel/src
       ${CudaCode}/GXTracking/GPTracking/Common/include
       ${CudaCode}/GXTracking/GPTracking/Geometry/include
       ${CudaCode}/GXTracking/GPTracking/Material/include
       ${CudaCode}/GXTracking/GPTracking/Transportation/include
       ${CudaCode}/GXTracking/GPTracking/Common/src
       ${CudaCode}/GXTracking/GPTracking/EMPhysics/src
       ${CudaCode}/GXTracking/GPTracking/Geometry/src
       ${CudaCode}/GXTracking/GPTracking/Material/src
       ${CudaCode}/GXTracking/GPTracking/Transportation/src
       ${CudaCode}/GXTracking/SimCore/Common/src
       ${CudaCode}/GXTracking/SimCore/EMPhysics/src
       ${CudaCode}/GXTracking/SimCore/Geometry/src
       ${CudaCode}/GXTracking/SimCore/Transportation/src
       ${CudaCode}/GXTracking/SimCore/Manager/src
       ${CudaCode}/GXTracking/FWCore/FrameWork/src
       ${CudaCode}/GXTracking/FWCore/Application/src
)

#CXXFLAGS += -I$(CUDA_INC) -I$(CUDA_INCDIR)

INCLUDE_DIRECTORIES( 
       ${CudaCode}/GXTracking/SimCore/Common/include
       ${CudaCode}/GXTracking/gpu/cuda/kernel/include
       ${CudaCode}/GXTracking/SimCore/EMPhysics/include
       ${CudaCode}/GXTracking/SimCore/Transportation/include
       ${CudaCode}/GXTracking/SimCore/Manager/include
       ${CudaCode}/GXTracking/GPTracking/Common/include
       ${CudaCode}/GXTracking/GPTracking/Transportation/include
       ${CudaCode}/GXTracking/GPTracking/EMPhysics/include
       ${CudaCode}/GXTracking/GPTracking/Geometry/include
       ${CudaCode}/GXTracking/GPTracking/Material/include
       ${CUDA_INCDIR}
)


set(CUDA_DATA_TARFILE ${PROJECT_SOURCE_DIR}/data/vpdata_v1.tar.gz)
set(CUDA_DATA_FILE ${PROJECT_SOURCE_DIR}/data/Lambda.eBrem.e-.asc)

set(dict_headers ${CUDA_INCDIR}/CoprocessorBroker.h)

set(cuda_cxx_sources
    ${CUDA_SRCDIR}/CoprocessorBroker.cxx
    G__GeantCuda.cxx
)


if (CUDA_SEPARABLE_COMPILATION)
set(cuda_sources
    ${CudaCode}/GXTracking/gpu/cuda/kernel/src/trackingTest2_kernel.cu
    ${CudaCode}/GXTracking/gpu/cuda/kernel/src/random_kernel.cu
    ${CudaCode}/GXTracking/gpu/cuda/kernel/src/sort_kernel.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXBrem.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXIoni.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXMsc.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXeBremsstrahlung.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXeMultipleScattering.cu
    ${CudaCode}/GXTracking/SimCore/EMPhysics/src/GXeIonisation.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXMagneticField.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXEquationOfMotion.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXClassicalRK4.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXPropagatorInField.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXTransportation.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXFieldTrack.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXFieldManager.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXChordFinder.cu
    ${CudaCode}/GXTracking/SimCore/Transportation/src/GXMagIntegratorDriver.cu
    ${CudaCode}/GXTracking/SimCore/Geometry/src/GXMultiLevelLocator.cu
    ${CudaCode}/GXTracking/SimCore//Common/src/GXPhysics2DVector.cu
    ${CudaCode}/GXTracking/SimCore//Common/src/GXPhysicsTable.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXeProcessManager.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXeSteppingManager.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXeTrackingManager.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXProcessManager.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXSteppingManager.cu
    ${CudaCode}/GXTracking/SimCore//Manager/src/GXTrackingManager.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPThreeVector.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPUtils.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPThreeVectorList.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPRotationMatrix.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPVParticleChange.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPStep.cu
    ${CudaCode}/GXTracking/GPTracking/Common/src/GPStepPoint.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPFieldTrack.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPChordFinder.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPChargeState.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPMagIntegratorDriver.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPEquationOfMotion.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPClassicalRK4.cu
    ${CudaCode}/GXTracking/GPTracking/Transportation/src/GPMagneticField.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPPhysicsTable.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPPhotonProcess.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPeMultipleScattering.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPeIonisation.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPeBremsstrahlung.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPSeltzerBergerRelModel.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPMollerBhabhaModel.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPPhotonModel.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPUniversalFluctuation.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPUrbanMscModel95.cu
    ${CudaCode}/GXTracking/GPTracking/EMPhysics/src/GPPhysics2DVector.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPAffineTransform.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPAuxiliaryNavServices.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPBox.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPCombinedNavigation.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPCons.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPLineSection.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPLogicalVolume.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPNavigationHistory.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPNavigationLevel.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPNavigator.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPNormalNavigation.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPOrb.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPSmartVoxelHeader.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPSmartVoxelNode.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPSmartVoxelProxy.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPTouchableHistory.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPTrd.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPTubs.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPUserGeometry.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPVoxelHeader.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPVoxelLimits.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPVoxelNavigation.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPVPhysicalVolume.cu
    ${CudaCode}/GXTracking/GPTracking/Geometry/src/GPVSolid.cu
    ${CudaCode}/GXTracking/GPTracking/Material/src/GPAtomicShells.cu
    ${CudaCode}/GXTracking/GPTracking/Material/src/GPElement.cu
    ${CudaCode}/GXTracking/GPTracking/Material/src/GPIonisParamMat.cu
    ${CudaCode}/GXTracking/GPTracking/Material/src/GPMaterial.cu
    ${CudaCode}/GXTracking/GPTracking/Material/src/GPSandiaTable.cu

    ${CUDA_SRCDIR}/tracking_kernel.cu
    ${CUDA_SRCDIR}/util_kernel.cu
    ${CUDA_SRCDIR}/electron_kernel.cu
)
else()
set(cuda_sources
#    ${CudaCode}/GXTracking/gpu/cuda/kernel/src/trackingTest2_kernel.cu
    ${CUDA_SRCDIR}/electron_kernel.cu
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
   list(REMOVE_ITEM  cuda_cxx_sources ${CMAKE_HOME_DIRECTORY}/vecprot_v2/src/GeantTrack.cxx)
   list(REMOVE_ITEM  cuda_cxx_sources ${CMAKE_HOME_DIRECTORY}/vecprot_v2/src/GeantTaskData.cxx)

   set (extra_dict_options -DGEANTCUDA_REPLACE)
endif()

ROOT_GENERATE_DICTIONARY(G__GeantCuda ${dict_headers}
                MODULE lib${GeantCudaLib}
                LINKDEF ${CUDA_INCDIR}/CudaLinkDef.h
                OPTIONS -I${CUDA_INCLUDE_DIRS} ${extra_dict_options} )

# add the command to generate the source code
add_custom_command (
  OUTPUT ${CUDA_DATA_FILE}
  COMMAND (cd data; tar xfz $(notdir $(CUDA_DATA_TARFILE)) ) && touch $@
  DEPENDS ${CUDA_DATA_TARFILE}
)

add_custom_command (
  OUTPUT ${CUDA_DATA_TARFILE}
  COMMAND mkdir -p data && wet https://oink.fnal.gov/perfanalysis/vp${OUTPUT} -O ${OUTPUT}
#  DEPENDS 
)


CUDA_ADD_LIBRARY( ${GeantCudaLib}
                  ${cuda_sources} ${cuda_cxx_sources}
                 SHARED
                 OPTIONS ${CUDA_ARCH}
)

if( VECGEOM_FOUND )
  target_link_libraries(${GeantCudaLib} ${ROOT_LIBRARIES} -lEG -lGeom -lGui  ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
else()
  target_link_libraries(${GeantCudaLib} ${ROOT_LIBRARIES} -lEG -lGeom -lGui )
endif()

add_custom_target(Geantlib_v DEPENDS ${GeantCudaLib})

if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${GeantCudaLib}.dylib DESTINATION lib)
else()
  install(TARGETS ${GeantCudaLib} DESTINATION lib)
endif()
