#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GXTrackLiason.h"
#include "GXPhysicsTable.h"
#include "GXPhysicsTableType.h"
#include "GXPhysics2DVector.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

//Material
#include "GPElement.h"
#include "GPMaterial.h"
#include "GPSandiaTable.h"

//Electron Processes/Models
#include "GXeBremsstrahlung.h"
#include "GXeIonisation.h"
#include "GXeMultipleScattering.h"

#include "GPSeltzerBergerRelModel.h"
#include "GPMollerBhabhaModel.h"
#include "GPIonisParamMat.h"
#include "GPUniversalFluctuation.h"
#include "GPUrbanMscModel95.h"

//Manager
#include "GXeProcessManager.h"
#include "GXeTrackingManager.h"
#include "GXeSteppingManager.h"

//Geometry
#include "GPGeomManager.h"
#include "GPVPhysicalVolume.h"
#include "GPTouchableHistory.h"
#include "GPNavigator.h"

//Transportation
#include "GXFieldMap.h"
#include "GXFieldTrack.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"
#include "GXMagIntegratorDriver.h"
#include "GXChordFinder.h"
#include "GXFieldManager.h"
#include "GXPropagatorInField.h"
#include "GXMultiLevelLocator.h"
#include "GXTransportation.h"

#include "stdio.h"

//-----------------------------------------------------------------------------
// Tracking DoIt Kernel
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Tracking GPIL Kernel
//-----------------------------------------------------------------------------

GLOBALFUNC
void elec_GPIL_kernel(curandState* devStates,
		      GXTrack *track, 
		      GXTrackLiason *liason, 
		      size_t nTrackSize,
		      GPGeomManager *geomManager,
		      GXFieldMap *magMap,
		      GXPhysicsTable* physicsTable, 
		      GXPhysics2DVector* sbData,
		      GXTrack *secTracks, G4int *stackSize,
		      G4int numStep) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GXeProcessManager aElectronProcessManager;

  //EM processes and models

  GXeBremsstrahlung bremProcess(devStates,tid,
				&physicsTable[kLambda_eBrem]);

  GXeIonisation ioniProcess(devStates,tid, &physicsTable[kLambda_eIoni], 
                            &physicsTable[kRange_eIoni], 
                            &physicsTable[kDEDX_eIoni], 
                            &physicsTable[kInverseRange_eIoni]);

  GXeMultipleScattering mscProcess(devStates,tid,&physicsTable[kLambda_msc]);

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXeSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXeTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  //Transportation
  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  GXMagneticField magField;
  GXMagneticField_Constructor(&magField,magMap);
  
  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
  
  GXClassicalRK4 rk4;
  GXClassicalRK4_Constructor(&rk4,&equaOfMotion);
  
  GXMagInt_Driver magDriver;
  GXMagInt_Driver_Constructor(&magDriver,1.0,&rk4);

  GXChordFinder chordFinder;
  GXChordFinder_Constructor(&chordFinder,&magDriver);

  GXFieldManager aFieldManager;
  GXFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GXMultiLevelLocator mLocator;

  //Propagator
  GXPropagatorInField propagatorInField;

  //Transporation
  GXTransportation transp;

  while (tid < nTrackSize) {

    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    aSteppingManager.SetSecondaryStack(secTracks,stackSize);

    aTrackingManager.ProcessGPIL(&track[tid]);

    //save material information for the DoIt step
    track[tid].id = tid;
    //    liason[tid].material = aSteppingManager.GetMaterial();

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for msc_kernel
//-----------------------------------------------------------------------------
void elec_GPIL_gpu(curandState* devStates,
		   GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GXPhysicsTable* physicsTable, 
		   GXPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep,
		   int blocksPerGrid, int threadsPerBlock,
		   cudaStream_t stream) 
{
  elec_GPIL_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, liason, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, numStep);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void elec_GPIL_cpu(GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GXPhysicsTable* physicsTable, 
		   GXPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep)
{
  //process manager
  GXeProcessManager aElectronProcessManager;

  //EM processes and models

  GXeBremsstrahlung bremProcess(0,-1, 
				&physicsTable[kLambda_eBrem]);

  GXeIonisation ioniProcess(0,-1, &physicsTable[kLambda_eIoni], 
			    &physicsTable[kRange_eIoni], 
			    &physicsTable[kDEDX_eIoni], 
			    &physicsTable[kInverseRange_eIoni]);

  GXeMultipleScattering mscProcess(0,-1,&physicsTable[kLambda_msc]);

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXeSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXeTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  //Transportation
  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  GXMagneticField magField;
  GXMagneticField_Constructor(&magField,magMap);
  
  GXEquationOfMotion equaOfMotion;
  GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
  
  GXClassicalRK4 rk4;
  GXClassicalRK4_Constructor(&rk4,&equaOfMotion);
  
  GXMagInt_Driver magDriver;
  GXMagInt_Driver_Constructor(&magDriver,1.0,&rk4);

  GXChordFinder chordFinder;
  GXChordFinder_Constructor(&chordFinder,&magDriver);

  GXFieldManager aFieldManager;
  GXFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GXMultiLevelLocator mLocator;

  //Propagator
  GXPropagatorInField propagatorInField;

  //Transporation
  GXTransportation transp;

  for (int tid = 0; tid < nTrackSize; tid++) {
    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    aSteppingManager.SetSecondaryStack(secTracks,stackSize);

    aTrackingManager.ProcessGPIL(&track[tid]);
    //save transite information for the DoIt step
    track[tid].id = tid;
//    liason[tid].material = aSteppingManager.GetMaterial();
  }
}

#ifndef GX_MULTI_OBJ_FILES

//-----------------------------------------------------------------------------
// Other kernels
//-----------------------------------------------------------------------------
#include "random_kernel.cu"
#include "sort_kernel.cu"

//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPThreeVector.cu"
#include "GPThreeVectorList.cu"
#include "GPRotationMatrix.cu"
#include "GPUtils.cu"
#include "GPVParticleChange.cu"
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"

#include "GXPhysicsTable.cu"
#include "GXPhysics2DVector.cu"

//-----------------------------------------------------------------------------
// Material
//-----------------------------------------------------------------------------
#include "GPElement.cu"
#include "GPMaterial.cu"
#include "GPStep.cu"
#include "GPStepPoint.cu"
#include "GPAtomicShells.cu"
#include "GPSandiaTable.cu"
#include "GPIonisParamMat.cu"

//-----------------------------------------------------------------------------
// Electron processes/models
//-----------------------------------------------------------------------------
#include "GXeBremsstrahlung.cu"
#include "GXeIonisation.cu"
#include "GXeMultipleScattering.cu"

#include "GPeBremsstrahlung.cu"
#include "GPeIonisation.cu"
#include "GPeMultipleScattering.cu"
#include "GPSeltzerBergerRelModel.cu"
#include "GPMollerBhabhaModel.cu"
#include "GPUniversalFluctuation.cu"
#include "GPUrbanMscModel95.cu"
//-----------------------------------------------------------------------------
// Manager
//-----------------------------------------------------------------------------
#include "GXeProcessManager.cu"
#include "GXeTrackingManager.cu"
#include "GXeSteppingManager.cu"

//-----------------------------------------------------------------------------
// GPGeometry
//-----------------------------------------------------------------------------
#include "GPAffineTransform.cu"
#include "GPAuxiliaryNavServices.cu"
#include "GPCombinedNavigation.cu"
#include "GPLineSection.cu"
#include "GPLogicalVolume.cu"
#include "GPNavigationLevel.cu"
#include "GPNavigator.cu"
#include "GPNavigationHistory.cu"
#include "GPNormalNavigation.cu"
#include "GPSmartVoxelHeader.cu"
#include "GPSmartVoxelNode.cu"
#include "GPSmartVoxelProxy.cu"
#include "GPTouchableHistory.cu"
#include "GPVoxelLimits.cu"
#include "GPVoxelNavigation.cu"
#include "GPVPhysicalVolume.cu"
#include "GPBox.cu"
#include "GPCons.cu"
#include "GPOrb.cu"
#include "GPTrd.cu"
#include "GPTubs.cu"
#include "GPVSolid.cu"
#include "GPChargeState.cu"
#include "GPUserGeometry.cu"
#include "GPVoxelHeader.cu"

//-----------------------------------------------------------------------------
// Magnetic Field and Transporation
//-----------------------------------------------------------------------------
#include "GXFieldTrack.cu"
#include "GXFieldManager.cu"
#include "GXMagneticField.cu"
#include "GXClassicalRK4.cu"
#include "GXEquationOfMotion.cu"
#include "GXMagIntegratorDriver.cu"
#include "GXChordFinder.cu"
#include "GXPropagatorInField.cu"
#include "GXMultiLevelLocator.cu"
#include "GXTransportation.cu"

//-----------------------------------------------------------------------------
// FrameWork and Application
//-----------------------------------------------------------------------------
#include "GXRunManager.cc"
#include "GXTrackHandler.cc"
#include "GXUserGeometry.cc"
#include "GXSimpleEcal.cc"
#include "GXGPUManager.cc"

#endif