#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPPhysicsTableType.h"
#include "GPPhysics2DVector.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

//Material
#include "GPElement.h"
#include "GPMaterial.h"
#include "GPSandiaTable.h"

//Electron Processes/Models
#include "GPeBremsstrahlung.h"
#include "GPSeltzerBergerRelModel.h"
#include "GPeIonisation.h"
#include "GPMollerBhabhaModel.h"
#include "GPIonisParamMat.h"
#include "GPUniversalFluctuation.h"
#include "GPeMultipleScattering.h"
#include "GPUrbanMscModel95.h"

//Photon Processes/Models
#include "GPPhotonModel.h"
#include "GPPhotonProcess.h"
#include "GPEmProcessType.h"

//Manager
#include "GPProcessManager.h"
#include "GPTrackingManager.h"
#include "GPSteppingManager.h"

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
// transportation Kernel - no physics processes
//-----------------------------------------------------------------------------

GLOBALFUNC
void trans_kernel(curandState* devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPGeomManager *geomManager,
		  GXFieldMap *magMap,
		  GPPhysicsTable* physicsTable, 
		  GPPhysics2DVector* sbData,
		  GXTrack *secTracks, G4int *stackSize, G4int *offset,
		  G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(devStates,tid,sbData);
  GPeBremsstrahlung bremProcess(devStates,tid,
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  GPMollerBhabhaModel ioniModel(devStates,tid);
  GPeIonisation ioniProcess(devStates,tid, &physicsTable[kLambda_eIoni], 
                            &physicsTable[kRange_eIoni], 
                            &physicsTable[kDEDX_eIoni], 
                            &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  GPUrbanMscModel95 mscModel(devStates,tid);
  GPeMultipleScattering mscProcess(devStates,tid,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(0,
					     0,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for msc_kernel
//-----------------------------------------------------------------------------
void trans_gpu(curandState* devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable, 
	       GPPhysics2DVector* sbData,
	       GXTrack *secTracks, G4int *stackSize, G4int *offset,
	       G4int numStep, G4int runType, 
	       int blocksPerGrid, int threadsPerBlock,
	       cudaStream_t stream) 
{
  trans_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, offset, numStep,runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void trans_cpu(GXTrack *track, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable, 
	       GPPhysics2DVector* sbData,
	       GXTrack *secTracks, G4int *stackSize,
	       G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(0,-1,sbData);
  GPeBremsstrahlung bremProcess(0,-1, 
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  GPMollerBhabhaModel ioniModel(0,-1);
  GPeIonisation ioniProcess(0,-1, &physicsTable[kLambda_eIoni], 
			    &physicsTable[kRange_eIoni], 
			    &physicsTable[kDEDX_eIoni], 
			    &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  GPUrbanMscModel95 mscModel(0,-1);
  GPeMultipleScattering mscProcess(0,-1,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(0,
					     0,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// brem Kernel - GPeBremsstrahlung/GPSeltzerBergerRelModel
//-----------------------------------------------------------------------------

GLOBALFUNC
void brem_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPGeomManager *geomManager,
		 GXFieldMap *magMap,
		 GPPhysicsTable* physicsTable, 
		 GPPhysics2DVector* sbData,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models

  GPSeltzerBergerRelModel bremModel(devStates,tid,sbData);
  GPeBremsstrahlung bremProcess(devStates,tid,
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  /*
  GPMollerBhabhaModel ioniModel(devStates,tid);
  GPeIonisation ioniProcess(devStates,tid, &physicsTable[kLambda_eIoni], 
                            &physicsTable[kRange_eIoni], 
                            &physicsTable[kDEDX_eIoni], 
                            &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  GPUrbanMscModel95 mscModel(devStates,tid);
  GPeMultipleScattering mscProcess(devStates,tid,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     0,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for brem_kernel
//-----------------------------------------------------------------------------
void brem_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int numStep, G4int runType, 
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream) 
{
  brem_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, offset, numStep,runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void brem_cpu(GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  GPSeltzerBergerRelModel bremModel(0,-1,sbData);
  GPeBremsstrahlung bremProcess(0,-1, 
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  /*
  GPMollerBhabhaModel ioniModel(0,-1);
  GPeIonisation ioniProcess(0,-1, &physicsTable[kLambda_eIoni], 
			    &physicsTable[kRange_eIoni], 
			    &physicsTable[kDEDX_eIoni], 
			    &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  GPUrbanMscModel95 mscModel(0,-1);
  GPeMultipleScattering mscProcess(0,-1,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     0,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}
//-----------------------------------------------------------------------------
// ioni Kernel - GPeIonisation/GPMollerBhabhaModel
//-----------------------------------------------------------------------------

GLOBALFUNC
void ioni_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPGeomManager *geomManager,
		 GXFieldMap *magMap,
		 GPPhysicsTable* physicsTable, 
		 GPPhysics2DVector* sbData,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(devStates,tid,sbData);
  GPeBremsstrahlung bremProcess(devStates,tid,
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);
  */

  GPMollerBhabhaModel ioniModel(devStates,tid);
  GPeIonisation ioniProcess(devStates,tid, &physicsTable[kLambda_eIoni], 
                            &physicsTable[kRange_eIoni], 
                            &physicsTable[kDEDX_eIoni], 
                            &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  /*
  GPUrbanMscModel95 mscModel(devStates,tid);
  GPeMultipleScattering mscProcess(devStates,tid,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(0,
					     &ioniProcess,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for ioni_kernel
//-----------------------------------------------------------------------------
void ioni_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int numStep, G4int runType, 
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream) 
{
  ioni_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, offset, numStep,runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void ioni_cpu(GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(0,-1,sbData);
  GPeBremsstrahlung bremProcess(0,-1, 
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);
  */

  GPMollerBhabhaModel ioniModel(0,-1);
  GPeIonisation ioniProcess(0,-1, &physicsTable[kLambda_eIoni], 
			    &physicsTable[kRange_eIoni], 
			    &physicsTable[kDEDX_eIoni], 
			    &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);

  /*
  GPUrbanMscModel95 mscModel(0,-1);
  GPeMultipleScattering mscProcess(0,-1,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);
  */

  aElectronProcessManager.AddElectronProcess(0,
					     &ioniProcess,
					     0);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// msc Kernel - GPeMultipleScattering/GPUrbanMscModel95
//-----------------------------------------------------------------------------

GLOBALFUNC
void msc_kernel(curandState* devStates,
		GXTrack *track, size_t nTrackSize,
		GPGeomManager *geomManager,
		GXFieldMap *magMap,
		GPPhysicsTable* physicsTable, 
		GPPhysics2DVector* sbData,
		GXTrack *secTracks, G4int *stackSize, G4int *offset,
		G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(devStates,tid,sbData);
  GPeBremsstrahlung bremProcess(devStates,tid,
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  GPMollerBhabhaModel ioniModel(devStates,tid);
  GPeIonisation ioniProcess(devStates,tid, &physicsTable[kLambda_eIoni], 
                            &physicsTable[kRange_eIoni], 
                            &physicsTable[kDEDX_eIoni], 
                            &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);
  */

  GPUrbanMscModel95 mscModel(devStates,tid);
  GPeMultipleScattering mscProcess(devStates,tid,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);

  aElectronProcessManager.AddElectronProcess(0,
					     0,
					     &mscProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for msc_kernel
//-----------------------------------------------------------------------------
void msc_gpu(curandState* devStates,
	     GXTrack *track, size_t nTrackSize,
	     GPGeomManager *geomManager,
	     GXFieldMap *magMap,
	     GPPhysicsTable* physicsTable, 
	     GPPhysics2DVector* sbData,
	     GXTrack *secTracks, G4int *stackSize, G4int *offset,
	     G4int numStep, G4int runType, 
	     int blocksPerGrid, int threadsPerBlock,
	     cudaStream_t stream) 
{
  msc_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, offset, numStep,runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void msc_cpu(GXTrack *track, size_t nTrackSize,
	     GPGeomManager *geomManager,
	     GXFieldMap *magMap,
	     GPPhysicsTable* physicsTable, 
	     GPPhysics2DVector* sbData,
	     GXTrack *secTracks, G4int *stackSize,
	     G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aElectronProcessManager;

  //EM processes and models
  /*
  GPSeltzerBergerRelModel bremModel(0,-1,sbData);
  GPeBremsstrahlung bremProcess(0,-1, 
				&physicsTable[kLambda_eBrem]);
  bremProcess.InitialiseProcess(&bremModel);

  GPMollerBhabhaModel ioniModel(0,-1);
  GPeIonisation ioniProcess(0,-1, &physicsTable[kLambda_eIoni], 
			    &physicsTable[kRange_eIoni], 
			    &physicsTable[kDEDX_eIoni], 
			    &physicsTable[kInverseRange_eIoni]);
  ioniProcess.InitialiseProcess(&ioniModel);
  */

  GPUrbanMscModel95 mscModel(0,-1);
  GPeMultipleScattering mscProcess(0,-1,&physicsTable[kLambda_msc]);
  mscProcess.InitialiseProcess(&mscModel);

  aElectronProcessManager.AddElectronProcess(0,
					     0,
					     &mscProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// compt Kernel 
//-----------------------------------------------------------------------------

GLOBALFUNC
void compt_kernel(curandState* devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPGeomManager *geomManager,
		  GXFieldMap *magMap,
		  GPPhysicsTable* physicsTable, 
		  GXTrack *secTracks, G4int *stackSize, G4int *offset,
		  G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aPhotonProcessManager;

  //Construct photon processes/models
  GPPhotonModel comptModel(devStates,tid);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(devStates,tid,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);

  /*
  GPPhotonModel convModel(devStates,tid);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(devStates,tid,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);

  GPPhotonModel photModel(devStates,tid);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(devStates,tid,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);
  */

  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void compt_gpu(curandState* devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable, 
	       GXTrack *secTracks, G4int *stackSize, G4int *offset,
	       G4int numStep, G4int runType,
	       int blocksPerGrid, int threadsPerBlock,
	       cudaStream_t stream)
{
  compt_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize, geomManager, magMap, physicsTable,
     secTracks, stackSize, offset, numStep, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void compt_cpu(GXTrack *track, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable, 
	       GXTrack *secTracks, G4int *stackSize,
	       G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aPhotonProcessManager;

  GPPhotonModel comptModel(0,-1);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(0,-1,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);

  /*
  GPPhotonModel convModel(0,-1);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(0,-1,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);

  //PhotoElectricEffect processes
  GPPhotonModel photModel(0,-1);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(0,-1,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);
  */

  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
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

  for (int tid = 0; tid < nTrackSize; ++tid) {
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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// conv Kernel 
//-----------------------------------------------------------------------------

GLOBALFUNC
void conv_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPGeomManager *geomManager,
		 GXFieldMap *magMap,
		 GPPhysicsTable* physicsTable, 
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aPhotonProcessManager;

  //Construct photon processes/models
  /*
  GPPhotonModel comptModel(devStates,tid);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(devStates,tid,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);
  */

  GPPhotonModel convModel(devStates,tid);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(devStates,tid,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);

  /*
  GPPhotonModel photModel(devStates,tid);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(devStates,tid,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);
  */

  //  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void conv_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int numStep, G4int runType,
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream)
{
  conv_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize, geomManager, magMap, physicsTable,
     secTracks, stackSize, offset, numStep, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void conv_cpu(GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aPhotonProcessManager;

  /*
  GPPhotonModel comptModel(0,-1);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(0,-1,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);
  */

  GPPhotonModel convModel(0,-1);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(0,-1,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);

  //PhotoElectricEffect processes
  /*
  GPPhotonModel photModel(0,-1);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(0,-1,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);
  */

  //  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
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

  for (int tid = 0; tid < nTrackSize; ++tid) {
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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// compt Kernel 
//-----------------------------------------------------------------------------

GLOBALFUNC
void phot_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPGeomManager *geomManager,
		 GXFieldMap *magMap,
		 GPPhysicsTable* physicsTable, 
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GPProcessManager aPhotonProcessManager;

  //Construct photon processes/models
  /*
  GPPhotonModel comptModel(devStates,tid);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(devStates,tid,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);

  GPPhotonModel convModel(devStates,tid);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(devStates,tid,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);
  */

  GPPhotonModel photModel(devStates,tid);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(devStates,tid,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);

  //  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void phot_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int numStep, G4int runType,
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream)
{
  phot_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, nTrackSize, geomManager, magMap, physicsTable,
     secTracks, stackSize, offset, numStep, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void phot_cpu(GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable, 
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType)
{
  //process manager
  GPProcessManager aPhotonProcessManager;

  /*
  GPPhotonModel comptModel(0,-1);
  comptModel.SetProcessType(kCompton);

  GPPhotonProcess comptProcess(0,-1,&physicsTable[kLambda_compt]);
  comptProcess.InitialiseProcess(kCompton,&comptModel);

  GPPhotonModel convModel(0,-1);
  convModel.SetProcessType(kConversion);

  GPPhotonProcess convProcess(0,-1,&physicsTable[kLambda_conv]);
  convProcess.InitialiseProcess(kConversion,&convModel);
  */

  //PhotoElectricEffect processes
  GPPhotonModel photModel(0,-1);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  GPPhotonProcess photProcess(0,-1,&physicsTable[kLambdaPrim_phot]);
  photProcess.InitialiseProcess(kPhotoElectric,&photModel);

  //  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  //  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
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

  for (int tid = 0; tid < nTrackSize; ++tid) {
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

    //    aSteppingManager.SetTransportationProcess(&transp);
    aSteppingManager.SetTransportation(&transp,&aNavigator);

    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}


//-----------------------------------------------------------------------------
// Other kernels
//-----------------------------------------------------------------------------
#include "random_kernel.cu"

//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPThreeVector.cu"
#include "GPThreeVectorList.cu"
#include "GPRotationMatrix.cu"
#include "GPUtils.cu"
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"
#include "GPVParticleChange.cu"

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
#include "GPeBremsstrahlung.cu"
#include "GPSeltzerBergerRelModel.cu"
#include "GPeIonisation.cu"
#include "GPMollerBhabhaModel.cu"
#include "GPUniversalFluctuation.cu"
#include "GPeMultipleScattering.cu"
#include "GPUrbanMscModel95.cu"

//-----------------------------------------------------------------------------
// Photon processes/models
//-----------------------------------------------------------------------------
#include "GPPhotonModel.cu"
#include "GPPhotonProcess.cu"

//-----------------------------------------------------------------------------
// Manager
//-----------------------------------------------------------------------------
#include "GPProcessManager.cu"
#include "GPTrackingManager.cu"
#include "GPSteppingManager.cu"

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
