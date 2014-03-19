#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GXTrackLiason.h"
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
#include "GXProcessManager.h"
#include "GXTrackingManager.h"
#include "GXSteppingManager.h"

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

GLOBALFUNC
void elec_doit_kernel(curandState* devStates,
		      GXTrack *track, 
		      GXTrackLiason *liason, 
		      size_t nTrackSize,
		      GPGeomManager *geomManager,
		      GXFieldMap *magMap,
		      GPPhysicsTable* physicsTable, 
		      GPPhysics2DVector* sbData,
		      GXTrack *secTracks, G4int *stackSize,
		      G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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
  
  //  GPeMultipleScattering *mscProcess = liason[track[tid].id].msc;

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);


  //  GXTrackingManager *aTrackingManager = liason[track[tid].id]->trackingManager;

  //Transportation
  /*
  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //Transporation
  GXTransportation transp;
  */

  while (tid < nTrackSize) {

    //    if(track[tid].proc != 0 && track[tid].proc != 1) {
    //      printf("!!!!Check This Event proc = %d\n",track[tid].proc);
    //    }

    //Geometry - Initialize Navigator and construct related structures
    /*
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);

    GXTransportation_Constructor2(&transp,&aNavigator,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    */
    aSteppingManager.SetSecondaryStack(secTracks,stackSize);

    //retrive information from GPIL
    aSteppingManager.SetMaterial((liason[track[tid].id]).material);

    aTrackingManager.ProcessDoIt(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for GPU DoIt
//-----------------------------------------------------------------------------
void elec_doit_gpu(curandState* devStates,
		   GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GPPhysicsTable* physicsTable, 
		   GPPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep, G4int runType, 
		   int blocksPerGrid, int threadsPerBlock,
		   cudaStream_t stream) 
{
  elec_doit_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, liason, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, numStep,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU DoIt
//-----------------------------------------------------------------------------

void elec_doit_cpu(GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GPPhysicsTable* physicsTable, 
		   GPPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep, G4int runType)
{
  int tid = -1;
  curandState* devStates = 0;

  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  /*
  //Transportation
  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //Transporation
  GXTransportation transp;
  */

  for (int tid = 0; tid < nTrackSize; tid++) {

    //    if(track[tid].proc != 0 && track[tid].proc != 1) {
    //      printf("!!!!Check This Event proc = %d\n",track[tid].proc);
    //    }

    //Geometry - Initialize Navigator and construct related structures
    /*
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);

    GXTransportation_Constructor2(&transp,&aNavigator,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    */

    aSteppingManager.SetSecondaryStack(secTracks,stackSize);

    //retrive information from GPIL
    aSteppingManager.SetMaterial((liason[track[tid].id]).material);
    //    aElectronProcessManager.ResetMultipleScattering((liason[track[tid].id]).msc);  

    aTrackingManager.ProcessDoIt(&track[tid]);

  }
}

//-----------------------------------------------------------------------------
// Tracking GPIL Kernel
//-----------------------------------------------------------------------------

GLOBALFUNC
void elec_GPIL_kernel(curandState* devStates,
		      GXTrack *track, 
		      GXTrackLiason *liason, 
                      int *logVolumeIndices,
                      int *physVolumeIndices,
		      size_t nTrackSize,
		      GPGeomManager *geomManager,
		      GXFieldMap *magMap,
		      GPPhysicsTable* physicsTable, 
		      GPPhysics2DVector* sbData,
		      GXTrack *secTracks, G4int *stackSize,
		      G4int numStep, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
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
    
    // We know where we are coming from, so we need to upgrade the
    // following code to use the location information we already have.
    GPVPhysicalVolume *where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);
    if (physVolumeIndices && logVolumeIndices) {
       if (where) {
          double x    = track[tid].x  ;
          double y    = track[tid].y  ;
          double z    = track[tid].z  ;
          if (where->fIndex != physVolumeIndices[tid])
             printf("ERROR: %d Found the track to be in phys volume node #%d expected %d at %lf %lf %lf\n",
                    tid,where->fIndex,physVolumeIndices[tid],
                    x,y,z);
          if (where->flogical->fIndex != logVolumeIndices[tid])
             printf("ERROR: %d Found the track to be in log volume node #%d expected %d at %lf %lf %lf\n",
                    tid,where->flogical->fIndex,logVolumeIndices[tid],
                    x,y,z);
       } else printf("ERROR: Could not find where the track is ...\n");
    }
    
    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    aSteppingManager.SetSecondaryStack(secTracks,stackSize);

    aTrackingManager.ProcessGPIL(&track[tid]);

    //save material information for the DoIt step
    track[tid].id = tid;
    liason[tid].material = aSteppingManager.GetMaterial();
    //aElectronProcessManager.GetMultipleScatteringProcess();
    
    //    printf("Got %d: %f %f %f exit: %d vol: %d %ld\n",tid,
    //           atrack.x,atrack.y,atrack.z,
    //           aNavigator.fExiting,
    //           aNavigator.fBlockedPhysicalVolume ? aNavigator.fBlockedPhysicalVolume->fIndex : -1,
    //           aNavigator.fBlockedPhysicalVolume);
    // if (aNavigator.fBlockedPhysicalVolume) {
    //    physVolumeIndices[tid] = aNavigator.fBlockedPhysicalVolume->fIndex;
    //    logVolumeIndices[tid] = aNavigator.fBlockedPhysicalVolume->flogical->fIndex;
    //    // printf("Moved from %d to in %d\n",logVolumeIndices[tid],aNavigator.fBlockedPhysicalVolume ? aNavigator.fBlockedPhysicalVolume->flogical->fIndex : -1);
    // } else {
    // humm so we have not moved ... let's check..
    if (physVolumeIndices && logVolumeIndices) {
       double x    = track[tid].x  ;
       double y    = track[tid].y  ;
       double z    = track[tid].z  ;
       where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),NULL,false,true);
       if (!where) {
          // printf("Moved out of the geometry\n");
          physVolumeIndices[tid] = -1;
          logVolumeIndices[tid] = -1;
       } else {
          //            if (where->flogical->fIndex != logVolumeIndices[tid] || where->fIndex != physVolumeIndices[tid]) {
          //                printf("Moved from %d:%d to %d:%d\n", logVolumeIndices[tid], physVolumeIndices[tid],
          //                      where->flogical->fIndex, where->fIndex);
          //            }
          physVolumeIndices[tid] = where->fIndex;
          logVolumeIndices[tid] = where->flogical->fIndex;
       }            
    }
    // }

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
		   GPPhysicsTable* physicsTable, 
		   GPPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep, G4int runType, 
		   int blocksPerGrid, int threadsPerBlock,
		   cudaStream_t stream) 
{
  elec_GPIL_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, liason, 0, 0, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, numStep,runType);
}

void elec_GPIL_gpu(curandState* devStates,
		   GXTrack *track, 
		   GXTrackLiason *liason, 
                   int *logVolumeIndices,
                   int *physVolumeIndices,
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GPPhysicsTable* physicsTable, 
		   GPPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep, G4int runType, 
		   int blocksPerGrid, int threadsPerBlock,
		   cudaStream_t stream) 
{
  elec_GPIL_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, liason, logVolumeIndices, physVolumeIndices, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, numStep, runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void elec_GPIL_cpu(GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GPPhysicsTable* physicsTable, 
		   GPPhysics2DVector* sbData,
		   GXTrack *secTracks, G4int *stackSize,
		   G4int numStep, G4int runType)
{
  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);  

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
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
    liason[tid].material = aSteppingManager.GetMaterial();
    //    liason[tid].msc = aElectronProcessManager.GetMultipleScatteringProcess();
  }
}

//-----------------------------------------------------------------------------
// Tracking Kernel (Original)
//-----------------------------------------------------------------------------

GLOBALFUNC
void elec_kernel(curandState* devStates,
		 GXTrack *track,
                 int *logVolumeIndices,
                 int *physVolumeIndices,
		 size_t nTrackSize,
		 GPGeomManager *geomManager,
		 GXFieldMap *magMap,
		 GPPhysicsTable* physicsTable,
		 GPPhysics2DVector* sbData,
		 GXTrack *secTracks, G4int *stackSize,
		 G4int numStep, G4int runType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
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

    // We know where we are coming from, so we need to upgrade the
    // following code to use the location information we already have.
    GPVPhysicalVolume *where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,
					  GPThreeVector_create(track[tid].x,
							       track[tid].y,
							       track[tid].z),
                                          NULL,false,true);
    if (physVolumeIndices && logVolumeIndices) {
      if (where) {
        double x    = track[tid].x  ;
        double y    = track[tid].y  ;
        double z    = track[tid].z  ;
        if (where->fIndex != physVolumeIndices[tid])
          printf("ERROR: %d Found the track to be in phys volume node #%d expected %d at %lf %lf %lf\n",
                 tid,where->fIndex,physVolumeIndices[tid],
                 x,y,z);
        if (where->flogical->fIndex != logVolumeIndices[tid])
          printf("ERROR: %d Found the track to be in log volume node #%d expected %d at %lf %lf %lf\n",
                 tid,where->flogical->fIndex,logVolumeIndices[tid],
                 x,y,z);
      } else printf("ERROR: Could not find where the track is ...\n");
    }

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);

    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                    &aFieldManager,&mLocator);

    GXTransportation_Constructor(&transp,&propagatorInField,0);

    aSteppingManager.SetTransportation(&transp,&aNavigator);
    aSteppingManager.SetSecondaryStack(secTracks,stackSize);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    if (physVolumeIndices && logVolumeIndices) {
       double x    = track[tid].x  ;
       double y    = track[tid].y  ;
       double z    = track[tid].z  ;
       where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),NULL,false,true);
       if (!where) {
          // printf("Moved out of the geometry\n");
          physVolumeIndices[tid] = -1;
          logVolumeIndices[tid] = -1;
       } else {
          physVolumeIndices[tid] = where->fIndex;
          logVolumeIndices[tid] = where->flogical->fIndex;
       }
    }

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for msc_kernel
//-----------------------------------------------------------------------------
void elec_gpu(curandState* devStates,
	      GXTrack *track,
	      size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType,
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream)
{
  elec_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, 0,0, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize,numStep,runType);
}
void elec_gpu(curandState* devStates,
	      GXTrack *track,
              int *logVolumeIndices,
              int *physVolumeIndices,
	      size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType,
	      int blocksPerGrid, int threadsPerBlock,
	      cudaStream_t stream)
{
  elec_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,track, logVolumeIndices, physVolumeIndices, nTrackSize,  geomManager, magMap, physicsTable,
     sbData, secTracks, stackSize, numStep, runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void elec_cpu(GXTrack *track,
	      size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int numStep, G4int runType)
{
  //process manager
  GXProcessManager aElectronProcessManager;

  //EM processes and models

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

  aElectronProcessManager.AddElectronProcess(&bremProcess,
					     &ioniProcess,
					     &mscProcess);

  GXSteppingManager aSteppingManager(geomManager,magMap);
  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GXTrackingManager aTrackingManager(&aElectronProcessManager,
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
    aTrackingManager.ProcessOneTrack(&track[tid]);
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
#include "GXProcessManager.cu"
#include "GXTrackingManager.cu"
#include "GXSteppingManager.cu"

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

#endif