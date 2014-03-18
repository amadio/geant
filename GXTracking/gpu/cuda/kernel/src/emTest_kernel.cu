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
#include "GPMultiLevelLocator.h"

#include "stdio.h"

//-----------------------------------------------------------------------------
// Electron Kernel
//-----------------------------------------------------------------------------

GLOBALFUNC
void electron_kernel(curandState* devStates,
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

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  while (tid < nTrackSize) {
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for electron_kernel
//-----------------------------------------------------------------------------
void electron_gpu(curandState* devStates,
                  GXTrack *track, size_t nTrackSize,
                  GPGeomManager *geomManager,
                  GXFieldMap *magMap,
		  GPPhysicsTable* physicsTable, 
                  GPPhysics2DVector* sbData,
                  GXTrack *secTracks, G4int *stackSize, G4int *offset,
                  G4int numStep, G4int runType, 
                  int blocksPerGrid, int threadsPerBlock) 
{
  electron_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		  nTrackSize,  geomManager, magMap, physicsTable,
                  sbData, secTracks, stackSize, offset, numStep,runType);
}
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void electron_cpu(GXTrack *track, size_t nTrackSize,
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

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aElectronProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  for (int tid = 0; tid < nTrackSize; tid++) {
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,0);
    aTrackingManager.ProcessOneTrack(&track[tid]);
  }
}

//-----------------------------------------------------------------------------
// Photon Kernel
//-----------------------------------------------------------------------------

GLOBALFUNC
void photon_kernel(curandState* devStates,
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

  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  while (tid < nTrackSize) {
    aSteppingManager.SetSecondaryStack(secTracks,stackSize,offset);
    aTrackingManager.ProcessOneTrack(&track[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void photon_gpu(curandState* devStates,
		GXTrack *track, size_t nTrackSize,
                GPGeomManager *geomManager,
                GXFieldMap *magMap,
		GPPhysicsTable* physicsTable, 
		GXTrack *secTracks, G4int *stackSize, G4int *offset,
		G4int numStep, G4int runType,
		int blocksPerGrid, int threadsPerBlock) 
{
  photon_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize, geomManager, magMap, physicsTable,
		       secTracks, stackSize, offset, numStep, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void photon_cpu(GXTrack *track, size_t nTrackSize,
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

  aPhotonProcessManager.AddPhotonProcess(comptProcess);  
  aPhotonProcessManager.AddPhotonProcess(convProcess);  
  aPhotonProcessManager.AddPhotonProcess(photProcess);  

  GPSteppingManager aSteppingManager(geomManager,magMap);

  GPStep aStep;
  aSteppingManager.SetStep(&aStep);

  GPTrackingManager aTrackingManager(&aPhotonProcessManager,
				     &aSteppingManager);
  aTrackingManager.SetNumberOfSteps(numStep);

  for (int tid = 0; tid < nTrackSize; ++tid) {
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
#include "GPMultiLevelLocator.cu"
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
#include "GPFieldTrack.cu"
#include "GPFieldManager.cu"
#include "GPMagneticField.cu"

#include "GPClassicalRK4.cu"
#include "GPEquationOfMotion.cu"
#include "GPMagIntegratorDriver.cu"
#include "GPChordFinder.cu"
#include "GPPropagatorInField.cu"
#include "GPTransportation.cu"
#include "GPTransportationManager.cu"
