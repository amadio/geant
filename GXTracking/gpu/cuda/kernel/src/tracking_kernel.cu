#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>

//Common

#include "GPThreeVector.h"
#include "GPThreeVectorList.h"
#include "GPRotationMatrix.h"
#include "GPUtils.h"

//Material
#include "GPMaterial.h"

//Geometry
#include "GPGeomManager.h"
#include "GPVPhysicalVolume.h"
#include "GPTouchableHistory.h"
#include "GPNavigator.h"
#include "GXMultiLevelLocator.h"

//Transportation
#include "GXFieldMap.h"
#include "GXTrack.h"

#include "GXFieldTrack.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"
#include "GXMagIntegratorDriver.h"
#include "GXChordFinder.h"
#include "GXFieldManager.h"
#include "GXPropagatorInField.h"
#include "GXTransportation.h"

#include "GXBrem.h"
#include "GXIoni.h"
#include "GXMsc.h"

#include "GPSTLVector.h"

#include "stdio.h"

#include "util_kernel.h"

//-----------------------------------------------------------------------------
//  Transportation engine
//-----------------------------------------------------------------------------
#ifdef __CUDACC__
GLOBALFUNC
void tracking_kernel(curandState* devStates,
                     GPGeomManager *geomManager,
		     GXFieldMap *magMap,
		     GXTrack *track, 
		     GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		     size_t nTrackSize) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  while(tid < nTrackSize) {
    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;

    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,&condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;
    
    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGetPhysicalInteractionLength(&transp,
									  &atrack,
									  1.0,
									  proposedStep,
									  safety,
									  &fGPILSelection);
    tid += blockDim.x * gridDim.x;
    
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_gpu(curandState* devStates, 
                  GPGeomManager *geomManager,
		  GXFieldMap *magMap,
		  GXTrack *track,  
		  GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		  size_t nTrackSize,
		  int NBLOCKS,
		  int NTHREADS,
		  cudaStream_t stream)
{
  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 
  //  int blocksPerGrid   = (int)((nTrackSize+threadsPerBlock-1)/threadsPerBlock);
  
  tracking_kernel<<< blocksPerGrid, threadsPerBlock, 0 , stream >>>(devStates,
								    geomManager,magMap,
								    track,eBrem_table,eIoni_table,msc_table,
								    nTrackSize);
}
#endif

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_cpu(GPGeomManager *geomManager,
		  GXFieldMap *magMap, 
		  GXTrack *track, 
		  GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		  size_t nTrackSize)
{
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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;
   
    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
					      &condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;

    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGetPhysicalInteractionLength(&transp,
									  &atrack,
									  1.0,
									  proposedStep,
									  safety,
									  &fGPILSelection);
  }

}

//-----------------------------------------------------------------------------
//  Transportation engine - Electron/Positron
//-----------------------------------------------------------------------------
#ifdef __CUDACC__
GLOBALFUNC
void tracking_electron_kernel(curandState* devStates,
                              GPGeomManager *geomManager,
			      GXFieldMap *magMap,
			      GXTrack *track, 
			      GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			      size_t nTrackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  while(tid < nTrackSize) {

    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;
    
    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
					      &condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;
    
    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGPIL_Electron(&transp,
							   &atrack,
							   1.0,
							   proposedStep,
							   safety,
							   &fGPILSelection);
    tid += blockDim.x * gridDim.x;
    
  }
  
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_electron_gpu(curandState* devStates,
                           GPGeomManager *geomManager,
			   GXFieldMap *magMap,
			   GXTrack *track,  
			   GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			   size_t nTrackSize,
			   int NBLOCKS,
			   int NTHREADS, 
			   cudaStream_t stream)
{
  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 
  

  tracking_electron_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
     (devStates,geomManager,magMap,track,eBrem_table,eIoni_table,msc_table,nTrackSize);
  
}
#endif

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_electron_cpu( GPGeomManager *geomManager,
			    GXFieldMap *magMap, 
			    GXTrack *track, 
			    GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			    size_t nTrackSize)
{
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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;
    
    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,&condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;

    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGPIL_Electron(&transp,
							   &atrack,
							   1.0,
							   proposedStep,
							   safety,
							   &fGPILSelection);
  }

}

//-----------------------------------------------------------------------------
//  Transportation engine - Photon
//-----------------------------------------------------------------------------

#ifdef __CUDACC__
GLOBALFUNC
void tracking_photon_kernel(curandState* devStates,
                            GPGeomManager *geomManager,
			    GXFieldMap *magMap,
			    GXTrack *track, 
			    GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			    size_t nTrackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  while(tid < nTrackSize) {

    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;

    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
					      &condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;
    
    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGPIL_Photon(&transp,
							 &atrack,
							 1.0,
							 proposedStep,
							 safety,
							 &fGPILSelection);
    tid += blockDim.x * gridDim.x;
    
  }
  
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_photon_gpu(curandState* devStates, 
                         GPGeomManager *geomManager,
			 GXFieldMap *magMap,
			 GXTrack *track,  
			 GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			 size_t nTrackSize,
			 int NBLOCKS,
			 int NTHREADS,
			 cudaStream_t stream)
{
  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 
  
  tracking_photon_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (devStates,geomManager,magMap,track,eBrem_table,eIoni_table,msc_table,nTrackSize);
}
#endif

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_photon_cpu(GPGeomManager *geomManager,
			 GXFieldMap *magMap, 
			 GXTrack *track, 
			 GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
			 size_t nTrackSize)
{
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

  //EMPhysics
  GXBrem brem;
  GXIoni ioni;
  GXMsc msc;
  GPForceCondition condition;
  GPGPILSelection  selection;

  EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
		 eBrem_table,eIoni_table,msc_table);

  G4double x,y,z,step;
  G4double proposedStep = 0.0;

  GXTrack atrack;

  //  #pragma omp parallel for
  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    step = atrack.s   ;

    //Geometry - Initialize Navigator and construct related structures
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
					  NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
    
    GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
				    &aFieldManager,&mLocator); 
    
    GXTransportation_Constructor(&transp,&propagatorInField,0);
    
    //EM Physics
    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
					      &condition,&selection);

    //Transportation
    proposedStep += step;
    GPGPILSelection  fGPILSelection;
    G4double safety = 0.0;

    //add energyLoss here as faking output
    track[tid].s = GXTransportation_AlongStepGPIL_Photon(&transp,
							 &atrack,
							 1.0,
							 proposedStep,
							 safety,
							 &fGPILSelection);
  }

}

//-----------------------------------------------------------------------------
// Other Kernels
//-----------------------------------------------------------------------------
#include "random_kernel.cu"
#include "util_kernel.cu"

//-----------------------------------------------------------------------------
// GPCommon
//-----------------------------------------------------------------------------

#include "GPThreeVector.cu"
#include "GPThreeVectorList.cu"
#include "GPRotationMatrix.cu"
#include "GPUtils.cu"

//-----------------------------------------------------------------------------
// GPMaterial
//-----------------------------------------------------------------------------
#include "GPElement.cu"
#include "GPMaterial.cu"

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

#include "GXMultiLevelLocator.cu"

#include "GPBox.cu"
#include "GPCons.cu"
#include "GPOrb.cu"
#include "GPTrd.cu"
#include "GPTubs.cu"
#include "GPVSolid.cu"

#include "GPChargeState.cu"

#include "GPUserGeometry.cu"
#include "GPVoxelHeader.cu"

//chordFiner
#include "GXFieldTrack.cu"
#include "GXClassicalRK4.cu"
#include "GXEquationOfMotion.cu"
#include "GXMagneticField.cu"
#include "GXMagIntegratorDriver.cu"
#include "GXChordFinder.cu"
#include "GXFieldManager.cu"
#include "GXPropagatorInField.cu"
#include "GXTransportation.cu"

//EMPhysics
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"
#include "GXBrem.cu"
#include "GXIoni.cu"
#include "GXMsc.cu"

#ifdef __CUDACC__
GLOBALFUNC
void random_testing_kernel(curandState* devStates, double *result)
{
   unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
 
   result[tid] = rand_wrapper(&devStates[tid]);
}

void random_testing_gpu(curandState* devStates, double *results, 
                        int  blocksPerGrid, int threadsPerBlock,
                        cudaStream_t stream)
{

   random_testing_kernel<<< blocksPerGrid, threadsPerBlock, 0 , stream >>>(devStates,results);
  
}
#endif

