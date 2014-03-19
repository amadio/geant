#include <iostream>
#include <memory>
#include <stdio.h>

//Geometry
#include "GPVPhysicalVolume.h"
#include "GPThreeVector.h"
#include "GPNavigator.h"
#include "GPTouchableHistory.h"
#include "GPMultiLevelLocator.h"

//Transportation
#include "GPFieldMap.h"
#include "GPTransportation.h"
#include "GPTransportationManager.h"

//EM Physics
#include "GPBrem.h"
#include "GPIoni.h"
#include "GPMsc.h"
#include "GPSTLVector.h"
#include "GPPhysicsTable.h"

#include "GPTrack.h"
#include "GPFieldTrack.h"


GLOBALFUNC 
void curand_setup_kernel(curandState * state, unsigned long seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &state[tid]);
}


GLOBALFUNC 
void tracking_kernel(GPVPhysicalVolume *world,
		     GPFieldMap *magMap,
		     curandState* devStates,
		     GPPhysicsTable* eBrem_table,
		     GPPhysicsTable* eIoni_table,
		     GPPhysicsTable* msc_table,
		     GPTrack *track, 
		     size_t nTrackSize) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMagneticField magField;    
  GPMagneticField_Constructor(&magField,magMap);

  GPEquationOfMotion equaOfMotion;
  GPEquationOfMotion_Constructor(&equaOfMotion,&magField);

  GPClassicalRK4 rk4;
  GPClassicalRK4_Constructor(&rk4,&equaOfMotion,6);

  GPMagInt_Driver magIntDriver;
  GPMagInt_Driver_Constructor(&magIntDriver,1.0,&rk4,6,0);

  GPChordFinder chordFinder;
  GPChordFinder_Constructor(&chordFinder,&magIntDriver);        

  GPFieldManager aFieldManager;
  GPFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GPMultiLevelLocator mLocator;
  GPMultiLevelLocator_Constructor(&mLocator, &aNavigator);

  //Propagator
  GPPropagatorInField propagatorInField;
  GPPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                  &aFieldManager,&mLocator); 

  //Transporation
  GPTransportation transp;
  GPTransportation_Constructor(&transp,&propagatorInField,0);

  //EM Physics

  GPBrem brem;
  GPIoni ioni;
  GPMsc msc;

  bool useIntegral = true;
  bool useLambdaTable = true;
  G4double energyLoss = 0;
  GPForceCondition condition;

  while(tid < nTrackSize) {

    //--->EM Physics

    brem.threadId = tid;
    brem.SetCurandStates(devStates);
    brem.UseIntegral(useIntegral);
    brem.SetLambdaTable(eBrem_table);
    brem.UseLambdaTable(useLambdaTable);
    track[tid].length = brem.PostStepGetPhysicalInteractionLength(track[tid].E, track[tid].length, &condition);

    ioni.threadId = tid;
    ioni.SetCurandStates(devStates);
    ioni.UseIntegral(useIntegral);
    ioni.SetLambdaTable(eIoni_table);
    ioni.UseLambdaTable(useLambdaTable);
    track[tid].length = ioni.PostStepGetPhysicalInteractionLength(track[tid].E, track[tid].length, &condition);

    msc.threadId = tid;
    msc.SetCurandStates(devStates);
    msc.UseIntegral(useIntegral);
    msc.SetLambdaTable(msc_table);
    msc.UseLambdaTable(useLambdaTable);
    track[tid].length = msc.PostStepGetPhysicalInteractionLength(track[tid].E, track[tid].length, &condition);

    energyLoss += brem.PostStepDoIt(&track[tid]);
    //    printf("eloss from device = %f\n",energyLoss);
    energyLoss += ioni.PostStepDoIt(&track[tid]);
    energyLoss += msc.PostStepDoIt(&track[tid]);

    //Transportation
    double proposedStep = 0.0;
    double safety = 0.0;
    GPGPILSelection  fGPILSelection;

    proposedStep = track[tid].step;

    track[tid].length = 
    GPTransportation_AlongStepGetPhysicalInteractionLength(&transp,
							   &track[tid],
							   1.0,
							   proposedStep,
							   safety,
							   &fGPILSelection);
    tid += blockDim.x * gridDim.x;

  }

  //synchronize threads in this block
  //  __syncthreads();

}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_gpu( GPVPhysicalVolume *world,	    
		   GPFieldMap *magMap,
		   GPPhysicsTable* eBrem_table,
		   GPPhysicsTable* eIoni_table,
		   GPPhysicsTable* msc_table,
		   GPTrack *track,
		   size_t nTrackSize,
		   int NBLOCKS,
		   int NTHREADS) 
{
  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 
  //imin(32, (nTrackSize+threadsPerBlock-1)/threadsPerBlock);

  curandState* devStates = 0;

  cudaMalloc(&devStates, blocksPerGrid * blocksPerGrid * sizeof(curandState));
  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,
							     time(NULL));

  tracking_kernel<<< blocksPerGrid, threadsPerBlock >>>(world,magMap,
							devStates,
							eBrem_table,
							eIoni_table,
							msc_table,
							track,
							nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "cmsCudaTransportation status = " << kstatus << "\n";

  //  cudaFree(devStates);

}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_cpu( GPVPhysicalVolume *world,
		   GPFieldMap *magMap, 
		   GPPhysicsTable* eBrem_table,
		   GPPhysicsTable* eIoni_table,
		   GPPhysicsTable* msc_table,
		   GPTrack *track, 
		   size_t nTrackSize)
{
  GPMagneticField magField;    
  GPMagneticField_Constructor(&magField,magMap);

  GPEquationOfMotion equaOfMotion;
  GPEquationOfMotion_Constructor(&equaOfMotion,&magField);

  GPClassicalRK4 rk4;
  GPClassicalRK4_Constructor(&rk4,&equaOfMotion,6);

  GPMagInt_Driver magIntDriver;
  GPMagInt_Driver_Constructor(&magIntDriver,1.0,&rk4,6,0);

  GPChordFinder chordFinder;
  GPChordFinder_Constructor(&chordFinder,&magIntDriver);        

  GPFieldManager aFieldManager;
  GPFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GPMultiLevelLocator mLocator;
  GPMultiLevelLocator_Constructor(&mLocator, &aNavigator);

  //Propagator
  GPPropagatorInField propagatorInField;
  GPPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                  &aFieldManager,&mLocator); 

  //Transporation
  GPTransportation transp;
  GPTransportation_Constructor(&transp,&propagatorInField,0);

  //EM Physics
  GPBrem brem;
  GPIoni ioni;
  GPMsc msc;

  bool useIntegral = true;
  bool useLambdaTable = true;
  G4double energyLoss = 0;
  GPForceCondition condition;

  srand(time(NULL));

  for (int i = 0 ; i < nTrackSize ; i++) {

    //--->EM Physics
    brem.threadId = i;
    brem.UseIntegral(useIntegral);
    brem.SetLambdaTable(eBrem_table);
    brem.UseLambdaTable(useLambdaTable);
    track[i].length = brem.PostStepGetPhysicalInteractionLength(track[i].E, track[i].length, &condition);
//    printf("%dth track length: %f\n",i,track[i].length);

    ioni.threadId = i;
    ioni.UseIntegral(useIntegral);
    ioni.SetLambdaTable(eIoni_table);
    ioni.UseLambdaTable(useLambdaTable);
    track[i].length = ioni.PostStepGetPhysicalInteractionLength(track[i].E, track[i].length, &condition);

    msc.threadId = i;
    msc.UseIntegral(useIntegral);
    msc.SetLambdaTable(msc_table);
    msc.UseLambdaTable(useLambdaTable);
    track[i].length = msc.PostStepGetPhysicalInteractionLength(track[i].E, track[i].length, &condition);

    energyLoss += brem.PostStepDoIt(&track[i]);
    //    printf("eloss from host = %f\n",energyLoss);
    energyLoss += ioni.PostStepDoIt(&track[i]);
    energyLoss += msc.PostStepDoIt(&track[i]);

    //Transportation
    double proposedStep = 0.0;
    double safety = 0.0;
    GPGPILSelection  fGPILSelection;
    proposedStep = track[i].step;

    track[i].length = 
    GPTransportation_AlongStepGetPhysicalInteractionLength(&transp,
							   &track[i],
							   1.0,
							   proposedStep,
							   safety,
							   &fGPILSelection);
  }
}
