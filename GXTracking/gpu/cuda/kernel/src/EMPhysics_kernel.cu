#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

#include "GXBrem.h"
#include "GXIoni.h"
#include "GXMsc.h"

#include "GPSTLVector.h"

GLOBALFUNC
void EMPhysics_kernel(curandState* devStates, GXTrack *track, 
		      size_t nTrackSize,
		      GPPhysicsTable* eBrem_table, 
		      GPPhysicsTable* eIoni_table, 
		      GPPhysicsTable* msc_table,
		      bool useIntegral, bool useLambdaTable) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXBrem brem;
  brem.threadId = tid;
  brem.SetCurandStates(devStates);
  brem.UseIntegral(useIntegral);
  brem.SetLambdaTable(eBrem_table);
  brem.UseLambdaTable(useLambdaTable);

  GXIoni ioni;
  ioni.threadId = tid;
  ioni.SetCurandStates(devStates);
  ioni.UseIntegral(useIntegral);
  ioni.SetLambdaTable(eIoni_table);
  ioni.UseLambdaTable(useLambdaTable);

  GXMsc msc;
  msc.threadId = tid;
  msc.SetCurandStates(devStates);
  msc.UseIntegral(useIntegral);
  msc.SetLambdaTable(msc_table);
  msc.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_brem = 0.0;
  G4double step_ioni = 0.0;
  G4double step_msc  = 0.0;
  G4double proposedStep = 0.0;
  
  while (tid < nTrackSize) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_brem = brem.PostStepGetPhysicalInteractionLength(E, step, &condition);
    step_ioni = ioni.PostStepGetPhysicalInteractionLength(E, step, &condition);
    step_msc = msc.PostStepGetPhysicalInteractionLength(E, step, &condition);

    //physics model defining the current step
    unsigned int model =0;
    if(step_brem > step_ioni) model = 1;   
    if(step_ioni > step_msc) model = 2;  

    switch (model) {
    case 0 : 
      proposedStep = step_brem;
      energyLoss += brem.PostStepDoIt(&track[tid]);
      break;
    case 1 :
      proposedStep = step_ioni;
      energyLoss += ioni.PostStepDoIt(&track[tid]);
      break;
    case 2 :
      proposedStep = step_msc;			       
      energyLoss += msc.PostStepDoIt(&track[tid]);
      break;
    default :
      break;
    }

    //alongstep
    step_msc = msc.AlongStepGetPhysicalInteractionLength(&track[tid], 
							 step, &selection);
    energyLoss += msc.AlongStepDoIt(&track[tid]);
    if(step_msc < proposedStep) proposedStep = step_msc;

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = proposedStep;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void EMPhysics_gpu(curandState* devStates, GXTrack *track, size_t nTrackSize,
		   GPPhysicsTable* eBrem_table, 
		   GPPhysicsTable* eIoni_table, 
		   GPPhysicsTable* msc_table,
		   bool useIntegral, bool useLambdaTable, 
		   int blocksPerGrid, int threadsPerBlock) 
{
  EMPhysics_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
    nTrackSize,eBrem_table,eIoni_table,msc_table,useIntegral,useLambdaTable);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void EMPhysics_cpu(GXTrack *track, size_t nTrackSize,
		   GPPhysicsTable* eBrem_table, 
		   GPPhysicsTable* eIoni_table, 
		   GPPhysicsTable* msc_table,
		   bool useIntegral, bool useLambdaTable) 
{
  GXBrem brem;
  brem.UseIntegral(useIntegral);
  brem.SetLambdaTable(eBrem_table);
  brem.UseLambdaTable(useLambdaTable);

  GXIoni ioni;
  ioni.UseIntegral(useIntegral);
  ioni.SetLambdaTable(eIoni_table);
  ioni.UseLambdaTable(useLambdaTable);

  GXMsc msc;
  msc.UseIntegral(useIntegral);
  msc.SetLambdaTable(msc_table);
  msc.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_brem = 0.0;
  G4double step_ioni = 0.0;
  G4double step_msc  = 0.0;
  G4double proposedStep = 0.0;

  for (int tid = 0; tid < nTrackSize; tid++) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;

    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_brem = brem.PostStepGetPhysicalInteractionLength(E, step, &condition);
    step_ioni = ioni.PostStepGetPhysicalInteractionLength(E, step, &condition);
    step_msc = msc.PostStepGetPhysicalInteractionLength(E, step, &condition);

    //physics model defining the current step
    unsigned int model =0;
    if(step_brem > step_ioni) model = 1;   
    if(step_ioni > step_msc) model = 2;  

    switch (model) {
    case 0 : 
      proposedStep = step_brem;
      energyLoss += brem.PostStepDoIt(&track[tid]);
      break;
    case 1 :
      proposedStep = step_ioni;
      energyLoss += ioni.PostStepDoIt(&track[tid]);
      break;
    case 2 :
      proposedStep = step_msc;			       
      energyLoss += msc.PostStepDoIt(&track[tid]);
      break;
    default :
      break;
    }

    //alongstep
    step_msc = msc.AlongStepGetPhysicalInteractionLength(&track[tid],
							 step, &selection);
    energyLoss += msc.AlongStepDoIt(&track[tid]);
    if(step_msc < proposedStep) proposedStep = step_msc;

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = proposedStep;

  }
}

//Bremsstralung model only
GLOBALFUNC
void brem_kernel(curandState* devStates, GXTrack *track, size_t nTrackSize,
		 GPPhysicsTable* eBrem_table, 
		 bool useIntegral, bool useLambdaTable) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  GXBrem brem;
  brem.threadId = tid;
  brem.SetCurandStates(devStates);
  brem.UseIntegral(useIntegral);
  brem.SetLambdaTable(eBrem_table);
  brem.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_brem = 0.0;

  while (tid < nTrackSize) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_brem = brem.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss -= brem.PostStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss;
    track[tid].s = step_brem;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void brem_gpu(curandState* devStates, GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      bool useIntegral, bool useLambdaTable, 
	      int blocksPerGrid, int threadsPerBlock)
{
  brem_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates,track,nTrackSize,eBrem_table,useIntegral,useLambdaTable);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------
void brem_cpu(GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table,
	      bool useIntegral, bool useLambdaTable) 
{
  GXBrem brem;
  brem.UseIntegral(useIntegral);
  brem.SetLambdaTable(eBrem_table);
  brem.UseLambdaTable(useLambdaTable);
  
  GPForceCondition condition;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_brem = 0.0;

  for (int tid = 0; tid < nTrackSize; tid++) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_brem = brem.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss -= brem.PostStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = step_brem;
  }
}

//Ionization model
GLOBALFUNC
void ioni_kernel(curandState* devStates, GXTrack *track, size_t nTrackSize,
		 GPPhysicsTable* eIoni_table, 
		 bool useIntegral, bool useLambdaTable)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXIoni ioni;
  ioni.threadId = tid;
  ioni.SetCurandStates(devStates);
  ioni.UseIntegral(useIntegral);
  ioni.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_ioni = 0.0;

  while (tid < nTrackSize) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_ioni = ioni.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss -= ioni.PostStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = step_ioni;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void ioni_gpu(curandState* devStates, GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eIoni_table, 
	      bool useIntegral, bool useLambdaTable,
	      int blocksPerGrid, int threadsPerBlock)
{ 
  ioni_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates,track,nTrackSize,eIoni_table,useIntegral,useLambdaTable);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------
void ioni_cpu(GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eIoni_table, 
	      bool useIntegral, bool useLambdaTable) 
{
  GXIoni ioni;
  ioni.UseIntegral(useIntegral);
  ioni.SetLambdaTable(eIoni_table);
  ioni.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_ioni = 0.0;

  for (int tid = 0; tid < nTrackSize; tid++) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_ioni = ioni.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss -= ioni.PostStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss;
    track[tid].s = step_ioni;

  }
}

// Multiple Scattering Only
GLOBALFUNC
void msc_kernel(curandState* devStates, GXTrack *track, size_t nTrackSize,
		GPPhysicsTable* msc_table,
		bool useIntegral, bool useLambdaTable)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  GXMsc msc;
  msc.threadId = tid;
  msc.SetCurandStates(devStates);
  msc.UseIntegral(useIntegral);
  msc.SetLambdaTable(msc_table);
  msc.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_msc  = 0.0;

  while (tid < nTrackSize) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_msc = msc.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss += msc.PostStepDoIt(&track[tid]);

    step_msc = msc.AlongStepGetPhysicalInteractionLength(&track[tid],
							 step, &selection);
    energyLoss += msc.AlongStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = step_msc;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void msc_gpu(curandState* devStates, GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* msc_table,
	     bool useIntegral, bool useLambdaTable,
	      int blocksPerGrid, int threadsPerBlock)
{
  msc_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates,track,nTrackSize,msc_table,useIntegral,useLambdaTable);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------
void msc_cpu(GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* msc_table,
	     bool useIntegral, bool useLambdaTable)
{
  GXMsc msc;
  msc.UseIntegral(useIntegral);
  msc.SetLambdaTable(msc_table);
  msc.UseLambdaTable(useLambdaTable);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_msc  = 0.0;

  for (int tid = 0; tid < nTrackSize; tid++) {

    px   = track[tid].px  ;
    py   = track[tid].py  ;
    pz   = track[tid].pz  ;
    step = track[tid].s   ;
    
    E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

    step_msc = msc.PostStepGetPhysicalInteractionLength(E, step, &condition);
    energyLoss += msc.PostStepDoIt(&track[tid]);

    step_msc = msc.AlongStepGetPhysicalInteractionLength(&track[tid], 
							 step, &selection);
    energyLoss += msc.AlongStepDoIt(&track[tid]);

    //use q as a temporary store
    track[tid].q = energyLoss ;
    track[tid].s = step_msc;

  }
}

#include "random_kernel.cu"

#include "GPThreeVector.cu"
#include "GPPhysicsTable.cu"
#include "GXBrem.cu"
#include "GXIoni.cu"
#include "GXMsc.cu"
