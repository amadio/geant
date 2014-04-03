#ifndef PHOTONTEST_KERNEL_H
#define PHOTONTEST_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"

#include <curand.h>
#include <curand_kernel.h>

//the maximum number of secondaries per step
const G4int maxSecondaryPerStep = 2;

// compt_kernel wrapper
void compt_gpu(curandState *devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table, 
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int *offset,
	       G4int isStrack,
	       G4int runType,
	       int nBlocks,
	       int nThreads); 

void compt_cpu(GXTrack *trackIn, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table,
	       GPPhysicsTable* gConv_table,
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int isStrack,
	       G4int runType);

// conv_kernel wrapper
void conv_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table, 
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int *offset,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads); 

void conv_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table,
	      GPPhysicsTable* gConv_table,
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// phot_kernel wrapper
void phot_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table, 
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int *offset,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads); 

void phot_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table,
	      GPPhysicsTable* gConv_table,
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// gamma_kernel wrapper
void gamma_gpu(curandState *devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table, 
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int *offset,
	       G4int isStrack,
	       G4int runType,
	       int nBlocks,
	       int nThreads); 

void gamma_cpu(GXTrack *trackIn, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table,
	       GPPhysicsTable* gConv_table,
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int isStrack,
	       G4int runType);

#endif
