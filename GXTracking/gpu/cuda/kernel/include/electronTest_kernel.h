#ifndef ELECTRONTEST_KERNEL_H
#define ELECTRONTEST_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPPhysics2DVector.h"

#include <curand.h>
#include <curand_kernel.h>

#include "dma_kernel.h"

// electron_kernel wrapper
void electron_gpu(curandState *devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table, 
		  GPPhysics2DVector* sbData,
		  GXTrack *secTrack,
		  G4int *stackSize,
		  G4int *offset,
		  G4int isStrack,
		  G4int runType,
		  int nBlocks,
		  int nThreads); 

void electron_cpu(GXTrack *trackIn, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table,
		  GPPhysicsTable* eIoni_table,
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData,
		  GXTrack *secTrack,
		  G4int *stackSize,
		  G4int isStrack,
		  G4int runType);

// brem_kernel wrapper
void brem_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* msc_table, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks,
	      G4int *stackSize,
	      G4int *offset,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads); 

void brem_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table,
	      GPPhysicsTable* eIoni_table,
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// brem_kernel wrapper
void brem_gpu_dma(curandState *devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* msc_table, 
		  GPPhysics2DVector* sbData,
		  GXSecContainer *secContainer_d,
		  G4int *stackSize,
		  G4int *offset,
		  int nBlocks,
		  int nThreads); 

void brem_cpu_dma(GXTrack *trackIn, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table,
		  GPPhysicsTable* eIoni_table,
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData);

// ioni_kernel wrapper
void ioni_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* eIoni_range, 
	      GPPhysicsTable* eIoni_dedx, 
	      GPPhysicsTable* eIoni_invr, 
	      GPPhysicsTable* msc_table, 
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks,
	      G4int *stackSize,
	      G4int *offset,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads); 

void ioni_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table,
	      GPPhysicsTable* eIoni_table,
	      GPPhysicsTable* eIoni_range, 
	      GPPhysicsTable* eIoni_dedx, 
	      GPPhysicsTable* eIoni_invr, 
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);


void ioni_gpu_dma(curandState *devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table, 
		  GPPhysics2DVector* sbData,
		  GXSecContainer *secContainer_d,
		  G4int *stackSize,
		  G4int *offset,
		  int nBlocks,
		  int nThreads); 

void ioni_cpu_dma(GXTrack *trackIn, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table,
		  GPPhysicsTable* eIoni_table,
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData);

// msc_kernel wrapper
void msc_gpu(curandState *devStates,
	     GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* eBrem_table, 
	     GPPhysicsTable* eIoni_table, 
	     GPPhysicsTable* msc_table, 
	     GPPhysics2DVector* sbData,
	     GXTrack *secTracks,
	     G4int *stackSize,
	     G4int *offset,
	     G4int isStrack,
	     G4int runType,
	     int nBlocks,
	     int nThreads); 

void msc_cpu(GXTrack *trackIn, size_t nTrackSize,
	     GPPhysicsTable* eBrem_table,
	     GPPhysicsTable* eIoni_table,
	     GPPhysicsTable* msc_table,
	     GPPhysics2DVector* sbData,
	     GXTrack *secTrack,
	     G4int *stackSize,
	     G4int isStrack,
	     G4int runType);

#endif
