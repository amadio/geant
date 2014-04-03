#ifndef EMPHYSICS_KERNEL_H
#define EMPHYSICS_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"

#include <curand.h>
#include <curand_kernel.h>

// EMPhysics_kernel wrapper
void EMPhysics_gpu(curandState *devStates,
		   GXTrack *track, size_t nTrackSize,
		   GPPhysicsTable* eBrem_table, 
		   GPPhysicsTable* eIoni_table, 
		   GPPhysicsTable* msc_table, 
		   bool useIntegral, bool useLambdaTable,
		   int nBlocks,
		   int nThreads); 

void EMPhysics_cpu(GXTrack *trackIn, size_t nTrackSize,
		   GPPhysicsTable* eBrem_table,
		   GPPhysicsTable* eIoni_table,
		   GPPhysicsTable* msc_table,
		   bool useIntegral, bool useLambdaTable);

// bream_kernel wrapper
void brem_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      bool useIntegral, bool useLambdaTable,
	      int nBlocks,
	      int nThreads); 

void brem_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table,
	      bool useIntegral, bool useLambdaTable);

// ioni_kernel wrapper
void ioni_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eIoni_table, 
	      bool useIntegral, bool useLambdaTable,
	      int nBlocks,
	      int nThreads); 

void ioni_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPPhysicsTable* eIoni_table,
	      bool useIntegral, bool useLambdaTable);

// msc_kernel wrapper
void msc_gpu(curandState *devStates,
	     GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* eBrem_table, 
	     bool useIntegral, bool useLambdaTable,
	     int nBlocks,
	     int nThreads); 

void msc_cpu(GXTrack *trackIn, size_t nTrackSize,
	     GPPhysicsTable* msc_table,
	     bool useIntegral, bool useLambdaTable);

#endif
