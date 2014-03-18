#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>
#include <cstdlib>

#include "GPBrem.h"
#include "GPIoni.h"

#include "GPSTLVector.h"

GLOBALFUNC void curand_setup_kernel(curandState * state, unsigned long seed);
/*
GLOBALFUNC void curand_setup_kernel(curandState * state, unsigned long seed) {
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	curand_init(seed, tid, 0, &state[tid]);
}
*/

GLOBALFUNC
void EMPhysics_kernel(curandState* devStates, GPTrack *track, size_t nTrackSize,
		       GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		       bool useIntegral, bool useLambdaTable) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  GPBrem brem;
  GPIoni ioni;
  GPForceCondition condition;

  while (tid < nTrackSize) {
    brem.threadId = tid;
    brem.SetCurandStates(devStates);
    brem.UseIntegral(useIntegral);
    brem.SetLambdaTable(eBrem_table);
    brem.UseLambdaTable(useLambdaTable);
    track[tid].length = brem.PostStepGetPhysicalInteractionLength(track[tid].E, track[tid].length, &condition);
    track[tid].E -= brem.PostStepDoIt(&track[tid]);

    ioni.threadId = tid;
    ioni.SetCurandStates(devStates);
    ioni.UseIntegral(useIntegral);
    ioni.SetLambdaTable(eIoni_table);
    ioni.UseLambdaTable(useLambdaTable);
    track[tid].length = ioni.PostStepGetPhysicalInteractionLength(track[tid].E, track[tid].length, &condition);
    track[tid].E -= ioni.PostStepDoIt(&track[tid]);

    tid += blockDim.x * gridDim.x;
  }

  //synchronize threads in this block
  __syncthreads();
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void EMPhysics_gpu(GPTrack *track, size_t nTrackSize,
		    GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		    bool useIntegral, bool useLambdaTable, int NBLOCKS, int NTHREADS) {
	int kstatus = 0;

	int threadsPerBlock = NTHREADS;
	int blocksPerGrid = NBLOCKS;

	curandState* devStates = 0;
	cudaMalloc(&devStates, blocksPerGrid * blocksPerGrid * sizeof(curandState));
	curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,time(NULL));

	EMPhysics_kernel<<< blocksPerGrid, threadsPerBlock >>>
	  (devStates,track,nTrackSize,eBrem_table,eIoni_table,msc_table,useIntegral,useLambdaTable);

	kstatus = cudaThreadSynchronize();
	if (kstatus)
		std::cout << "cmsCudaTransportation status = " << kstatus << "\n";

	cudaFree(devStates);

}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void EMPhysics_cpu(GPTrack *track, size_t nTrackSize,
		    GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
		    bool useIntegral, bool useLambdaTable) {
  GPBrem brem;
  GPIoni ioni;
  GPForceCondition condition;

  for (int i = 0; i < nTrackSize; i++) {
    brem.threadId = i;
    brem.UseIntegral(useIntegral);
    brem.SetLambdaTable(eBrem_table);
    brem.UseLambdaTable(useLambdaTable);
    track[i].length = brem.PostStepGetPhysicalInteractionLength(track[i].E, track[i].length, &condition);
    track[i].E -= brem.PostStepDoIt(&track[i]);

    ioni.threadId = i;
    ioni.UseIntegral(useIntegral);
    ioni.SetLambdaTable(eIoni_table);
    ioni.UseLambdaTable(useLambdaTable);
    track[i].length = ioni.PostStepGetPhysicalInteractionLength(track[i].E, track[i].length, &condition);
    track[i].E -= ioni.PostStepDoIt(&track[i]);
  }

}

