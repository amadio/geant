#include "stdio.h"
#include <stdlib.h>

#include "GPEmProcessType.h"
#include "GXTrack.h"

//-----------------------------------------------------------------------------
// count tracks by the physics process 
//-----------------------------------------------------------------------------

__global__
void count_by_process_kernel(G4int nTracks, GXTrack *tracks, 
			     G4int *nbrem, G4int *nioni)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < nTracks) {
    //offset is a global counter for the last array position of secondaries 
    if(tracks[tid].proc == 0 ) {
      atomicAdd(nbrem,1);
    }    
    if(tracks[tid].proc == 1 ) {
      atomicAdd(nioni,1);
    }    
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// wrapper for count
//-----------------------------------------------------------------------------
void count_by_process_gpu(G4int nTracks, GXTrack *tracks, 
			 G4int *nbrem, G4int *nioni, 
			 int blocksPerGrid, int threadsPerBlock,
                          cudaStream_t stream)
{
   count_by_process_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (nTracks, tracks, nbrem, nioni);
}

void count_by_process_cpu(G4int nTracks, GXTrack *tracks,
			  G4int *nbrem, G4int *nioni)
{
  //generate secondaries
  for(int tid = 0 ; tid < nTracks ; tid++) {
    //offset is a global counter for the last array position of secondaries 
    if(tracks[tid].proc == 0 ) {
      ++(*nbrem);
    }
    if(tracks[tid].proc == 1 ) {
      ++(*nioni);
    }
  }
}

//-----------------------------------------------------------------------------
// atomic add per block basis (by Azamat Mametjano - ANL)
//-----------------------------------------------------------------------------

__device__
void warpReduce(volatile int *sdata, unsigned int tid) {
  // smallest exec group is a warp of 32 threads
  // no need to sync threads within a warp
  sdata[tid] += sdata[tid + 32];
  sdata[tid] += sdata[tid + 16];
  sdata[tid] += sdata[tid +  8];
  sdata[tid] += sdata[tid +  4];
  sdata[tid] += sdata[tid +  2];
  sdata[tid] += sdata[tid +  1];
}

__global__
void count_by_process_block_kernel(G4int nTracks, GXTrack *tracks, 
				   G4int *nbrem, G4int *nioni)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  // declare and initialize
  extern __shared__ int sb[]; //shared memory nbrem counter
  int *si = &sb[blockDim.x];  //shared memory nioni counter
  sb[threadIdx.x] = 0;
  si[threadIdx.x] = 0;
  __syncthreads();

  // accumulate per-thread counters within a block's shared memory
  while(tid < nTracks) {
    if(tracks[tid].proc == 0 ) {
      sb[threadIdx.x] += 1;
    }    
    if(tracks[tid].proc == 1 ) {
      si[threadIdx.x] += 1;
    }    
    tid += blockDim.x * gridDim.x;
  }
  __syncthreads();

  // tree-wise reduce within a block, global memory atomicAdd across blocks
  // max threads per block is 1024 in any arch including Kepler (sm_35)
  if (blockDim.x > 512 && threadIdx.x < 512 && threadIdx.x+512 < blockDim.x) {
    sb[threadIdx.x]+=sb[threadIdx.x+512];
    si[threadIdx.x]+=si[threadIdx.x+512];
    __syncthreads();
  }
  if (blockDim.x > 256 && threadIdx.x < 256 && threadIdx.x+256 < blockDim.x) {
    sb[threadIdx.x]+=sb[threadIdx.x+256];
    si[threadIdx.x]+=si[threadIdx.x+256];
    __syncthreads();
  }
  if (blockDim.x > 128 && threadIdx.x < 128 && threadIdx.x+128 < blockDim.x) {
    sb[threadIdx.x]+=sb[threadIdx.x+128];
    si[threadIdx.x]+=si[threadIdx.x+128];
    __syncthreads();
  }
  if (blockDim.x >  64 && threadIdx.x <  64 && threadIdx.x+ 64 < blockDim.x) {
    sb[threadIdx.x]+=sb[threadIdx.x+ 64];
    si[threadIdx.x]+=si[threadIdx.x+ 64];
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    warpReduce(sb, threadIdx.x);
    warpReduce(si, threadIdx.x);
  }
  if (threadIdx.x == 0) {
    atomicAdd(nbrem, sb[0]);
    atomicAdd(nioni, si[0]);
  }
}

//-----------------------------------------------------------------------------
// wrapper for the block basis count
//-----------------------------------------------------------------------------
void count_by_process_block_gpu(G4int nTracks, GXTrack *tracks, 
				G4int *nbrem, G4int *nioni, 
				int blocksPerGrid, int threadsPerBlock,
                                cudaStream_t stream)
{
  count_by_process_block_kernel<<< blocksPerGrid, threadsPerBlock, 
     2*threadsPerBlock*sizeof(G4int), stream >>> (nTracks, tracks, nbrem, nioni);
}

//-----------------------------------------------------------------------------
// sort tracks by the physics process
//-----------------------------------------------------------------------------
__global__
void sort_by_process_kernel(G4int nTracks, GXTrack *itracks, GXTrack *stracks, 
			    G4int *nbrem, G4int *stackSize_brem,
			    G4int *nioni, G4int *stackSize_ioni)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < nTracks) {
    //offset is a global counter for the last array position of secondaries 
    if(itracks[tid].proc == 0 ) {
      stracks[atomicAdd(stackSize_brem,1)] = itracks[tid];
    }    
    else if(itracks[tid].proc == 1 ) {
      stracks[*nbrem + atomicAdd(stackSize_ioni,1)] = itracks[tid];
    }    
    else {
      printf("@@@ Missing Counter for in GPU %d @@@\n",itracks[tid].proc);
    }
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// sort tracks by the physics process
//-----------------------------------------------------------------------------
__global__
void sort_by_process_kernel_inplace_part1(G4int nTracks, GXTrack *itracks, G4int *index, 
			    G4int *nbrem, G4int *stackSize_brem,
			    G4int *nioni, G4int *stackSize_ioni)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  int loc_nbrem = *nbrem;
  //int loc_nioni = *nioni;

  // Record misplaced tracks
  while(tid < nTracks) {
    if(itracks[tid].proc == 0 && tid > loc_nbrem ) {
      index[atomicAdd(stackSize_brem,1)] = tid;
    }    
    else if(itracks[tid].proc == 1 && tid < loc_nbrem) {
      index[loc_nbrem + atomicAdd(stackSize_ioni,1)] = tid;
    }    
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void sort_by_process_kernel_inplace_part2(G4int nTracks, GXTrack *itracks, G4int *index, 
			    G4int *nbrem, G4int *stackSize_brem,
			    G4int *nioni, G4int *stackSize_ioni)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  int loc_nbrem = *nbrem;
  //int loc_nioni = *nioni;

  // We must have *stackSize_brem == *stackSize_ioni (i.e. equal number to be moved)
  // and *stackSizeBrem < nTracks
  int nProcess = *stackSize_brem;

  // Move misplaced tracks
  while(tid < nProcess) {
    int mis_brem = index[loc_nbrem + tid];
    GXTrack temp = itracks[ mis_brem ];   // Save brem track
    itracks[ mis_brem ] = itracks[ tid ] ; // Move ioni track in place
    itracks[ tid ] = temp;                 // Put brem track in place

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// wrapper for sort
//-----------------------------------------------------------------------------

void sort_by_process_gpu(G4int nTracks, GXTrack *itracks, GXTrack *stracks, 
			 G4int *nbrem, G4int *stackSize_brem, 
			 G4int *nioni, G4int *stackSize_ioni, 
			 int blocksPerGrid, int threadsPerBlock,
                         cudaStream_t stream)
{
   sort_by_process_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
    (nTracks, itracks, stracks, nbrem, stackSize_brem, nioni, stackSize_ioni);
}

//-----------------------------------------------------------------------------
// sort on CPU
//-----------------------------------------------------------------------------

void sort_by_process_cpu(G4int nTracks, GXTrack *itracks, GXTrack *stracks,
			 G4int *nbrem, G4int *stackSize_brem,
			 G4int *nioni, G4int *stackSize_ioni)
{
  for(int tid = 0 ; tid < nTracks ; tid++) {
    //offset is a global counter for the last array position of secondaries 
    if(itracks[tid].proc == 0 ) {
      stracks[*stackSize_brem] = itracks[tid];
      ++(*stackSize_brem);
    }
    else if(itracks[tid].proc == 1 ) {
      stracks[*nbrem + *stackSize_ioni] = itracks[tid];
      ++(*stackSize_ioni);
    }
    else {
      printf("@@@ Missing Counter in CPU for %d @@@\n",itracks[tid].proc);
    }
  }
}

