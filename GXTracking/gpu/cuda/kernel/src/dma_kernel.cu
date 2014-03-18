#include "stdio.h"
#include <stdlib.h>

#include <curand_kernel.h>

#include "dma_kernel.h"

//-----------------------------------------------------------------------------
// utility functions
//-----------------------------------------------------------------------------

__device__ int pmutex = 0;

__device__ void lock(int *pmutex)
{
  while(atomicCAS(pmutex, 0, 1) != 0);
}

__device__ void unlock(int *pmutex)
{
  atomicExch(pmutex, 0);
}

__device__
void fill_track_gpu(curandState *localState, GXTrack *aTrack) 
{
  aTrack->x  = 300*(2.0*curand_uniform(localState)-1.0);
  aTrack->y  = 300*(2.0*curand_uniform(localState)-1.0);
  aTrack->z  = 300*(2.0*curand_uniform(localState)-1.0);
  aTrack->s  = 100*(2.0*curand_uniform(localState)-1.0);
  aTrack->px = 100*(2.0*curand_uniform(localState)-1.0);
  aTrack->py = 100*(2.0*curand_uniform(localState)-1.0);
  aTrack->pz = 100*(2.0*curand_uniform(localState)-1.0);
  aTrack->q  = -1;
}

void fill_track_cpu(GXTrack *aTrack) 
{
  aTrack->x  = 300*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->y  = 300*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->z  = 300*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->s  = 100*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->px = 100*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->py = 100*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->pz = 100*(2.0*((G4double)rand()/RAND_MAX)-1.0);
  aTrack->q  = -1;
}

//-----------------------------------------------------------------------------
// test pre-allocated memory on the device with a fixed size and
// write each secondary particle directly to it from each thread
//-----------------------------------------------------------------------------

__global__
void fixed_wof_kernel(curandState *devStates, G4int nTracks, 
		      GXTrack *secTracks, G4int *stackSize, G4int *offset)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  curandState localState = devStates[tid];

  G4int nsec = 0;
  GXTrack aTrack;

  while(tid < nTracks) {

    nsec = (maxSecondaryPerStep+1)*curand_uniform(&localState);
     
    //generate secondaries and directly store them to the stack on the heap
    for(int isec = 0 ; isec < nsec ; isec++) {

      fill_track_gpu(&localState, &aTrack); 

      //offset is a global counter for the last array position of secondaries 
      *offset = atomicAdd(stackSize,1);
      secTracks[*offset] = aTrack;
    }
    
    tid += blockDim.x * gridDim.x;
  }
  
}

//-----------------------------------------------------------------------------
// wrapper for testing global memory allocation with a fixed size
//-----------------------------------------------------------------------------

void fixed_wof_gpu(curandState *devStates, G4int nTracks, 
		   GXTrack *secTracks, G4int *stackSize, G4int *offset,
		   int blocksPerGrid, int threadsPerBlock)
{
  fixed_wof_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates, nTracks, secTracks, stackSize, offset);
}

//-----------------------------------------------------------------------------
// test fixed_wof on CPU
//-----------------------------------------------------------------------------

void fixed_wof_cpu(G4int nTracks, GXTrack *secTracks, G4int *stackSize)
{
  G4int nsec = 0;
  GXTrack aTrack;

  //generate secondaries
  for(int tid = 0 ; tid < nTracks ; tid++) {

    nsec = (maxSecondaryPerStep+1)*((G4double)rand()/RAND_MAX);
     
    //generate secondaries and directly store them to the stack on the heap
    for(int isec = 0 ; isec < nsec ; isec++) {

      fill_track_cpu(&aTrack); 

      //offset is a global counter for the last array position of secondaries 
      secTracks[*stackSize] = aTrack;
      ++(*stackSize);
    }
  }
}


//-----------------------------------------------------------------------------
// test pre-allocated memory on the device with a fixed size, and store 
// secondary particles locally on the thread and copy them to the global memory
//-----------------------------------------------------------------------------

__global__
void fixed_sac_kernel(curandState *devStates, G4int nTracks, 
		      GXTrack *secTracks, G4int *stackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  curandState localState = devStates[tid];

  G4int nsec = 0;
  GXTrack aTrack;

  while(tid < nTracks) {

    nsec = (maxSecondaryPerStep+1)*curand_uniform(&localState);
     
    // do this even for the case that there is no secondary to avoid divergence
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));
     
    //generate secondaries and store them locally on the thread
    for(int isec = 0 ; isec < nsec ; isec++) {

      fill_track_gpu(&localState, &aTrack); 

      secondaries[isec] = aTrack;

    }
    
    //offset returns old value of stackSize after stackSize += nsec
    int offset = atomicAdd(stackSize,nsec);

    //copy all secondaries produced by this thread to the global memory
    memcpy(secTracks+offset,secondaries,nsec*sizeof(GXTrack));

    //free the thread local storage
    free(secondaries);
     
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// wrapper for testing a fixed global memory allocation 
//-----------------------------------------------------------------------------

void fixed_sac_gpu(curandState *devStates, G4int nTracks, 
		   GXTrack *secTracks, G4int *stackSize,
		   int blocksPerGrid, int threadsPerBlock)
{
  fixed_sac_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates, nTracks, secTracks, stackSize);
}

//-----------------------------------------------------------------------------
// test fixed_sac on CPU
//-----------------------------------------------------------------------------

void fixed_sac_cpu(G4int nTracks, GXTrack *secTracks, G4int *stackSize)
{
  G4int nsec = 0;
  GXTrack aTrack;

  for(int tid = 0 ; tid < nTracks ; tid++) {

    nsec = (maxSecondaryPerStep+1)*((G4double)rand()/RAND_MAX);
     
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));
     
    //generate secondaries and store them locally on the thread
    for(int isec = 0 ; isec < nsec ; isec++) {

      //generate and fill each secondary
      fill_track_cpu(&aTrack); 

      secondaries[isec] = aTrack;
    }
    
    //copy all secondaries produced by this thread to the global memory
    memcpy(secTracks+(*stackSize),secondaries,nsec*sizeof(GXTrack));
    ++(*stackSize);

    //free the thread local storage
    free(secondaries);
  }
}

//-----------------------------------------------------------------------------
// test dynamic memory allocation 
//-----------------------------------------------------------------------------

__global__
void dma_alloc_kernel(curandState *devStates, G4int nTracks, 
		      GXSecContainer *container, G4int *stackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  curandState localState = devStates[tid];

  G4int nsec = 0;
  GXTrack aTrack;

  while(tid < nTracks) {

    nsec = (maxSecondaryPerStep+1)*curand_uniform(&localState);
     
    // do this for the case that there is no secondaries to avoid divergence
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));
     
    //generate secondaries
    for(int isec = 0 ; isec < nsec ; isec++) {

      //generate and fill each secondary
      fill_track_gpu(&localState, &aTrack); 

      secondaries[isec] = aTrack;
    }
     
    //store the size and the pointer for secondaries
    container[tid].size = nsec;
    container[tid].addr = secondaries;
    
    //add the number of secondaries produced by all threads 
    //and send the total number of secondaries to the host
    atomicAdd(stackSize,nsec);

    tid += blockDim.x * gridDim.x;

  }
}

//-----------------------------------------------------------------------------
// wrapper for the test dynamic memory allocation 
//-----------------------------------------------------------------------------

void dma_alloc_gpu(curandState *devStates, G4int nTracks, 
		   GXSecContainer *container, G4int *stackSize,
		   int blocksPerGrid, int threadsPerBlock)
{
  dma_alloc_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates, nTracks, container, stackSize);
}

//-----------------------------------------------------------------------------
// reallocate dynamic memory allocation for a single stack
//-----------------------------------------------------------------------------

__global__
void dma_realloc_kernel(int nTracks, GXSecContainer *container, 
			GXTrack *secTracks, G4int *stackSize, G4int *offset) {

  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nTracks) {

    if( container[tid].size > 0 ) {
      *offset = atomicAdd(stackSize,container[tid].size);
      memcpy(secTracks+(*offset),container[tid].addr,
	     container[tid].size*sizeof(GXTrack));
    }

    tid += blockDim.x * gridDim.x;
  }

}
//-----------------------------------------------------------------------------
// wrapper for dma_realloc_kernel
//-----------------------------------------------------------------------------

void dma_realloc_gpu(G4int nTracks, GXSecContainer *container, 
		     GXTrack *secTracks, G4int *stackSize, G4int *offset,
		     int blocksPerGrid, int threadsPerBlock)
{
  dma_realloc_kernel<<< blocksPerGrid, threadsPerBlock>>>
    (nTracks,container,secTracks, stackSize, offset);
}

//-----------------------------------------------------------------------------
// free the dynamically allocated memory
//-----------------------------------------------------------------------------

__global__
void dma_free_kernel(int nTracks, GXSecContainer *container) {

  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
  while (tid < nTracks) {
    free(container[tid].addr);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// wrapper for dma_free_kernel
//-----------------------------------------------------------------------------

void dma_free_gpu(int nTracks, GXSecContainer *container, 
		  int blocksPerGrid, int threadsPerBlock)
{
  dma_free_kernel<<< blocksPerGrid, threadsPerBlock >>>(nTracks,container);
}

//-----------------------------------------------------------------------------
// test dynamic memory allocation on CPU
//-----------------------------------------------------------------------------

void dma_alloc_cpu(G4int nTracks)
{
  G4int stackSize = 0;
  GXTrack aTrack;

  GXSecContainer *container 
    = (GXSecContainer*) malloc(nTracks*sizeof(GXSecContainer));

  for(int tid = 0 ; tid < nTracks ; tid++) {

    G4int nsec = 0;
    nsec = (maxSecondaryPerStep+1)*((G4double)rand()/RAND_MAX);
     
    // do this for the case that there is no secondaries to avoid divergence
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));
     
    //generate secondaries
    for(int isec = 0 ; isec < nsec ; isec++) {

      fill_track_cpu(&aTrack); 
      secondaries[isec] = aTrack;
    }
     
    container[tid].size = nsec;
    container[tid].addr = secondaries;

    stackSize += nsec;
  }

  GXTrack *secTracks = (GXTrack*)malloc(stackSize*sizeof(GXTrack));

  G4int offset = 0;

  for(int it = 0 ; it < nTracks ; it++) {
    memcpy(secTracks+offset,container[it].addr,
	   container[it].size*sizeof(GXTrack));
    free(container[it].addr);
    offset += container[it].size; 
  }

  free(secTracks);
  
}

//-----------------------------------------------------------------------------
// print 
//-----------------------------------------------------------------------------

__global__
void dma_print_kernel(G4int stackSize, GXTrack *secTracks) {

  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
  while (tid < stackSize) {
    printf("GPU Stack (x,y,z,s)= (%f,%f,%f,%f)\n",secTracks[tid].x,
	   secTracks[tid].y,secTracks[tid].z,secTracks[tid].s);
    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
// wrapper for dma_print_kernel
//-----------------------------------------------------------------------------

void dma_print_gpu(G4int stackSize, GXTrack *secTracks, 
		   int blocksPerGrid, int threadsPerBlock)
{
  dma_print_kernel<<< blocksPerGrid, threadsPerBlock >>>(stackSize,secTracks);
}

//-----------------------------------------------------------------------------
// extended tests: 4. dynamic memory alloc/realloc per thread basis
//                 5. dynamic memory alloc/realloc per block basis
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// 4. test dynamic memory allocation per thread basis
//-----------------------------------------------------------------------------

__global__
void dma_alloc_perthread_kernel(curandState *devStates, G4int nTracks, 
				GXSecContainer *container, G4int *stackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  curandState localState = devStates[tid];

  // estimated the maximum number of secondaries produced in this thread = 
  // n tracks processed in this thread * maximum number of secondaries per track
  // = [nTracks/(blockDim.x*gridDim.x) + 1] * maxSecondaryPerStep

  G4int maxSec = (nTracks%(blockDim.x*gridDim.x) !=0 ) ?
    (nTracks/(blockDim.x*gridDim.x) + 1)*maxSecondaryPerStep : 
    (nTracks/(blockDim.x*gridDim.x))*maxSecondaryPerStep ;

  //temporary buffer on the global (heap) memory
  GXTrack* secBuffer = (GXTrack*)malloc(maxSec*sizeof(GXTrack));

  G4int ibuffer = 0;
  G4int nsec = 0;
  G4int tsec = 0;

  GXTrack aTrack;
  G4int threadIndex = tid;

  while(tid < nTracks) {

    //generate secondaries and store in the temporary buffer
    nsec = (maxSecondaryPerStep+1)*curand_uniform(&localState);
     
    for(int isec = 0 ; isec < nsec ; isec++) {
      fill_track_gpu(&localState, &aTrack); 
      secBuffer[ibuffer] = aTrack;
      ibuffer++;
    }
    tsec += nsec;
    tid += blockDim.x * gridDim.x;
  }

  atomicAdd(stackSize,tsec);

  //allocate and copy the size of secondaries at the pre-allocated memory 
  container[threadIndex].size = tsec;
  container[threadIndex].addr = (GXTrack*)malloc(tsec*sizeof(GXTrack));
  memcpy(container[threadIndex].addr,secBuffer,tsec*sizeof(GXTrack));

  //free the local memory
  free(secBuffer);
}

//-----------------------------------------------------------------------------
// wrapper for testing dynamic memory allocation per thread basis
//-----------------------------------------------------------------------------

void dma_alloc_perthread_gpu(curandState *devStates, G4int nTracks, 
			     GXSecContainer *container, G4int *stackSize,
			     int blocksPerGrid, int threadsPerBlock)
{
  dma_alloc_perthread_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates, nTracks, container, stackSize);
}

//-----------------------------------------------------------------------------
// reallocate the stored secondary arrays into a single stack per thread basis
//-----------------------------------------------------------------------------

__global__
void dma_realloc_perthread_kernel(int nTracks, GXSecContainer *container, 
				  GXTrack *secTracks, G4int *stackSize, 
				  G4int *offset) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if( container[tid].size > 0 ) {
    *offset = atomicAdd(stackSize,container[tid].size);
    memcpy(secTracks+(*offset),container[tid].addr,
	   container[tid].size*sizeof(GXTrack));
  }
}

//-----------------------------------------------------------------------------
// wrapper for dma_realloc_perthread_kernel
//-----------------------------------------------------------------------------

void dma_realloc_perthread_gpu(G4int nTracks, GXSecContainer *container, 
			       GXTrack *secTracks, G4int *stackSize, 
			       G4int *offset,
			       int blocksPerGrid, int threadsPerBlock)
{
  dma_realloc_perthread_kernel<<< blocksPerGrid, threadsPerBlock>>>
    (nTracks,container,secTracks, stackSize, offset);
}

//-----------------------------------------------------------------------------
// free the dynamically allocated memory
//-----------------------------------------------------------------------------

__global__
void dma_free_perthread_kernel(int nTracks, GXSecContainer *container) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  free(container[tid].addr);
}

//-----------------------------------------------------------------------------
// wrapper for dma_free_perthread_kernel
//-----------------------------------------------------------------------------

void dma_free_perthread_gpu(int nTracks, GXSecContainer *container, 
			    int blocksPerGrid, int threadsPerBlock)
{
  dma_free_perthread_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (nTracks,container);
}

//-----------------------------------------------------------------------------
// 5. test dynamic memory allocation per block basis
//-----------------------------------------------------------------------------

__global__
void dma_alloc_perblock_kernel(curandState *devStates, G4int nTracks, 
			       GXSecContainer *container, G4int *stackSize)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  curandState localState = devStates[tid];

  G4int blockIndex;

  // estimated the maximum number of secondaries produced in this block = 
  // n tracks processed in this block * maximum number of secondaries per track
  // = [nTracks/(gridDim.x) + 1] * maxSecondaryPerStep

  G4int maxSec = (nTracks%gridDim.x !=0 ) ?
    (nTracks/gridDim.x + 1)*maxSecondaryPerStep : 
    (nTracks/gridDim.x)*maxSecondaryPerStep ;

  //temporary buffer on the shared memory 
  //maximum amount of shared memory per multiprocess (Fermi) = 48KB
  //maxSec*sizeof(GXTrack) < 48KB/??

  __shared__ GXTrack* secBuffer;
  __shared__ int* tsec;

  if(threadIdx.x ==0 ) {
    blockIndex = blockIdx.x;

    secBuffer = (GXTrack*)malloc(maxSec*sizeof(GXTrack));
    tsec = (int*)malloc(sizeof(int));
    *tsec = 0;
  }
  __syncthreads();

  G4int nsec = 0;
  GXTrack aTrack;

  while(tid < nTracks) {

    //generate secondaries and store in the temporary buffer
    nsec = (maxSecondaryPerStep+1)*curand_uniform(&localState);
     
    for(int isec = 0 ; isec < nsec ; isec++) {

      fill_track_gpu(&localState, &aTrack); 

      //atomicAdd(tsec,1) will return the last position of the shared array

      secBuffer[atomicAdd(tsec,1)] = aTrack;
    }
    tid += blockDim.x * gridDim.x;
  }

  //ensure that all threads complete before copying
  __syncthreads();

  //store the size/address of secondaries and copy them to the global memory
  if(threadIdx.x ==0 ) {
    atomicAdd(stackSize,*tsec);

    container[blockIndex].size = *tsec;
    container[blockIndex].addr = (GXTrack*)malloc((*tsec)*sizeof(GXTrack));
    memcpy(container[blockIndex].addr,secBuffer,(*tsec)*sizeof(GXTrack));
 
    //free the shared memory
    free(secBuffer);
    free(tsec);
  }
}

//-----------------------------------------------------------------------------
// wrapper for testing dynamic memory allocation per block basis
//-----------------------------------------------------------------------------

void dma_alloc_perblock_gpu(curandState *devStates, G4int nTracks, 
			    GXSecContainer *container, G4int *stackSize,
			    int blocksPerGrid, int threadsPerBlock)
{
  dma_alloc_perblock_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (devStates, nTracks, container, stackSize);
}

//-----------------------------------------------------------------------------
// reallocate the stored secondary arrays into a single stack per block basis
//-----------------------------------------------------------------------------

__global__
void dma_realloc_perblock_kernel(int nTracks, GXSecContainer *container, 
				 GXTrack *secTracks, G4int *stackSize, 
				 G4int *offset) 
{
  unsigned int tid = blockIdx.x;
  if( container[tid].size > 0 ) {
    *offset = atomicAdd(stackSize,container[tid].size);
    memcpy(secTracks+(*offset),container[tid].addr,
	   container[tid].size*sizeof(GXTrack));
  }
}

//-----------------------------------------------------------------------------
// wrapper for dma_realloc_perblock_kernel
//-----------------------------------------------------------------------------

void dma_realloc_perblock_gpu(G4int nTracks, GXSecContainer *container, 
			      GXTrack *secTracks, G4int *stackSize, 
			      G4int *offset,
			      int blocksPerGrid, int threadsPerBlock)
{
  dma_realloc_perblock_kernel<<< blocksPerGrid, threadsPerBlock>>>
    (nTracks,container,secTracks, stackSize, offset);
}

//-----------------------------------------------------------------------------
// free the dynamically allocated memory
//-----------------------------------------------------------------------------

__global__
void dma_free_perblock_kernel(int nTracks, GXSecContainer *container) 
{
  unsigned int tid = blockIdx.x;
  free(container[tid].addr);
}

//-----------------------------------------------------------------------------
// wrapper for dma_free_perblock_kernel
//-----------------------------------------------------------------------------

void dma_free_perblock_gpu(int nTracks, GXSecContainer *container, 
			   int blocksPerGrid, int threadsPerBlock)
{
  dma_free_perblock_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (nTracks,container);
}

//-----------------------------------------------------------------------------
// include files
//-----------------------------------------------------------------------------

#include "random_kernel.cu"
