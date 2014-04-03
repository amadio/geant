#ifndef DMA_KERNEL_H
#define DMA_KERNEL_H

#include <curand_kernel.h>

#include "GPTypeDef.h"
#include "GXTrack.h"

//the maximum number of secondaries per step
const G4int maxSecondaryPerStep = 2;

typedef struct { 
  size_t size;
  GXTrack* addr; 
} GXSecContainer;


extern void fixed_wof_gpu(curandState *devStates,
			  G4int nTracks,
			  GXTrack *secTrack,
			  G4int *stackSize,
			  G4int *offset,
			  int nBlocks,
			  int nThreads); 

extern void fixed_wof_cpu(G4int nTracks,
			  GXTrack *secTrack,
			  G4int *stackSize); 

extern void fixed_sac_gpu(curandState *devStates,
			  G4int nTracks,
			  GXTrack *secTrack,
			  G4int *stackSize,
			  int nBlocks,
			  int nThreads); 

extern void fixed_sac_cpu(G4int nTracks,
			  GXTrack *secTrack,
			  G4int *stackSize);

extern void dma_alloc_gpu(curandState *devStates,
			  G4int nTracks,
			  GXSecContainer *container,
			  G4int *stackSize,
			  int nBlocks,
			  int nThreads); 

extern void dma_realloc_gpu(G4int nTracks,
			    GXSecContainer *container,
			    GXTrack *secTrack,
			    G4int *stackSize,
			    G4int *offset,
			    int nBlocks,
			    int nThreads); 

extern void dma_free_gpu(G4int nTracks,
			 GXSecContainer *container,
			 int nBlocks,
			 int nThreads); 

extern void dma_alloc_cpu(int nTracks); 

extern void dma_print_gpu(G4int stackSize,
			  GXTrack *secTracks,
			  int nBlocks,
			  int nThreads); 

//extended tests: per thread basis
extern void dma_alloc_perthread_gpu(curandState *devStates,
				    G4int nTracks,
				    GXSecContainer *container,
				    G4int *stackSize,
				    int nBlocks,
				    int nThreads); 

extern void dma_realloc_perthread_gpu(G4int nTracks,
				      GXSecContainer *container,
				      GXTrack *secTrack,
				      G4int *stackSize,
				      G4int *offset,
				      int nBlocks,
				      int nThreads); 

extern void dma_free_perthread_gpu(G4int nTracks,
				   GXSecContainer *container,
				   int nBlocks,
				   int nThreads); 

//extended tests: per block basis
extern void dma_alloc_perblock_gpu(curandState *devStates,
				   G4int nTracks,
				   GXSecContainer *container,
				   G4int *stackSize,
				   int nBlocks,
				   int nThreads); 

extern void dma_realloc_perblock_gpu(G4int nTracks,
				     GXSecContainer *container,
				     GXTrack *secTrack,
				     G4int *stackSize,
				     G4int *offset,
				     int nBlocks,
				     int nThreads); 

extern void dma_free_perblock_gpu(G4int nTracks,
				  GXSecContainer *container,
				  int nBlocks,
				  int nThreads); 

#endif
