//----------------------------------------------------------------------------
// a test code for handling global memory for secondary tracks produced on GPU
// March 7, 2013
//----------------------------------------------------------------------------
#include <iostream>
#include <ctime>
#include "stdio.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "dma_kernel.h"
#include "random_kernel.h"

#include "GXTrack.h"

int main(int argc, char* argv[])
{
  //argument
  int kernelType = 1; //1 fixed write-on-the-fly 
                      //2 fixed store-and-copy   
                      //3 dynamic memory allocation per track bisis
                      //4 dynamic memory allocation per thread bisis
                      //5 dynamic memory allocation per block bisis
  int testType = 1;   //1 default 2 test

  if(argc >= 2) kernelType = atoi(argv[1]);
  if(argc >= 3) testType = atoi(argv[2]);

  if(kernelType < 1 || kernelType > 5 || testType < 1 || testType > 2) {
    std::cout << "Usage: dma [1=wof|2=sac|3=pertrack|4=perthread|5=perblock] " 
	      << "[1=notest|2=test]... " 
	      << std::endl;
    return 0;
  }

  //cudaEvent for performance measurement
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  //maximum size of dynamic memory allocation on the device
  //set 128 megabytes on the heap (global mememory) 
  cudaThreadSetLimit(cudaLimitMallocHeapSize, 128*1024*1024);

  //thread organization: default threads and blocks
  int theNBlocks  =  32;
  int theNThreads = 128;
  if(testType ==2) {
    theNBlocks  =  2;
    theNThreads =  4;
  }

  //event-loop
  int nevent = 10;

  for (int i = 0 ; i < nevent ; i++) {

    //number of input tracks
    int nTracks = 100000;
    if(testType==2) nTracks = 8;
    
    //atomic counter for the number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    //atomic counter for the last array position of scondaries
    G4int *offset_d;
    cudaMalloc((void**)&offset_d, sizeof(G4int));
    cudaMemset(offset_d,0,sizeof(G4int));

    //construct default random states
    curandState* devStates = 0;
    cudaMalloc(&devStates,theNBlocks*theNThreads*sizeof(curandState));
    curand_setup_gpu(devStates, time(NULL)+i, theNBlocks, theNThreads);
    
    //reset clock
    float elapsedTimeAlloc = 0.0;
    float elapsedTimeRealloc = 0.0;
    float elapsedTimeCPU = 0.0;

    G4int stackSize = 0;
    GXTrack *secTracks_d;

    //start time for GPU kernel(s)
    cudaEventRecord (start,0);

    if(kernelType ==1) {

      //allocate fixed memory size on the device 
      cudaMalloc((void**)&secTracks_d, 
		 maxSecondaryPerStep*nTracks*sizeof(GXTrack));

      //kernel
      fixed_wof_gpu(devStates,nTracks,secTracks_d,stackSize_d,
		    offset_d,theNBlocks,theNThreads);

      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);
      printf("Number of Secondaries Created on GPU = %d\n",stackSize);

      if(testType==2) {
	dma_print_gpu(stackSize,secTracks_d,theNBlocks,theNThreads);
      }
    }
    else if(kernelType ==2) {
      //allocate fixed memory size on the device 
      cudaMalloc((void**)&secTracks_d, 
		 maxSecondaryPerStep*nTracks*sizeof(GXTrack));

      //kernel
      fixed_sac_gpu(devStates,nTracks,secTracks_d,stackSize_d,
		    theNBlocks,theNThreads);

      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);
      printf("Number of Secondaries Created on GPU = %d\n",stackSize);

      if(testType==2) {
	dma_print_gpu(stackSize,secTracks_d,theNBlocks,theNThreads);
      }

    }
    else if(kernelType ==3) {

      //memory allocation for a containter holding secondaries information
      GXSecContainer *secContainer_d;
      cudaMalloc((void**)&secContainer_d, nTracks*sizeof(GXSecContainer));

      //create and fill secondaries on dymamically allocated memory
      //<--------------------------------------------------------
      dma_alloc_gpu(devStates,nTracks,secContainer_d,stackSize_d,
		    theNBlocks,theNThreads);
      //<--------------------------------------------------------
      cudaThreadSynchronize();
      
      //get the number of secondaries created on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      //stop time for kernel
      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
      
      printf("Number of Secondaries Created on GPU = %d\n",stackSize);
      
      //reallocate memory 
      cudaEventRecord (start,0);

      cudaMalloc((void**)&secTracks_d, stackSize*sizeof(GXTrack));

      //reset the size of stackSize
      stackSize = 0;
      cudaMemset(stackSize_d,0,sizeof(G4int));

      //reallocate secondaries into one array
      //<------------------------------------------------------------
      dma_realloc_gpu(nTracks,secContainer_d,secTracks_d,stackSize_d,
      		      offset_d,theNBlocks,theNThreads);
      //<------------------------------------------------------------
      cudaThreadSynchronize();

      //get the number of reallocated secondaries on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      if(testType==2) {
	dma_print_gpu(stackSize,secTracks_d,theNBlocks,theNThreads);
      }

      dma_free_gpu(nTracks,secContainer_d,theNBlocks,theNThreads);

      cudaFree(secContainer_d);
    }
    else if(kernelType ==4) {

      //memory allocation for a containter holding secondaries per thread 
      GXSecContainer *secContainer_d;
      cudaMalloc((void**)&secContainer_d, 
		 theNBlocks*theNThreads*sizeof(GXSecContainer));

      //create and fill secondaries on dymamically allocated memory
      //<------------------------------------------------------------------
      dma_alloc_perthread_gpu(devStates,nTracks,secContainer_d,stackSize_d,
			      theNBlocks,theNThreads);
      //<------------------------------------------------------------------
      cudaThreadSynchronize();
      
      //get the number of secondaries created on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      //stop time for kernel
      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
      
      printf("Number of Secondaries Created on GPU = %d\n",stackSize);
      
      //reallocate memory 
      cudaEventRecord (start,0);

      cudaMalloc((void**)&secTracks_d, stackSize*sizeof(GXTrack));

      //reset the size of stackSize
      stackSize = 0;
      cudaMemset(stackSize_d,0,sizeof(G4int));

      //reallocate secondaries into one array
      //<----------------------------------------------------------------------
      dma_realloc_perthread_gpu(nTracks,secContainer_d,secTracks_d,stackSize_d,
				offset_d,theNBlocks,theNThreads);
      //<----------------------------------------------------------------------
      cudaThreadSynchronize();

      //get the number of reallocated secondaries on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      if(testType==2) {
	dma_print_gpu(stackSize,secTracks_d,theNBlocks,theNThreads);
      }

      dma_free_perthread_gpu(nTracks,secContainer_d,theNBlocks,theNThreads);
      cudaFree(secContainer_d);
    }
    else if(kernelType ==5) {

      //memory allocation for a containter holding secondaries per block 
      GXSecContainer *secContainer_d;
      cudaMalloc((void**)&secContainer_d, 
		 theNBlocks*sizeof(GXSecContainer));

      //create and fill secondaries on dymamically allocated memory
      //<-----------------------------------------------------------------
      dma_alloc_perblock_gpu(devStates,nTracks,secContainer_d,stackSize_d,
			     theNBlocks,theNThreads);
      //<-----------------------------------------------------------------
      cudaThreadSynchronize();
      
      //get the number of secondaries created on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      //stop time for kernel
      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
      
      printf("Number of Secondaries Created on GPU = %d\n",stackSize);
      
      //reallocate memory 
      cudaEventRecord (start,0);

      cudaMalloc((void**)&secTracks_d, stackSize*sizeof(GXTrack));

      //reset the size of stackSize
      stackSize = 0;
      cudaMemset(stackSize_d,0,sizeof(G4int));

      //reallocate secondaries into one array
      //<---------------------------------------------------------------------
      dma_realloc_perblock_gpu(nTracks,secContainer_d,secTracks_d,stackSize_d,
			       offset_d,theNBlocks,1);
      //<---------------------------------------------------------------------
      cudaThreadSynchronize();

      //get the number of reallocated secondaries on GPU
      cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),cudaMemcpyDeviceToHost);

      if(testType==2) {
	dma_print_gpu(stackSize,secTracks_d,theNBlocks,theNThreads);
      }

      dma_free_perblock_gpu(nTracks,secContainer_d,theNBlocks,1);
      cudaFree(secContainer_d);
    }

    //stop time for reallocation
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeRealloc,start,stop);
    
    //Repeat on CPU 

    stackSize = 0;
    GXTrack *secTracks_h 
      = (GXTrack*) malloc(maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    cudaEventRecord (start,0);

    if(kernelType ==1) {
      fixed_wof_cpu(nTracks, secTracks_h, &stackSize);
    }
    else if(kernelType ==2) {
      fixed_sac_cpu(nTracks, secTracks_h, &stackSize);
    }
    else if(kernelType ==3 || kernelType ==4 || kernelType ==5) {
      dma_alloc_cpu(nTracks);
    }

    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeCPU,start,stop);

    //Print Performance

    if(kernelType ==1 || kernelType ==2) {
      printf("Elapsed Time: GPU CPU [msec]: %6.3f %6.3f\n",
	     elapsedTimeRealloc,elapsedTimeCPU);
    }
    else if(kernelType ==3 || kernelType ==4 || kernelType ==5) {
      printf("Elapsed Time: Alloc Realloc GPU CPU [msec]: %6.3f %6.3f %6.3f %6.3f\n",
	     elapsedTimeAlloc,elapsedTimeRealloc,
	     elapsedTimeAlloc+elapsedTimeRealloc,elapsedTimeCPU);
    }

    //clean up    
    cudaFree(devStates);
    cudaFree(stackSize_d);
    cudaFree(offset_d);
    cudaFree(secTracks_d);
    free(secTracks_h);

  }
  //end of event-loop
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

}
