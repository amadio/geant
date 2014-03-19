#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "random_kernel.h"
#include "rngTest_kernel.h"

#include <curand_kernel.h>
#include <curand_mtgp32_host.h>
#include <curand_mtgp32.h>
#include <curand_mtgp32dc_p_11213.h>

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument
  int rngType = 1; //1 XORWOW, 2 MRG32k3a, 3MTgt32

  if(argc >= 2) rngType = atoi(argv[1]);
  if(rngType>3) {
    std::cout << "Wrong rngType (>3): ... Bailing out ... "<< std::endl;
    return 0;
  }

 //cuda interface
  bool  cudaFirst;
  char* cudaEnabled;
  char* cudaNBlocks;
  char* cudaNThreads;

  int theNBlocks  =  32;
  int theNThreads = 128;

  cudaDeviceReset();

  //event-loop
  int nevent = 10;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 100000;

    float elapsedTimeInit = 0.0;
    float elapsedTimeGen  = 0.0;
    float elapsedTimeGPU  = 0.0;

    //prepare curand states
    curandState* devStates = 0;
    curandStateMRG32k3a* devStatesMRG32k3a = 0;

    curandStateMtgp32* devStatesMtgp32 = 0;
    mtgp32_kernel_params *kp;

    //return values
    int *result_h;
    int *result_c;
    int *result_d;

    result_h = (int*) calloc(theNBlocks*theNThreads, sizeof(int));
    result_c = (int*) calloc(theNBlocks*theNThreads, sizeof(int));
    cudaMalloc((void**)&result_d,theNBlocks*theNThreads*sizeof(int));

    cudaMemset(result_d,0,theNBlocks*theNThreads*sizeof(int));

    //start time for init
    cudaEvent_t start, stop;
    cudaEventCreate (&start);
    cudaEventCreate (&stop);
    cudaEventRecord (start,0);

    //initialize random engines on the device
    if(rngType == 1) {
      cudaMalloc(&devStates,theNBlocks*theNThreads*sizeof(curandState));
      curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);
    }
    else if(rngType == 2) {
      cudaMalloc(&devStatesMRG32k3a,
		 theNBlocks*theNThreads*sizeof(curandStateMRG32k3a));
      curand_setup_gpu(devStatesMRG32k3a, time(NULL), theNBlocks, theNThreads);
    }    
    else if(rngType == 3) {
      //MTGP32 allows thread-safe generation and state update for up to 256
      //concurrent threads (within a single block) for each of the 200 sequences
      //theNBlocks <= 200 && theNThreads <= 256
      cudaMalloc((void**)&devStatesMtgp32,theNBlocks*sizeof(curandStateMtgp32));

      cudaMalloc((void**)&kp,sizeof(mtgp32_kernel_params));
      curandMakeMTGP32Constants(mtgp32dc_params_fast_11213,kp);
      curandMakeMTGP32KernelState(devStatesMtgp32,mtgp32dc_params_fast_11213,kp,
				  theNBlocks,time(NULL));
    }

    //stop time for init
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop); 
    cudaEventElapsedTime (&elapsedTimeInit,start,stop);

    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(rngType == 1) {
      XORWOW_gpu(devStates,result_d,nTracks,theNBlocks,theNThreads);
    }
    else if(rngType == 2) {
      MRG32k3a_gpu(devStatesMRG32k3a,result_d,nTracks,theNBlocks,theNThreads);
    }
    else if(rngType == 3) {
      Mtgp32_gpu(devStatesMtgp32,result_d,nTracks,theNBlocks,theNThreads);
    }
  
    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeGen,start,stop);
    
    elapsedTimeGPU = elapsedTimeInit + elapsedTimeGen;
    
    printf("PNG Time on GPU : %6.3f %6.3f %6.3f ms\n",
	   elapsedTimeInit,elapsedTimeGen,elapsedTimeGPU);

    //copy results 
    cudaMemcpy(result_h,result_d,theNBlocks*theNThreads*sizeof(int),
	       cudaMemcpyDeviceToHost);
    
    printf("Executing the host code in CPU\n");
    
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    rand_cpu(result_c,nTracks);
    //<<<---------------------------------------------------------->>>
    
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    
    float elapsedTimeCPU = 0.0;
    cudaEventElapsedTime (&elapsedTimeCPU,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTimeCPU);
    
    printf("Ratio of Time Elapsed on CPU/GPU : %5.2f %5.2f \n",
	   elapsedTimeCPU/elapsedTimeGPU,elapsedTimeCPU/elapsedTimeGen);
    
    //clean up: destory cuda event and free memory on device and host
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaFree(devStates);
    cudaFree(devStatesMtgp32);
    cudaFree(devStatesMRG32k3a);

    cudaFree(result_d);
    free(result_h);
    free(result_c);
  }
  //end of event-loop

}
