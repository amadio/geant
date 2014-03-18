#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "EMPhysics_kernel.h"
#include "random_kernel.h"

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument: kernel type 
  int emModel  = 0; //0 em kernel, 1 brem, 2 ioni, 3 msc, 4 1-2-3-serial

  if(argc >= 2) emModel = atoi(argv[1]);

  if(emModel < 0 || emModel > 4) {
    std::cout << "Usage: EMPhysics mode=[0-4]: ... " << std::endl;
    return 0;
  }

  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    std::cout << "Waning: No Cuda Capable Device ... " << std::endl;
  }

  // set the default number of threads and thread blocks
  int theNBlocks  =  32;
  int theNThreads = 128;

  // prepare physics table
  bool useSpline = true;

  GPPhysicsTable eBrem_table;
  readTable(&eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eBrem_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_table;
  readTable(&eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable msc_table;
  readTable(&msc_table,"data/Lambda.msc.e-.asc");
  nv = msc_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    msc_table.physicsVectors[j].SetSpline(useSpline);
  }

  printf("Copying Physics data tables from host to device\n");

  GPPhysicsTable* eBrem_table_d;
  GPPhysicsTable* eIoni_table_d;
  GPPhysicsTable* msc_table_d;

  if(cudaEnabled) {
    cudaMalloc((void**)&eBrem_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eBrem_table_d, &eBrem_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_table_d, &eIoni_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&msc_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(msc_table_d, &msc_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);
  }
  
  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  //event-loop
  int nevent = 10;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 100000;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    // populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);

    //repeat on the CPU side, but with malloc
    GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));

    //randomize position and momenta    
    for(size_t i = 0 ; i < nTracks ; i++){
      track_h[i].x   = track_c[i].x   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].y   = track_c[i].y   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].z   = track_c[i].z   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].s   = track_c[i].s   = 10*(1.0*rand()/RAND_MAX);
      track_h[i].px  = track_c[i].px  = 10.+2000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].py  = track_c[i].py  = 10.+2000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].pz  = track_c[i].pz  = 10.+2000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].q   = track_c[i].q   = -1.0;
    }
    
    float elapsedTimeAlloc = 0.0;
    float elapsedTimeUp = 0.0;
    float elapsedTimeKernel = 0.0;
    float elapsedTimeDown = 0.0;
    float elapsedTimeGPU = 0.0;
    
    printf("Copying track data from device to host\n");
    
    //start time memory allocation on GPU
    cudaEvent_t start, stop;
    cudaEventCreate (&start);
    cudaEventCreate (&stop);

    //start time for allocation
    cudaEventRecord (start,0);
    
    GXTrack *track_d;
    cudaMalloc((void**)&track_d, nTracks*sizeof(GXTrack));
    
    //stop time for cudaMalloc
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
    
    //start time for Memcpy (H2D)
    cudaEventRecord (start,0);

    cudaMemcpy(track_d, track_h, nTracks*sizeof(GXTrack), 
	       cudaMemcpyHostToDevice);
    
    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeUp,start,stop);

    //prepare random engines on the device
    curandState* devStates = 0;
    cudaMalloc(&devStates, theNBlocks*theNThreads* sizeof(curandState));
    curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);

    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(cudaEnabled) {
      if(emModel == 0 ) {
	EMPhysics_gpu(devStates, track_d, nTracks,
		      eBrem_table_d, eIoni_table_d, msc_table_d,
		      true,true,
		      theNBlocks,theNThreads);
      }
      else if (emModel == 1) {
	brem_gpu(devStates, track_d, nTracks, eBrem_table_d, true, true,
		 theNBlocks, theNThreads);
      }
      else if (emModel == 2) {
	ioni_gpu(devStates, track_d, nTracks,eIoni_table_d, true, true,
		 theNBlocks, theNThreads);
      }
      else if (emModel == 3) {
	msc_gpu(devStates, track_d, nTracks, msc_table_d, true, true,
		theNBlocks, theNThreads);
      }
      else if (emModel == 4) {
	brem_gpu(devStates, track_d, nTracks, eBrem_table_d, true, true,
		 theNBlocks, theNThreads);
	ioni_gpu(devStates, track_d, nTracks, eIoni_table_d, true, true,
		 theNBlocks, theNThreads);
	msc_gpu(devStates, track_d, nTracks, msc_table_d, true, true,
		theNBlocks, theNThreads);
	//the physics process with the shortest proposed step length will be
	//chosen as the process governing the step, which requires separated 
	//kernel calls for PostStep and AlongStep process of each physics 
	//model - implementation details is not attempted for this test
      }
    }
    //<<<---------------------------------------------------------->>>
  
    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);
    
    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
    if(cudaEnabled) {
    	cudaMemcpy(track_h, track_d,nTracks*sizeof(GXTrack), 
		   cudaMemcpyDeviceToHost);
    }

    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeDown,start,stop);
    
    elapsedTimeGPU = elapsedTimeUp + elapsedTimeKernel + elapsedTimeDown;
    
    printf("Time Elapsed on GPU : %6.3f %6.3f %6.3f %6.3f %6.3f ms\n",
	   elapsedTimeAlloc, elapsedTimeUp,elapsedTimeKernel,
	   elapsedTimeDown,elapsedTimeGPU);
    
    printf("Executing the host code in CPU\n");
    
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(emModel == 0 ) {
      EMPhysics_cpu(track_c,nTracks,
		    &eBrem_table,&eIoni_table,&msc_table,
		    true,true);
    }
    else if (emModel == 1) {
      brem_cpu(track_c,nTracks, &eBrem_table, true,true);
    }
    else if (emModel == 2) {
      ioni_cpu(track_c,nTracks, &eIoni_table, true,true);
    }
    else if (emModel == 3) {
      msc_cpu(track_c,nTracks, &msc_table, true,true);
    }
    else if (emModel == 4) {
      brem_cpu(track_c,nTracks, &eBrem_table, true,true);
      ioni_cpu(track_c,nTracks, &eIoni_table, true,true);
      msc_cpu(track_c,nTracks, &msc_table, true,true);
    }
    //<<<---------------------------------------------------------->>>
    
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    
    float elapsedTime2 = 0.0;
    cudaEventElapsedTime (&elapsedTime2,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTime2);
    
    printf("Ratio of Time Elapsed on CPU/GPU : %5.2f %5.2f \n",
	   elapsedTime2/elapsedTimeGPU,elapsedTime2/elapsedTimeKernel);
    
    //clean up: destory cuda event and free memory on device and host
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    if(cudaEnabled) {
      cudaFree(track_d);
      cudaFree(devStates);
    }
    cudaFreeHost(track_h);
    free(track_c);
  }

  //end of event-loop
  if(cudaEnabled) {
    cudaFree(eBrem_table_d);
    cudaFree(eIoni_table_d);
    cudaFree(msc_table_d);
  }

}

