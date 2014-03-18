#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GPEMPhysicsUtils.h"
#include "photonTest_kernel.h"
#include "random_kernel.h"

#include <fstream>
#include <iomanip>

using namespace std;

int main (int argc, char* argv[]) 
{
  //arguments: [kernel type] [run type] [store stack] 
  int emModel  = 0; //0 gamma kernel, 1 compton, 2 conversion, 3 phot-electric
  int runType  = 0; //0 normal 1 test
  int isStack  = 0; //0 stack secondaries 1 no secondary stack

  if(argc >= 2) emModel = atoi(argv[1]);
  if(argc >= 3) runType = atoi(argv[2]);
  if(argc >= 4) isStack = atoi(argv[3]);

  if(emModel < 0 || emModel > 3) {
    std::cout << "Usage: GammaPhysics mode=[0-3]: ... " << std::endl;
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
  if(runType == 1 ) { 
    theNBlocks  =  2;
    theNThreads =  4;
  }

  //@@@pre-allocated memory for Secondaries
  //maximum size of dynamic memory allocation on the device
  //set 256 megabytes on the heap (global mememory) 
  cudaThreadSetLimit(cudaLimitMallocHeapSize, 256*1024*1024);

  // prepare physics table
  bool useSpline = true;

  GPPhysicsTable gCompt_table;
  readTable(&gCompt_table,"data/Lambda.compt.gamma.asc");
  G4int nv = gCompt_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    gCompt_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable gConv_table;
  readTable(&gConv_table,"data/Lambda.conv.gamma.asc");
  nv = gConv_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    gConv_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable gPhot_table;
  readTable(&gPhot_table,"data/LambdaPrim.phot.gamma.asc");
  nv = gPhot_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    gPhot_table.physicsVectors[j].SetSpline(useSpline);
  }

  printf("Copying Physics data tables from host to device\n");

  GPPhysicsTable* gCompt_table_d;
  GPPhysicsTable* gConv_table_d;
  GPPhysicsTable* gPhot_table_d;

  if(cudaEnabled) {
    cudaMalloc((void**)&gCompt_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(gCompt_table_d, &gCompt_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&gConv_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(gConv_table_d, &gConv_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&gPhot_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(gPhot_table_d, &gPhot_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);
  }
  
  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  //event-loop
  int nevent = 10;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 20000;
    if(runType == 1 ) nTracks = 16;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    // populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);

    //repeat on the CPU side, but with malloc
    GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));

    //generate tracks with random position and momentum
    G4double photonFraction = 1.0;
    const G4double minP = 20.;
    const G4double maxP = 1000.;

    G4double rho, z, p, mass;
    G4double theta, phi;

    for(size_t i = 0 ; i < nTracks ; i++){

      rho = ecalRmim + (ecalRmax-ecalRmim)*Rndm();
      p = minP + (maxP - minP)*Rndm();
      z = ecalZmax*(2*Rndm()-1.0);
      phi = twopi*Rndm();
      theta = std::atan(rho/z);

      track_h[i].status = 0;
      track_h[i].q = ( Rndm() < photonFraction ) ? 0.0 : -1.0;
      track_h[i].x = rho*std::cos(phi);
      track_h[i].y = rho*std::sin(phi); 
      track_h[i].z = z; 
      track_h[i].s = 10.*Rndm(); 

      track_h[i].px = p*std::sin(theta)*std::cos(phi);
      track_h[i].py = p*std::sin(theta)*std::sin(phi);
      track_h[i].pz = p*std::cos(theta);

      mass = electron_mass_c2*track_h[i].q*track_h[i].q;
      track_h[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);

      CopyTrack(&track_h[i],&track_c[i]);
    }

    float elapsedTimeAlloc = 0.0;
    float elapsedTimeUp = 0.0;
    float elapsedTimeKernel = 0.0;
    float elapsedTimeDown = 0.0;
    float elapsedTimeGPU = 0.0;
    
    printf("Copying track data from host to device\n");
    
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

    //@@@pre-allocated memory for secondaries
    G4int stackSize = 0;
    GXTrack *secTracks_d;

    cudaMalloc((void**)&secTracks_d, 
	       maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //atomic counter for the total number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    //atomic counter for the last array position of scondaries
    G4int *offset_d;
    cudaMalloc((void**)&offset_d, sizeof(G4int));
    cudaMemset(offset_d,0,sizeof(G4int));

    GXTrack *secTracks_h
      = (GXTrack*) malloc(maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(cudaEnabled) {
      if(emModel == 0 ) {
	gamma_gpu(devStates, track_d, nTracks,
		  gCompt_table_d, gConv_table_d, gPhot_table_d,
		  secTracks_d,stackSize_d,offset_d,isStack,runType,
		  theNBlocks,theNThreads);

        cudaThreadSynchronize();
        cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                   cudaMemcpyDeviceToHost);

      }
      else if(emModel == 1 ) {
	compt_gpu(devStates, track_d, nTracks,
		  gCompt_table_d, gConv_table_d, gPhot_table_d,
		  secTracks_d,stackSize_d,offset_d,isStack,runType,
		  theNBlocks,theNThreads);

	if(isStack==0) {
	  cudaThreadSynchronize();
	  cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		     cudaMemcpyDeviceToHost);
	}
      }
      else if(emModel == 2 ) {
	conv_gpu(devStates, track_d, nTracks,
		 gCompt_table_d, gConv_table_d, gPhot_table_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);

	if(isStack==0) {
	  cudaThreadSynchronize();
	  cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		     cudaMemcpyDeviceToHost);
	}	  
      }
      else if(emModel == 3 ) {
	phot_gpu(devStates, track_d, nTracks,
		 gCompt_table_d, gConv_table_d, gPhot_table_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);

	if(isStack==0) {
	  cudaThreadSynchronize();
	  cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		     cudaMemcpyDeviceToHost);
	}
      }
    }
    //<<<---------------------------------------------------------->>>

    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);

    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
    printf("Copying track data from device to host\n");
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

    //reset stackSize for CPU
    stackSize = 0;
    
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(emModel == 0 ) {
      gamma_cpu(track_c,nTracks,
		&gCompt_table,&gConv_table,&gPhot_table,
		secTracks_h,&stackSize,isStack,runType);
    }
    else if(emModel == 1 ) {
      compt_cpu(track_c,nTracks,
		&gCompt_table,&gConv_table,&gPhot_table,
		secTracks_h,&stackSize,isStack,runType);
    }
    else if(emModel == 2 ) {
      conv_cpu(track_c,nTracks,
      	       &gCompt_table,&gConv_table,&gPhot_table,
	       secTracks_h,&stackSize,isStack,runType);
    }
    else if(emModel == 3 ) {
      phot_cpu(track_c,nTracks,
	       &gCompt_table,&gConv_table,&gPhot_table,
	       secTracks_h,&stackSize,isStack,runType);
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
    free(secTracks_h);

    cudaFree(stackSize_d);
    cudaFree(offset_d);
    cudaFree(secTracks_d);
  }

  //end of event-loop
  if(cudaEnabled) {
    cudaFree(gCompt_table_d);
    cudaFree(gConv_table_d);
    cudaFree(gPhot_table_d);
  }

}
