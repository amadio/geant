#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GPEMPhysicsUtils.h"
#include "electronTest_kernel.h"
#include "random_kernel.h"
#include "dma_kernel.h"

#include "GPPhysics2DVector.h"

#include <fstream>
#include <iomanip>

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument: kernel type 
  int emModel = 0; //0 elctron 1 brem, 2 ioni, 3 msc
                   //21 brem, 22 ioni to stack secondaries with a dynamic memory
  int runType = 0; //0 normal, 1 test
  int isStack = 0; //0 stack secondaries 1 no secondary stack

  if(argc >= 2) emModel = atoi(argv[1]);
  if(argc >= 3) runType = atoi(argv[2]);
  if(argc >= 4) isStack = atoi(argv[3]);

  if(emModel < 0 || emModel > 40) {
    std::cout << "Usage: electronTest model=[0-3,21-22] [0-1] [0-1]: ... " 
	      << std::endl;
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

  GPPhysicsTable eBrem_table;
  readTable(&eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eBrem_table.physicsVectors[j].SetSpline(useSpline);
  }

  //Ionisation

  GPPhysicsTable eIoni_table;
  readTable(&eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_table.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_range;
  readTable(&eIoni_range,"data/Range.eIoni.e-.asc");
  nv = eIoni_range.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_range.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_dedx;
  readTable(&eIoni_dedx,"data/DEDX.eIoni.e-.asc");
  nv = eIoni_dedx.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_dedx.physicsVectors[j].SetSpline(useSpline);
  }

  GPPhysicsTable eIoni_invr;
  readTable(&eIoni_invr,"data/InverseRange.eIoni.e-.asc");
  nv = eIoni_invr.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_invr.physicsVectors[j].SetSpline(useSpline);
  }

  //msc

  GPPhysicsTable msc_table;
  readTable(&msc_table,"data/Lambda.msc.e-.asc");
  nv = msc_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    msc_table.physicsVectors[j].SetSpline(useSpline);
  }

  printf("Copying Physics data tables from host to device\n");

  GPPhysicsTable* eBrem_table_d;

  GPPhysicsTable* eIoni_table_d;
  GPPhysicsTable* eIoni_range_d;
  GPPhysicsTable* eIoni_dedx_d;
  GPPhysicsTable* eIoni_invr_d;

  GPPhysicsTable* msc_table_d;

  if(cudaEnabled) {
    cudaMalloc((void**)&eBrem_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eBrem_table_d, &eBrem_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_table_d, &eIoni_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_range_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_range_d, &eIoni_range, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_dedx_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_dedx_d, &eIoni_dedx, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&eIoni_invr_d, sizeof(GPPhysicsTable));
    cudaMemcpy(eIoni_invr_d, &eIoni_invr, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&msc_table_d, sizeof(GPPhysicsTable));
    cudaMemcpy(msc_table_d, &msc_table, sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);
  }
  
  //G4SeltzerBergerModel data
  G4int maxZ = 92;
  GPPhysics2DVector* sbData =
    (GPPhysics2DVector*) malloc(maxZ*sizeof(GPPhysics2DVector));

  char sbDataFile[256];
  for(G4int iZ = 0 ; iZ < maxZ ; iZ++) {  
    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
    if(!check) {
      printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }

  GPPhysics2DVector* sbData_d;
  cudaMalloc((void**)&sbData_d,maxZ*sizeof(GPPhysics2DVector));
  cudaMemcpy(sbData_d, sbData, maxZ*sizeof(GPPhysics2DVector),
	     cudaMemcpyHostToDevice);

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

    G4double photonFraction = 0.0;
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
      track_h[i].s = 10.0*Rndm(); 

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

    //fixed size of memory for secondaries on host
    GXTrack *secTracks_h
      = (GXTrack*) malloc(maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(cudaEnabled) {
      if(emModel == 0 ) {

        electron_gpu(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d, 
		     eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		     msc_table_d, sbData_d,
		     secTracks_d,stackSize_d,offset_d,isStack,runType,
		     theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }         
      }
      else if(emModel == 1 ) {

	brem_gpu(devStates, track_d, nTracks,
		 eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);

        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }         
      }
      else if(emModel == 2 ) {

	ioni_gpu(devStates, track_d, nTracks,
		 eBrem_table_d, eIoni_table_d, 
		 eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		 msc_table_d, sbData_d,
		 secTracks_d,stackSize_d,offset_d,isStack,runType,
		 theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }       
      }
      else if(emModel == 3 ) {

	msc_gpu(devStates, track_d, nTracks,
		eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		secTracks_d,stackSize_d,offset_d,isStack,runType,
		theNBlocks,theNThreads);
        if(isStack==0) {
          cudaThreadSynchronize();
          cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                     cudaMemcpyDeviceToHost);
        }         
      }
      //<<<---------------------------------------------------------->>>
      else if(emModel == 11 ) {

	//memory allocation for a containter holding secondaries per block 
	GXSecContainer *secContainer_d;
	cudaMalloc((void**)&secContainer_d, 
		   theNBlocks*sizeof(GXSecContainer));

	brem_gpu_dma(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d, msc_table_d, sbData_d,
		     secContainer_d,stackSize_d,offset_d,
		     theNBlocks,theNThreads);

	cudaThreadSynchronize();

	//get the number of secondaries created on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		   cudaMemcpyDeviceToHost);

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
	dma_free_perblock_gpu(nTracks,secContainer_d,theNBlocks,1);
	cudaFree(secContainer_d);

      }
      else if(emModel == 12 ) {

	//memory allocation for a containter holding secondaries per block 
	GXSecContainer *secContainer_d;
	cudaMalloc((void**)&secContainer_d, 
		   theNBlocks*sizeof(GXSecContainer));

	ioni_gpu_dma(devStates, track_d, nTracks,
		     eBrem_table_d, eIoni_table_d, 
		     eIoni_range_d, eIoni_dedx_d, eIoni_invr_d, 
		     msc_table_d, sbData_d,
		     secContainer_d,stackSize_d,offset_d,
		     theNBlocks,theNThreads);

	cudaThreadSynchronize();

	//get the number of secondaries created on GPU
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		   cudaMemcpyDeviceToHost);

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
	dma_free_perblock_gpu(nTracks,secContainer_d,theNBlocks,1);
	cudaFree(secContainer_d);

      }
    }

    cudaThreadSynchronize();
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

    //reset stackSize for CPU
    stackSize = 0;
   
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(emModel == 0 ) {
      electron_cpu(track_c,nTracks,
		   &eBrem_table,&eIoni_table,
		   &eIoni_range, &eIoni_dedx, &eIoni_invr,&msc_table,
		   sbData,secTracks_h,&stackSize,
		   isStack,runType);
    }
    else if(emModel == 1 ) {
      brem_cpu(track_c,nTracks,
	       &eBrem_table,&eIoni_table,&msc_table,sbData,
	       secTracks_h,&stackSize,
	       isStack,runType);
    }
    else if(emModel == 2 ) {
      ioni_cpu(track_c,nTracks,
	       &eBrem_table,&eIoni_table,
	       &eIoni_range, &eIoni_dedx, &eIoni_invr, &msc_table,
	       sbData,secTracks_h,&stackSize,
	       isStack,runType);
    }
    else if(emModel == 3 ) {
      msc_cpu(track_c,nTracks,
	      &eBrem_table,&eIoni_table,&msc_table,sbData,
	      secTracks_h,&stackSize,isStack,runType);
    }
    //<<<---------------------------------------------------------->>>
    else if(emModel == 11 ) {
      brem_cpu_dma(track_c,nTracks,
	       &eBrem_table,&eIoni_table,&msc_table,sbData);
    }
    else if(emModel == 12 ) {
      ioni_cpu_dma(track_c,nTracks,
	       &eBrem_table,&eIoni_table,
	       &eIoni_range, &eIoni_dedx, &eIoni_invr, 
	       &msc_table,sbData);
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
    cudaFree(eBrem_table_d);
    cudaFree(eIoni_table_d);
    cudaFree(eIoni_range_d);
    cudaFree(eIoni_dedx_d);
    cudaFree(eIoni_invr_d);
    cudaFree(msc_table_d);
  }

  free(sbData);
  cudaFree(sbData_d);
 

}

