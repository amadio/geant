#include <iostream>
#include <fstream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "stepper_kernel.h"

#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument
  int stepperType = 1; //1 rk45, 2 rkf45, 3 nrk4
  if(argc >= 2) stepperType = atoi(argv[1]);

  if(stepperType < 1 || stepperType > 4) {
    std::cout << "Usage: stepper [1|2|3|4] ... " <<
      "1=RK4, 2=Felhberg, 3=Nystrom, 4=Dummy" << std::endl;
    return 0;
  }

  cudaDeviceReset();

  //default threads and blocks
  int theNBlocks  =  32;
  int theNThreads = 128;

  char* cudaNBlocks;
  char* cudaNThreads;

  cudaNBlocks = getenv("GP_CUDA_NBLOCKS");
  if(cudaNBlocks) theNBlocks = atoi(cudaNBlocks);
  cudaNThreads = getenv("GP_CUDA_NTHREADS");
  if(cudaNThreads) theNThreads = atoi(cudaNThreads);

  std::cout << "...stepper_kernel<<<" << theNBlocks << "," 
	    << theNThreads <<">>> (...) ..." << std::endl;

  //2. Read magnetic field map

  GXFieldMap *bmap_d;
  GXFieldMap *bmap_h;
  GXFieldMap** fieldMap;

  fieldMap = (GXFieldMap **) malloc (nbinZ*sizeof(GXFieldMap *));
  for (int j = 0 ; j < nbinZ ; j++) {
    fieldMap[j] = (GXFieldMap *) malloc (nbinR*sizeof(GXFieldMap));
  } 

  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";

  std::ifstream ifile(fieldMapFile, ios::in | ios::binary | ios::ate);
  
  if (ifile.is_open()) {

    //field map structure
    GXFieldMapData fd;

    ifstream::pos_type fsize = ifile.tellg();
    size_t dsize = sizeof(GXFieldMapData);    

    long int ngrid = fsize/dsize;
    ifile.seekg (0, ios::beg);
    
    std::cout << "... transportation ... Loading magnetic field map: " 
	      << fieldMapFile << std::endl;

    for(int i = 0 ; i < ngrid ; i++) {
      ifile.read((char *)&fd, sizeof(GXFieldMapData));
      
      //check validity of input data
      if(abs(fd.iz) > noffZ || fd.ir > nbinR) {
	std::cout << " Field Map Array Out of Range" << std::endl;
      }
      else {
        fieldMap[noffZ+fd.iz][fd.ir].Bz = fd.Bz; 
        fieldMap[noffZ+fd.iz][fd.ir].Br = fd.Br;
      }
    }
    ifile.close();
  }

  //3. create magnetic field on the device

  printf("Creating magnetic field map on the GPU device\n");
  //prepare fieldMap array: fieldMap[nbinZ][nbinR];
  
  //cuda enabled
  bmap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));

  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
      bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
    }
  }
  
  printf("Copying data from host to device\n");
  
  cudaMalloc((void**)&bmap_d, nbinZ*nbinR*sizeof(GXFieldMap)) ;
  cudaMemcpy(bmap_d,bmap_h,nbinZ*nbinR*sizeof(GXFieldMap),
	     cudaMemcpyHostToDevice);
   
  //4. get the list of tracks from the stack
  printf("Preparing a bundle tracks at host\n");

  //event-loop
  int nevent = 10;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 100000;

    printf("Event %d> stepper=%d for %d tracks\n",i,stepperType,nTracks);

    //3. populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
		  cudaHostAllocDefault);

    //repeat on the CPU side, but with malloc
    GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));

    //randomize position and momenta    
    for(size_t i = 0 ; i < nTracks ; i++){
      track_h[i].x      = track_c[i].x     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].y      = track_c[i].y     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].z      = track_c[i].z     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].s      = track_c[i].s     = 1.+100*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].px     = track_c[i].px    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].py     = track_c[i].py    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].pz     = track_c[i].pz    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].q      = track_c[i].q     = -1.0;
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
    
    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(stepperType == 1) {
      rk4_gpu(bmap_d,track_d,nTracks,theNBlocks,theNThreads);
    }
    else if(stepperType == 2) {
      rkf45_gpu(bmap_d,track_d,nTracks,theNBlocks,theNThreads);
    }
    else if(stepperType == 3) {
      nrk4_gpu(bmap_d,track_d,nTracks,theNBlocks,theNThreads);
    }
    //<<<---------------------------------------------------------->>>
  
    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);
    
    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
    cudaMemcpy(track_h, track_d,nTracks*sizeof(GXTrack), 
	       cudaMemcpyDeviceToHost);
  
    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeDown,start,stop);
    
    elapsedTimeGPU= elapsedTimeUp + elapsedTimeKernel + elapsedTimeDown;
    
    printf("Time Elapsed on GPU : %6.3f %6.3f %6.3f %6.3f %6.3f ms\n",
	   elapsedTimeAlloc, elapsedTimeUp,elapsedTimeKernel,
	   elapsedTimeDown,elapsedTimeGPU);
    
    printf("Executing the host code in CPU\n");
    
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(stepperType == 1)         rk4_cpu(bmap_h,track_c,nTracks);
    else if (stepperType == 2) rkf45_cpu(bmap_h,track_c,nTracks);
    else if (stepperType == 3)  nrk4_cpu(bmap_h,track_c,nTracks);
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
    
    cudaFree(track_d);
    cudaFreeHost(track_h);
    free(track_c);
  }
  //end of event-loop

  //clean up: free memory on device and host
  cudaFree(bmap_d);
  free(bmap_h);
}
