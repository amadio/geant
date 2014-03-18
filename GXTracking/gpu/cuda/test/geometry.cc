#include <iostream>
#include <fstream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "geometry_kernel.h"
#include "GXTrack.h"

//Magnetic FieldMap
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPSimpleCMS.h"
#include "GPVPhysicalVolume.h"

using namespace std;

int main (int argc, char* argv[]) 
{
  //argument
  int kernelType = 1; // 1. Navigator,  
                      // 2. MultiLevelLocator
  int geomType = 1;   // 1. SimpleCalo, 2. SimpleCMS

  if(argc >= 2) kernelType = atoi(argv[1]);
  if(argc >= 3) geomType = atoi(argv[2]);

  if(kernelType>2 || geomType>2 ) {
    std::cout << "Usage: geometry [kernelType=1-2] [geomType=1-2] ..." 
	      << std::endl;
    return 0;
  }

  //check whether device overlap is available
  int whichDevice;
  cudaGetDevice(&whichDevice);
  if(whichDevice > 0) {
    cudaDeviceReset();
  }

  //default number of threads and blocks
  int theNBlocks  =  32;
  int theNThreads = 128;

  char* cudaNBlocks;
  char* cudaNThreads;

  cudaNBlocks = getenv("GP_CUDA_NBLOCKS");
  if(cudaNBlocks) theNBlocks = atoi(cudaNBlocks);
  cudaNThreads = getenv("GP_CUDA_NTHREADS");
  if(cudaNThreads) theNThreads = atoi(cudaNThreads);
  
  std::cout << "...transporation_kernel<<<" << theNBlocks << "," 
	    << theNThreads <<">>> (...) ..." << std::endl;
  
  //1. Constructor geometry
  GPVGeometry *geom;

  if(geomType == 1) {
    int nphi =  4;
    int nz   =  3;
    double density = 8.28;

    geom = new GPSimpleEcal(nphi,nz,density);
  }
  else if(geomType == 2) {
    geom = new GPSimpleCMS();
  }

  geom->create();
 
  GPVGeometry::byte *geom_h = (GPVGeometry::byte *) malloc (geom->size()) ;
  geom->relocate( geom_h );
  memcpy(geom_h,geom->getBuffer(),geom->size());  

  GPVGeometry::byte *geom_d;
  cudaMalloc( (void**)&geom_d, geom->size() );
  geom->relocate( geom_d );
  cudaMemcpy(geom_d, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice);

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
    
    std::cout << "... chordfinder ... Loading magnetic field map: " 
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

  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 100000;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    //3. populate tracks with track position, momentum, energy, steplength
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
      track_h[i].s   = track_c[i].s   =  1.0+1*(2.0*rand()/RAND_MAX-1.0);

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
    
    printf("Size of Track = %d and GPNavigator = %d\n",sizeof(GXTrack),sizeof(GPNavigator));
    //start time for kernel
    cudaEventRecord (start,0);
    
   //<<<---------------------------------------------------------->>>
    if(kernelType == 1) {
      navigator_gpu((GPGeomManager*)geom_d,bmap_d,track_d,nTracks,
		     theNBlocks,theNThreads);
    }
    else if(kernelType == 2) {
      mllocator_gpu((GPGeomManager*)geom_d,bmap_d,track_d,nTracks,
		    theNBlocks,theNThreads);
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
    if(kernelType == 1) {
      navigator_cpu((GPGeomManager*)geom_h,bmap_h,track_c,nTracks);
    }
    else if(kernelType == 2) {
      mllocator_cpu((GPGeomManager*)geom_h,bmap_h,track_c,nTracks);
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
    cudaFree(track_d);
    cudaFreeHost(track_h);
    free(track_c);
  }
  //end of event-loop
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

}
