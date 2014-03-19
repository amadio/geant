#include <iostream>
#include <fstream>

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>

#include "tracking_kernel.h"

#include "GPConstants.h"
#include "GPFieldMap.h"
#include "GPFieldMapData.h"

//Geometry
#include "GPTrack.h"
#include "GPThreeVector.h"

#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPSimpleCMS.h"
#include "GPVPhysicalVolume.h"

#include "GPConstants.h"
#include "EMPhysics_kernel.h"

using namespace std;

int main () 
{

  //GPThreeVector a = GPThreeVector_create(0.,0.,0.);

  int nDevice;
  cudaGetDeviceCount(&nDevice);
  bool cudaEnabled = false;

  int theNBlocks  =  32;
  int theNThreads = 128;
  // int theNBlocks  =  2;
  // int theNThreads = 8;

  if(nDevice > 0) {
	  cudaDeviceReset();
	  cudaEnabled = true;
  }

  //1. Constructor geometry
  int nphi = 4;
  int nz   = 3;
  double density = 8.28;

  GPVGeometry *geom = new GPSimpleEcal(nphi,nz,density);

  geom->create();
 
  GPVGeometry::byte *geom_h = (GPVGeometry::byte *) malloc (geom->size()) ;
  geom->relocate( geom_h );
  memcpy(geom_h,geom->getBuffer(),geom->size());  

  GPVGeometry::byte *geom_d;
  cudaMalloc( (void**)&geom_d, geom->size() );
  geom->relocate( geom_d );
  cudaMemcpy(geom_d, geom->getBuffer(), geom->size(), cudaMemcpyHostToDevice);

  //2. Read magnetic field map

  GPFieldMap *bmap_d;
  GPFieldMap *bmap_h;
  GPFieldMap** fieldMap;

  fieldMap = (GPFieldMap **) malloc (nbinZ*sizeof(GPFieldMap *));
  for (int j = 0 ; j < nbinZ ; j++) {
    fieldMap[j] = (GPFieldMap *) malloc (nbinR*sizeof(GPFieldMap));
  } 

  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";

  std::ifstream ifile(fieldMapFile, ios::in | ios::binary | ios::ate);
  
  if (ifile.is_open()) {

    //field map structure
    GPFieldMapData fd;

    ifstream::pos_type fsize = ifile.tellg();
    size_t dsize = sizeof(GPFieldMapData);    

    long int ngrid = fsize/dsize;
    ifile.seekg (0, ios::beg);
    
    std::cout << "... transportation ... Loading magnetic field map: " 
	      << fieldMapFile << std::endl;

    for(int i = 0 ; i < ngrid ; i++) {
      ifile.read((char *)&fd, sizeof(GPFieldMapData));
      
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
  bmap_h = (GPFieldMap *) malloc (nbinZ*nbinR*sizeof(GPFieldMap));

  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
      bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
    }
  }

  if(cudaEnabled) {
	  cudaMalloc((void**)&bmap_d, nbinZ*nbinR*sizeof(GPFieldMap));
	  cudaMemcpy(bmap_d,bmap_h,nbinZ*nbinR*sizeof(GPFieldMap),cudaMemcpyHostToDevice);
  }
   
  // 4 prepare physics table
  bool useSpline = true;

  GPPhysicsTable eBrem_table;
  readTable(&eBrem_table,"data/Lambda.eBrem.e-.asc");
  int nv = eBrem_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eBrem_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eBrem_table);

  GPPhysicsTable eIoni_table;
  readTable(&eIoni_table,"data/Lambda.eIoni.e-.asc");
  nv = eIoni_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    eIoni_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&eIoni_table);

  GPPhysicsTable msc_table;
  readTable(&msc_table,"data/Lambda.msc.e-.asc");
  nv = msc_table.nPhysicsVector;
  for(int j=0; j < nv; j++){
    msc_table.physicsVectors[j].SetSpline(useSpline);
  }
  //  GPPhysicsTable_Print(&msc_table);


  printf("Copying data from host to device\n");

  GPPhysicsTable* eBrem_table_d;
  GPPhysicsTable* eIoni_table_d;
  GPPhysicsTable* msc_table_d;

  if(cudaEnabled) {
	  cudaMalloc((void**)&eBrem_table_d, sizeof(GPPhysicsTable));
	  cudaMemcpy(eBrem_table_d, &eBrem_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);

	  cudaMalloc((void**)&eIoni_table_d, sizeof(GPPhysicsTable));
	  cudaMemcpy(eIoni_table_d, &eIoni_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);

	  cudaMalloc((void**)&msc_table_d, sizeof(GPPhysicsTable));
	  cudaMemcpy(msc_table_d, &msc_table, sizeof(GPPhysicsTable),cudaMemcpyHostToDevice);
  }
  
  //5. get the list of tracks from the stack
  printf("Preparing a bundle tracks at host\n");

  //event-loop
  int nevent = 2;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 5*32*128;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    //3. populate tracks with track position, momentum, energy, steplength
    GPTrack *track_h = (GPTrack *) malloc (nTracks*sizeof(GPTrack));
    //repeat on the CPU side
    
    GPTrack *track_c = (GPTrack *) malloc (nTracks*sizeof(GPTrack));

    //randomize position and momenta    
    for(size_t i = 0 ; i < nTracks ; i++){
      //barrel only
      track_h[i].x      = track_c[i].x     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].y      = track_c[i].y     = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].z      = track_c[i].z     = 300*(2.0*rand()/RAND_MAX-1.0);

      track_h[i].px     = track_c[i].px    = 10.+2000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].py     = track_c[i].py    = 10.+2000*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].pz     = track_c[i].pz    = 10.+2000*(2.0*rand()/RAND_MAX-1.0);

      track_h[i].E      = track_c[i].E     = 
	sqrt(track_h[i].px*track_h[i].px + track_h[i].py*track_h[i].py + 
      	     track_h[i].pz*track_h[i].pz + 140*140);

      track_h[i].q      = track_c[i].q      = -1.0;
      track_h[i].step   = track_c[i].step   =  1.0;
      track_h[i].length = track_c[i].length = 10.0;
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
    
    GPTrack *track_d;
    cudaMalloc((void**)&track_d, nTracks*sizeof(GPTrack));
    
    //stop time for cudaMalloc
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
    
    //start time for Memcpy (H2D)
    cudaEventRecord (start,0);

    cudaMemcpy(track_d, track_h, nTracks*sizeof(GPTrack), 
	       cudaMemcpyHostToDevice);
    
    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeUp,start,stop);
    
    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<---------------------------------------------------------->>>
    if(cudaEnabled) {
    	tracking_gpu((GPVPhysicalVolume*)geom_d,bmap_d,
		 eBrem_table_d, eIoni_table_d, msc_table_d,
		 track_d,nTracks,theNBlocks,theNThreads);
    }
    //<<<---------------------------------------------------------->>>
  
    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);
    
    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
    if(cudaEnabled) {
    	cudaMemcpy(track_h, track_d,nTracks*sizeof(GPTrack), cudaMemcpyDeviceToHost);
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
    tracking_cpu((GPVPhysicalVolume*)geom_h,bmap_h,
		 &eBrem_table,&eIoni_table,&msc_table,
		 track_c,nTracks);
    //<<<---------------------------------------------------------->>>
    
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    
    float elapsedTime2 = 0.0;
    cudaEventElapsedTime (&elapsedTime2,start,stop);
    printf("Time Elapsed on CPU : %6.3f ms\n",elapsedTime2);
    
    printf("Ratio of Time Elapsed on CPU/GPU : %5.2f %5.2f \n",
	   elapsedTime2/elapsedTimeGPU,elapsedTime2/elapsedTimeKernel);
    
    //Compare results from GPU and CPU
    long int miscount = 0;
    for (size_t i = 0 ; i < nTracks ; i++) {
      if( track_h[i].x != track_c[i].x || track_h[i].y != track_c[i].y || 
	  track_h[i].z != track_c[i].z ) {
    	  std::cout << "Compare results: check track position from gpu cpu "
    			  << "(" << track_h[i].x << "," << track_h[i].y << ","
    			  <<track_h[i].z << ")" << " (" << track_c[i].x << ","
    			  << track_c[i].y <<","<<track_c[i].z << ")" << std::endl;
    	  miscount++;
      }
    }
    printf("Number of Cuda MisMatched Tracks = %ld\n",miscount);

    //clean up: destory cuda event and free memory on device and host
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    if(cudaEnabled) cudaFree(track_d);
    free(track_h);
    free(track_c);
  }

  //end of event-loop
  if(cudaEnabled) {
	  cudaFree(geom_d);
	  cudaFree(bmap_d);
	  cudaFree(eBrem_table_d);
	  cudaFree(eIoni_table_d);
	  cudaFree(msc_table_d);
  }
  free(geom_h);
  free(bmap_h);

  std::cout << "Done!" << std::endl;
}

