#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GPEMPhysicsUtils.h"
#include "GPPhysicsTableType.h"
#include "emTest_kernel.h"
#include "random_kernel.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPSimpleCMS.h"
#include "GPVPhysicalVolume.h"

//Magnetic Field
#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include <fstream>
#include <iomanip>

using namespace std;

int main (int argc, char* argv[]) 
{
  //arguments: [kernel type] [run type] [store stack] 
  int emModel  = 0; //0 photon kernel, 1 electron kernel
  int runType  = 0; //0 normal 1 test
  int numStep  = 1; //number of steps

  if(argc >= 2) emModel = atoi(argv[1]);
  if(argc >= 3) runType = atoi(argv[2]);
  if(argc >= 4) numStep = atoi(argv[3]);

  if(emModel < 0 || emModel > 3) {
    std::cout << "Usage: emTest mode=[0-1]: ... " << std::endl;
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

  //1. Construct geometry
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
  //3. Create magnetic field on the device
  
  printf("Creating magnetic field map on the GPU device\n");
  
  bmap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));
  
  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
      bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
    }
  }
  
  cudaMalloc((void**)&bmap_d, nbinZ*nbinR*sizeof(GXFieldMap)) ;
  cudaMemcpy(bmap_d,bmap_h,nbinZ*nbinR*sizeof(GXFieldMap),
             cudaMemcpyHostToDevice);


  //@@@pre-allocated memory for Secondaries
  //maximum size of dynamic memory allocation on the device
  //set 256 megabytes on the heap (global mememory) 
  cudaThreadSetLimit(cudaLimitMallocHeapSize, 256*1024*1024);

  // prepare physics table
  GPPhysicsTable physicsTable[kNumberPhysicsTable];

  char filename[256];
  for(int it = 0 ; it < kNumberPhysicsTable ; ++it) {
    sprintf(filename,"data/%s",GPPhysicsTableName[it]);
    readTableAndSetSpline(&physicsTable[it],filename);
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

  printf("Copying Physics tables from host to device\n");

  GPPhysicsTable* physicsTable_d;
  GPPhysics2DVector* sbData_d;

  if(cudaEnabled) {

    cudaMalloc((void**)&physicsTable_d, 
	       kNumberPhysicsTable*sizeof(GPPhysicsTable));
    cudaMemcpy(physicsTable_d, &physicsTable, 
	       kNumberPhysicsTable*sizeof(GPPhysicsTable),
	       cudaMemcpyHostToDevice);

    cudaMalloc((void**)&sbData_d,maxZ*sizeof(GPPhysics2DVector));
    cudaMemcpy(sbData_d, sbData, maxZ*sizeof(GPPhysics2DVector),
	       cudaMemcpyHostToDevice);

  }

  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  //event-loop
  int nevent = 10;
  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 4096*24;
    if(runType == 1 ) nTracks = 16;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    // populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);

    //repeat on the CPU side, but with malloc
    GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));

    //randomize position and momenta    
    for(size_t i = 0 ; i < nTracks ; i++){
      track_h[i].status   = track_c[i].status   = 0;
      track_h[i].x   = track_c[i].x   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].y   = track_c[i].y   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].z   = track_c[i].z   = 300*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].s   = track_c[i].s   = 10*(1.0*rand()/RAND_MAX);
      track_h[i].px  = track_c[i].px  = 10.+100*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].py  = track_c[i].py  = 10.+100*(2.0*rand()/RAND_MAX-1.0);
      track_h[i].pz  = track_c[i].pz  = 10.+100*(2.0*rand()/RAND_MAX-1.0);
      if(emModel == 0 ) {
	track_h[i].q   = track_c[i].q   = 0.0;
      }
      if(emModel == 1 ) {
	track_h[i].q   = track_c[i].q   = -1.0;
      }
      G4double p = sqrt(track_h[i].x*track_h[i].x +
			track_h[i].y*track_h[i].y + 
		        track_h[i].z*track_h[i].z);
      G4double mass = electron_mass_c2*track_h[i].q*track_h[i].q;
      G4double ekin = p*p/(sqrt(p*p + mass*mass) + mass);

      track_h[i].E  = track_c[i].E   = ekin;
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
	       numStep*maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //atomic counter for the total number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    //atomic counter for the last array position of scondaries
    G4int *offset_d;
    cudaMalloc((void**)&offset_d, sizeof(G4int));
    cudaMemset(offset_d,0,sizeof(G4int));

    GXTrack *secTracks_h
      = (GXTrack*) malloc(numStep*maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //start time for kernel
    cudaEventRecord (start,0);
    
    //<<<----------------------------------------------------------------->>>
    if(cudaEnabled) {
      if(emModel == 0 ) {
	photon_gpu(devStates, track_d, nTracks, 
		   (GPGeomManager*)geom_d, bmap_d, physicsTable_d,
		   secTracks_d, stackSize_d, offset_d, numStep, runType,
		   theNBlocks, theNThreads);

        cudaThreadSynchronize();
        cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
                   cudaMemcpyDeviceToHost);
      }
      else if(emModel == 1 ) {
	electron_gpu(devStates, track_d, nTracks, 
		     (GPGeomManager*)geom_d, bmap_d, physicsTable_d, sbData_d,
		     secTracks_d, stackSize_d, offset_d, numStep, runType,
		     theNBlocks, theNThreads);
	cudaThreadSynchronize();
	cudaMemcpy(&stackSize,stackSize_d,sizeof(G4int),
		   cudaMemcpyDeviceToHost);
      }         
    }
    //<<<----------------------------------------------------------------->>>

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
      photon_cpu(track_c, nTracks, 
		 (GPGeomManager*)geom_h, bmap_h, physicsTable, 
		 secTracks_h,&stackSize, numStep,runType);
    }
    else if(emModel == 1 ) {
      electron_cpu(track_c, nTracks, 
		   (GPGeomManager*)geom_h, bmap_h, physicsTable, sbData, 
		   secTracks_h, &stackSize, numStep, runType);
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

    cudaFree(secTracks_d);
    cudaFree(stackSize_d);
    cudaFree(offset_d);
    cudaFreeHost(track_h);

    free(secTracks_h);
    free(track_c);
  }

  //end of event-loop
  if(cudaEnabled) {
    cudaFree(physicsTable_d);
    cudaFree(sbData_d);
  }
  free(sbData);

  cudaFree(geom_d);
  free(geom_h);
  free(geom);

  cudaFree(bmap_d);
  free(bmap_h);
  free(fieldMap);
}
