#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GXTrackLiason.h"
#include "GPEMPhysicsUtils.h"
#include "GPPhysicsTableType.h"
#include "trackingTest2_kernel.h"
#include "sort_kernel.h"
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
  int kernelType  = 0 ; //0 standard, 1 split-GPU, 2 split-GPU-CPU
  int runType  = 0; //0 normal 1 test
  int numStep  = 1; //number of steps

  if(argc >= 2) kernelType = atoi(argv[1]);
  if(argc >= 3) runType = atoi(argv[2]);
  if(argc >= 4) numStep = atoi(argv[3]);

  if(kernelType < 0 || kernelType > 15) {
    std::cout << "Usage: kernelTest mode=[0-4][11-14]: ... " << std::endl;
    return 0;
  }

  int nDevice;

  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
  }
  else {
    std::cout << "No Cuda Capable Device ... " << std::endl;
    return 0;
  }

  // set the default number of threads and thread blocks
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

  if(runType == 1 ) { 
    theNBlocks  =  4;
    theNThreads =  16;
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

  cudaMalloc((void**)&physicsTable_d, 
	     kNumberPhysicsTable*sizeof(GPPhysicsTable));
  cudaMemcpy(physicsTable_d, &physicsTable, 
	     kNumberPhysicsTable*sizeof(GPPhysicsTable),
	     cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&sbData_d,maxZ*sizeof(GPPhysics2DVector));
  cudaMemcpy(sbData_d, sbData, maxZ*sizeof(GPPhysics2DVector),
	     cudaMemcpyHostToDevice);

  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  //initialize the stream
  cudaStream_t stream0;
  cudaStream_t stream1;
  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);

  //event-loop
  int nevent = 10;
  int ndivarg = -1;
  G4double photonFraction = 0.2;

  for (int i = 0 ; i < nevent ; i++) {

    int nTracks = 4096*16;
    if(runType == 1 ) nTracks = theNBlocks*theNThreads;

    int nDiv = ndivarg > 0 ? ndivarg : 2;
    int Nchunk  =       0;
    int NchunkG  =      0;
    int NchunkE  =      0;

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    // populate tracks with track position, momentum, energy, steplength
    GXTrack *track_h;
    cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);

    //generate tracks with random position and momentum

    const G4double minP = 20.;
    const G4double maxP = 1000.;

    G4double rho, z, p, mass;
    G4double theta, phi;

    G4int nPhotons = 0;
    G4int nElectrons = 0;

    for(size_t i = 0 ; i < nTracks ; i++){

      rho = ecalRmim + (ecalRmax-ecalRmim)*Rndm();
      p = minP + (maxP - minP)*Rndm();
      z = ecalZmax*(2*Rndm()-1.0);
      phi = twopi*Rndm();
      theta = std::atan(rho/z);

      track_h[i].status = 0;
      track_h[i].proc = -1;
      track_h[i].q = ( Rndm() < photonFraction ) ? 0.0 : -1.0;
      track_h[i].x = rho*std::cos(phi);
      track_h[i].y = rho*std::sin(phi); 
      track_h[i].z = z; 
      track_h[i].s = 1.0; 

      track_h[i].px = p*std::sin(theta)*std::cos(phi);
      track_h[i].py = p*std::sin(theta)*std::sin(phi);
      track_h[i].pz = p*std::cos(theta);

      mass = electron_mass_c2*track_h[i].q*track_h[i].q;
      track_h[i].E  = p*p/(sqrt(p*p + mass*mass) + mass);

      if(track_h[i].q == 0.0) ++nPhotons;
      else ++nElectrons;
    }

    GXTrack *electron_h;
    GXTrack *photon_h;

    cudaHostAlloc((void**) &photon_h, nPhotons*sizeof(GXTrack),
		  cudaHostAllocDefault);
    GXTrack *photon_c = (GXTrack *) malloc (nPhotons*sizeof(GXTrack));

    cudaHostAlloc((void**) &electron_h, nElectrons*sizeof(GXTrack),
		  cudaHostAllocDefault);

    GXTrack *electron_c = (GXTrack *) malloc (nElectrons*sizeof(GXTrack));
    GXTrack *electron_sc = (GXTrack *) malloc (nElectrons*sizeof(GXTrack));

    G4int iphoton = 0;
    G4int ielectron = 0;
    
    for(size_t i = 0 ; i < nTracks ; i++){
      if(track_h[i].q == 0.0 && iphoton < nPhotons) {
	CopyTrack(&track_h[i],&photon_h[iphoton]);
	CopyTrack(&track_h[i],&photon_c[iphoton]);
	++iphoton;
      }
      if(track_h[i].q == -1.0 && ielectron < nElectrons) {
	CopyTrack(&track_h[i],&electron_h[ielectron]);     
	CopyTrack(&track_h[i],&electron_c[ielectron]);     
	++ielectron;
      }
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
    
    int nthreads_total = theNBlocks * theNThreads;

    GXTrack *photon_d;
    GXTrack *photon_d0;
    GXTrack *photon_d1;

    GXTrack *electron_d;
    GXTrack *electron_d0;
    GXTrack *electron_d1;

    NchunkG = (nPhotons/nDiv);
    if(NchunkG >= nthreads_total) {
      NchunkG -= ( (nPhotons/nDiv) % nthreads_total);
    } else {
      fprintf(stderr,"Warning: too many division for photons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkG,nthreads_total);
    }
    
    NchunkE = (nElectrons/nDiv);
    if (NchunkE >= nthreads_total) {
      NchunkE -= ( (nElectrons/nDiv) % nthreads_total); 
    } else {
      fprintf(stderr,"Warning: too many division for electrons (%d) resulting in low chunk size %d (< %d)\n",nDiv,NchunkE,nthreads_total);
    }
    
    cudaMalloc((void**)&photon_d, nPhotons*sizeof(GXTrack));
    cudaMalloc((void**)&photon_d0, NchunkG*sizeof(GXTrack));
    cudaMalloc((void**)&photon_d1, NchunkG*sizeof(GXTrack));
    
    cudaMalloc((void**)&electron_d, nElectrons*sizeof(GXTrack));
    cudaMalloc((void**)&electron_d0, NchunkE*sizeof(GXTrack));
    cudaMalloc((void**)&electron_d1, NchunkE*sizeof(GXTrack));
    
    //create associate GXTrackLiason
    GXTrackLiason *liason_d;
    cudaMalloc((void**)&liason_d, nElectrons*sizeof(GXTrackLiason));

    GXTrackLiason *liason_c 
      = (GXTrackLiason *) malloc (nElectrons*sizeof(GXTrackLiason));

    //stop time for cudaMalloc
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeAlloc,start,stop);
    
    //start time for Memcpy (H2D)
    cudaEventRecord (start,0);

    //stop time for Memcpy
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeUp,start,stop);

    //prepare random engines on the device
    curandState* devStates = 0;
    cudaMalloc(&devStates, theNBlocks*theNThreads* sizeof(curandState));
    curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);

    curandState* devStates1 = 0;
    cudaMalloc(&devStates1, theNBlocks*theNThreads* sizeof(curandState));
    curand_setup_gpu(devStates1, time(NULL), theNBlocks, theNThreads);

    //@@@pre-allocated memory for secondaries
    G4int stackSize = 0;
    GXTrack *secTracks_d;

    cudaMalloc((void**)&secTracks_d, 
	       numStep*maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //atomic counter for the total number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    GXTrack *secTracks_h
      = (GXTrack*) malloc(numStep*maxSecondaryPerStep*nTracks*sizeof(GXTrack));

    //<<<----------------------------------------------------------------->>>
    printf("Processing (nElectron,nPhoton) ");
    printf("= (%d,%d)\n", nElectrons,nPhotons);

    //start time for kernel
    cudaEventRecord (start,0);
    
    //electron kernel
    cudaMemcpyAsync(electron_d, electron_h, nElectrons*sizeof(GXTrack), 
		    cudaMemcpyHostToDevice, stream0);

    GXTrack *electron_t;
    cudaMalloc((void**)&electron_t, nElectrons*sizeof(GXTrack));

    GXTrack *electron_hs;
    cudaHostAlloc((void**) &electron_hs, nElectrons*sizeof(GXTrack),
		  cudaHostAllocDefault);

    G4int  nbrem_h = -1;
    G4int *nbrem_d;

    cudaMalloc((void**)&nbrem_d, sizeof(G4int));
    cudaMemset(nbrem_d,0,sizeof(G4int));
    
    G4int  nioni_h = -1;
    G4int *nioni_d;
    cudaMalloc((void**)&nioni_d, sizeof(G4int));
    cudaMemset(nioni_d,0,sizeof(G4int));
    
    G4int  stackSize_brem_h = -1;
    G4int *stackSize_brem;
    cudaMalloc((void**)&stackSize_brem, sizeof(G4int));
    cudaMemset(stackSize_brem,0,sizeof(G4int));
    
    G4int  stackSize_ioni_h = -1;
    G4int *stackSize_ioni;
    cudaMalloc((void**)&stackSize_ioni, sizeof(G4int));
    cudaMemset(stackSize_ioni,0,sizeof(G4int));
    
    if(kernelType == 0) {
      elec_gpu(devStates, electron_d, nElectrons, 
	       (GPGeomManager*)geom_d, bmap_d, physicsTable_d,sbData_d,
	       secTracks_d, stackSize_d, numStep, runType,
	       theNBlocks, theNThreads, stream0);
    }
    if(kernelType == 1 || kernelType == 2) {

      elec_GPIL_gpu(devStates, electron_d, liason_d, nElectrons, 
		    (GPGeomManager*)geom_d, bmap_d, physicsTable_d,sbData_d,
		    secTracks_d, stackSize_d, numStep, runType,
		    theNBlocks, theNThreads, stream0);

      //atomic counter for the last array position of physics processes
      //cudaThreadSynchronize();

      count_by_process_gpu(nElectrons, electron_d,
			   nbrem_d, nioni_d, theNBlocks, theNThreads, stream0);

      // cudaMemcpyAsync(&nbrem_h,nbrem_d,sizeof(G4int),cudaMemcpyDeviceToHost, stream0);
      // cudaMemcpyAsync(&nioni_h,nioni_d,sizeof(G4int),cudaMemcpyDeviceToHost, stream0);

      sort_by_process_gpu(nElectrons, electron_d, electron_t,
			  nbrem_d, stackSize_brem,
			  nioni_d, stackSize_ioni,
			  theNBlocks, theNThreads, stream0);

      // cudaThreadSynchronize();

      elec_doit_gpu(devStates, electron_t, liason_d, nElectrons, 
		    (GPGeomManager*)geom_d, bmap_d, physicsTable_d,sbData_d,
		    secTracks_d, stackSize_d, numStep, runType,
		    theNBlocks, theNThreads, stream0);
    }
    
    cudaMemcpyAsync(electron_h, electron_d, nElectrons*sizeof(GXTrack), 
		    cudaMemcpyDeviceToHost, stream0);
    
    //copy secondaries D2H 
    cudaThreadSynchronize();
    cudaMemcpyAsync(&stackSize,stackSize_d,sizeof(G4int),
		    cudaMemcpyDeviceToHost, stream0);
    
    //<<<----------------------------------------------------------------->>>

    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeKernel,start,stop);

    //start time for Memcpy (D2H)
    cudaEventRecord (start,0);
    
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

    nbrem_h = 0;
    nioni_h = 0;
    stackSize_brem_h = 0;
    stackSize_ioni_h = 0;
    
    cudaEventRecord (start,0);
    
    if(kernelType == 0 || kernelType == 1 ) {
      elec_cpu(electron_c, nElectrons, 
	       (GPGeomManager*)geom_h, bmap_h, physicsTable, sbData, 
	       secTracks_h, &stackSize, numStep, runType);
    }
    if(kernelType == 2 ) {
      elec_GPIL_cpu(electron_c, liason_c, nElectrons, 
		    (GPGeomManager*)geom_h, bmap_h, physicsTable, sbData, 
		    secTracks_h, &stackSize, numStep, runType);
      
      count_by_process_cpu(nElectrons, electron_c, &nbrem_h,&nioni_h); 
      
      sort_by_process_cpu(nElectrons, electron_c, electron_sc,
			  &nbrem_h, &stackSize_brem_h,
			  &nioni_h, &stackSize_ioni_h);
      
      elec_doit_cpu(electron_sc, liason_c, nElectrons, 
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

    cudaFree(devStates);
    cudaFree(devStates1);

    cudaFreeHost(track_h);
    cudaFreeHost(photon_h);
    cudaFreeHost(electron_h);
    free(photon_c);
    free(electron_c);
    free(electron_sc);

    cudaFree(photon_d);
    cudaFree(photon_d0);
    cudaFree(photon_d1);
    cudaFree(electron_d);
    cudaFree(electron_d0);
    cudaFree(electron_d1);

    cudaFree(electron_t);
    cudaFreeHost(electron_hs);
    
    cudaFree(secTracks_d);
    cudaFree(stackSize_d);
    free(secTracks_h);

    cudaFree(nbrem_d);
    cudaFree(nioni_d);
    cudaFree(stackSize_brem);
    cudaFree(stackSize_ioni);

    //deletec liason
    cudaFree(liason_d);
    free(liason_c);

  }

  //end of event-loop
  cudaStreamDestroy(stream0);
  cudaStreamDestroy(stream1);

  cudaFree(physicsTable_d);
  cudaFree(sbData_d);
  free(sbData);

  cudaFree(geom_d);
  free(geom_h);
  free(geom);

  cudaFree(bmap_d);
  free(bmap_h);
  free(fieldMap);
}
