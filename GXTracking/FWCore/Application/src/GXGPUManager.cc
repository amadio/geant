#include "GXGPUManager.h"
#include "gxtracking_kernel.h"

GXGPUManager::GXGPUManager(int nb, int nt) : GXVCoprocessorManager("GPU"),
  fNBlocks(nb), fNThreads(nt)
{
  fNTracks = 0;
  fNSteps = 1;
  fRandomStates = 0;
  fPerformance = true;
  fTrackHandler = new GXTrackHandler();
}

GXGPUManager::~GXGPUManager()
{
  delete fTrackHandler;
  cudaStreamDestroy(stream);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

void GXGPUManager::Initialize()
{
  int nDevice;
  cudaGetDeviceCount(&nDevice);
  
  if(nDevice > 0) {
    cudaDeviceReset();
    std::cout << "... CUDA Kernel<<<" << fNBlocks << "," 
	      << fNThreads <<">>> (...) ..." << std::endl;
  }
  else {
    std::cout << "No Cuda Capable Device ... " << std::endl;
  }    
  
  //initialize the stream
  cudaStreamCreate(&stream);
  
  //create cuda events
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // maximum size of dynamic memory allocation on the device
  // set 1024 megabytes limit on the global mememory 
  fLimitHeapSize = 1024*1024*1024;
  cudaThreadSetLimit(cudaLimitMallocHeapSize, fLimitHeapSize);
}

void GXGPUManager::SetBlockThread(int nblocks, int nthreads)
{
  fNBlocks  = nblocks;
  fNThreads = nthreads;
} 

void GXGPUManager::SetNumberOfSteps(int nsteps) 
{
  fNSteps = nsteps;
}

void GXGPUManager::SetLimitHeapSize(size_t sizeBytes)
{
  fLimitHeapSize = sizeBytes;
}

void GXGPUManager::StartTimer() 
{
  cudaEventRecord (start,0);
}

float GXGPUManager::StopTimer() 
{
  float elapsedTime = 0.0;
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&elapsedTime,start,stop);
  return elapsedTime;
}

void GXGPUManager::AllocateDeviceMemory(/* GXTaskData* taskData_h */
					GXVGeometry *geom,
					GXFieldMap** fieldMap2D,
					GXPhysicsTable* physicsTable,
					GXPhysics2DVector* sbData)
{
  if(fPerformance) StartTimer();

  //Magnetic Field Map
  fieldMap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));
  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      fieldMap_h[i+j*nbinZ].Bz = fieldMap2D[i][j].Bz;
      fieldMap_h[i+j*nbinZ].Br = fieldMap2D[i][j].Br;
    }
  }

  // Geometry
  cudaMalloc( (void**)&geom_d, geom->size() );
  geom->relocate( geom_d );
  cudaMemcpy(geom_d, geom->getBuffer(), geom->size(), 
	     cudaMemcpyHostToDevice);

  cudaMalloc((void**)&fieldMap_d, nbinZ*nbinR*sizeof(GXFieldMap)) ;
  cudaMemcpy(fieldMap_d,fieldMap_h,nbinZ*nbinR*sizeof(GXFieldMap),
             cudaMemcpyHostToDevice);

  //Physics Tables
  cudaMalloc((void**)&physicsTable_d, 
             kNumberPhysicsTable*sizeof(GXPhysicsTable));
  cudaMemcpy(physicsTable_d, physicsTable, 
             kNumberPhysicsTable*sizeof(GXPhysicsTable),
             cudaMemcpyHostToDevice);

  //G4SeltzerBergerModel data
  cudaMalloc((void**)&sbData_d,maxElements*sizeof(GPPhysics2DVector));
  cudaMemcpy(sbData_d, sbData, maxElements*sizeof(GPPhysics2DVector),
             cudaMemcpyHostToDevice);

  //--->@@@G4FWP - host memeory allocation for CPU tests (will be removed)
  geom_h = (GPVGeometry::byte *) malloc (geom->size()) ;
  geom->relocate( geom_h );
  memcpy(geom_h,geom->getBuffer(),geom->size());  

 //physics tables
  physicsTable_h = 
    (GXPhysicsTable*) malloc (kNumberPhysicsTable*sizeof(GXPhysicsTable));
  memcpy(physicsTable_h,physicsTable,
	 kNumberPhysicsTable*sizeof(GXPhysicsTable));
 
  //DB data
  sbData_h = (GXPhysics2DVector*) malloc(maxElements*sizeof(GXPhysics2DVector));
  memcpy(sbData_h,sbData,sizeof(maxElements*sizeof(GXPhysics2DVector)));

  if(fPerformance) {
    float elapsedTime = StopTimer();
    printf("Time for AllocateDeviceMemory = %6.3f ms\n",elapsedTime);
  }
}

void GXGPUManager::DeallocateDeviceMemory()
{
  cudaFree(geom_d);
  cudaFree(fieldMap_d);
  cudaFree(physicsTable_d);
  cudaFree(sbData_d);

  free(geom_h);
  free(fieldMap_h);
  free(physicsTable_h);
  free(sbData_h);
}

void GXGPUManager::UploadTaskData()
{
  fNTracks = fTrackHandler->GetNumberOfTracks();

  if(fNTracks > 0) {
    GXTrack* tracks = fTrackHandler->GetTracks();

    //allocate and copy track information (@@@temporary)
    cudaHostAlloc((void**) &tracks_h, fNTracks*sizeof(GXTrack),
                  cudaHostAllocDefault);
    for(size_t i = 0 ; i < fNTracks ; i++){
       fTrackHandler->CopyTrack(&tracks[i],&tracks_h[i]);
    }

    //allocate and copy to GPU
    if(fPerformance) StartTimer();

    cudaMalloc((void**)&tracks_d, fNTracks*sizeof(GXTrack));
    cudaMemcpyAsync(tracks_d, tracks_h, fNTracks*sizeof(GXTrack), 
                    cudaMemcpyHostToDevice,stream);

    //@@@pre-allocated memory for secondaries 
    cudaMalloc((void**)&secTracks_d, 
	       fNSteps*maxSecondaryPerStep*fNTracks*sizeof(GXTrack));
    
    //atomic counter for the total number of secondaries
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));
    
    if(fPerformance) fElapsedTimeH2D = StopTimer();
  }
  else {
    printf("GXGPUManager::UploadTaskData: No tracks to upload to GPU\n");
  }
}

void GXGPUManager::DeallocateTaskData() {
  cudaFreeHost(tracks_h);
  cudaFree(tracks_d);

  cudaFree(secTracks_d);
  cudaFree(stackSize_d);
}

void GXGPUManager::LaunchTask() {

  if(fPerformance) StartTimer();

  //prepare random engines on the device
  cudaMalloc(&fRandomStates, fNBlocks*fNThreads* sizeof(curandState));
  curand_setup_gpu(fRandomStates, time(NULL), fNBlocks, fNThreads);

  elec_GPIL_gpu(fRandomStates, tracks_d, 0, fNTracks, 
		(GPGeomManager*)geom_d, fieldMap_d, physicsTable_d,sbData_d,
		secTracks_d, stackSize_d, fNSteps,
		fNBlocks, fNThreads, stream);
  
  //atomic counter for the last array position of physics processes
  cudaThreadSynchronize();

  cudaFree(fRandomStates);

  if(fPerformance) fElapsedTimeGPU = StopTimer();
}

void GXGPUManager::DownloadTaskData()
{
  if(fPerformance) StartTimer();

  //this is a dummy container
  GXTrack *tracks_o;
  cudaMalloc((void**)&tracks_o, fNTracks*sizeof(GXTrack));

  cudaMemcpyAsync(tracks_o, tracks_d, fNTracks*sizeof(GXTrack), 
		  cudaMemcpyDeviceToHost, stream);

  cudaFree(tracks_o);

  if(fPerformance) fElapsedTimeD2H = StopTimer();
}


void GXGPUManager::LaunchCPUTask() 
{
  //counter for the total number of secondaries
  G4int stackSize_h = 0;

  //@@@pre-allocated memory for secondaries
  GXTrack *secTracks_h
    = (GXTrack*) malloc(fNSteps*maxSecondaryPerStep*fNTracks*sizeof(GXTrack));

  if(fPerformance) StartTimer();

  elec_GPIL_cpu(tracks_h, 0, fNTracks, 
		(GPGeomManager*)geom_h, fieldMap_h, physicsTable_h, sbData_h,
		secTracks_h, &stackSize_h, fNSteps);
  
  if(fPerformance) fElapsedTimeCPU = StopTimer();

  free(secTracks_h);

}

void GXGPUManager::PrintPerformance(int taskId) 
{
  if(fPerformance) {
    printf("Task %d> Elapsed on H2D GPU D2H CPU : %6.3f %6.3f %6.3f %6.3f ms\n",
	   taskId,
	   fElapsedTimeH2D,fElapsedTimeGPU,fElapsedTimeD2H,fElapsedTimeCPU);
    if( fElapsedTimeCPU > 0.0) {
    printf("Task %d> Ratio GPU/CPU TotalGPU/CPU : %6.3f %6.3f\n",taskId,
	     fElapsedTimeGPU/fElapsedTimeCPU,
	     (fElapsedTimeH2D+fElapsedTimeGPU+fElapsedTimeD2H)/fElapsedTimeCPU);
    }
  }
  else {
    printf("Performance Flag = %d, Nothing to Report\n",fPerformance);
  }
}

