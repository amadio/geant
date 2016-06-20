#include "ProcessBenchmarker.h"
#include "ProcessBenchmarker_gpu.h"

#include "base/Global.h"
#include "base/Stopwatch.h"

#include "PhotonProcess.h"
#include "ElectronProcess.h"

#include "GUTrackHandler.h"
#include "GUTrack.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

#include "GUCurand.h"

namespace vecphys {

void ProcessBenchmarker::RunCuda()
{
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    printf("Waning: No Cuda Capable Device ...");
  }

  //cuda event timing
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

#ifdef VECPHYS_ROOT
  GUHistogram* histogram = new GUHistogram("process_cuda.root", fMaxP);
#endif

  int *targetElements = new int [fNtracks];
  int *targetElements_d;

  //prepare table - this step may be move to the physics list later

  fMaterialHandler->BuildMaterialTable();

  PhotonProcess *gammaProcess = new PhotonProcess(0,-1);
  ElectronProcess *electronProcess = new ElectronProcess(0,-1);

  int nbins = gammaProcess->GetNumberOfEnergyBin() * gammaProcess->GetNumberOfMaterialBin();

  CrossSectionData** data_h = (CrossSectionData**) malloc(sizeof(CrossSectionData*)*kNumberPhysicsProcess); 
  data_h[0] =  gammaProcess->GetCrossSectionData();
  data_h[1] =  electronProcess->GetCrossSectionData();

  CrossSectionData* data_d[kNumberPhysicsProcess];
  cudaMalloc((void**)&data_d,kNumberPhysicsProcess*sizeof(CrossSectionData*));

  for(int i = 0 ; i < kNumberPhysicsProcess ; ++i) {
    cudaMalloc((void**)&data_d[i],sizeof(CrossSectionData)*nbins);
    cudaMemcpy(data_d[i],data_h[i],sizeof(CrossSectionData)*nbins,cudaMemcpyHostToDevice);
  }


  GUTrack* itrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  //allocate memory for input/output tracks
  GUTrack *itrack_d;
  cudaMalloc((void**)&itrack_d, fNtracks*sizeof(GUTrack));

  //set the default number of threads and thread blocks - should be setable
  int theNBlocks  =  26;
  int theNThreads = 192;

  //prepare random engines on the device
   Random_t* randomStates = 0;
   cudaMalloc(&randomStates, theNBlocks*theNThreads* sizeof(curandState));
   GUCurand_Init(randomStates, time(NULL), theNBlocks, theNThreads);

  float elapsedTotal[kNumberPhysicsProcess];
  float elapsedEventTime[kNumberPhysicsProcess];

  for (int k = 0; k < kNumberPhysicsProcess; ++k)
    elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {

    fMaterialHandler->PrepareMaterialIndex(targetElements, fNtracks, fMaterialMode);

    //Host2Device
    cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
    cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

    fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);

    GUTrack* track_aos = fTrackHandler->GetAoSTracks();
    fTrackHandler->SortAoSTracksByEnergy(track_aos,fNtracks);

    for (int k = 0; k < kNumberPhysicsProcess; ++k) {
      if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {

        fTrackHandler->CopyAoSTracks(track_aos,itrack_aos,fNtracks);
        cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GUTrack), cudaMemcpyHostToDevice);
	
        elapsedEventTime[k] = 0.0;
	
        if(cudaEnabled) {
          cudaEventRecord (start,0);

          //call CUDA kernels
          CudaKernelFunc[k](theNBlocks, theNThreads, randomStates,
                            data_d[k], fNtracks, itrack_d, targetElements_d);

          cudaEventRecord (stop,0);
          cudaEventSynchronize (stop);
          cudaEventElapsedTime (&elapsedEventTime[k],start,stop);
        }
        elapsedTotal[k] += elapsedEventTime[k]/1000.; //ms to second
       
        cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GUTrack), cudaMemcpyDeviceToHost);
	
#ifdef VECPHYS_ROOT
        histogram->RecordTime(k,elapsedEventTime[k]);
        for(int i = 0 ; i < fNtracks ; ++i) {
          histogram->RecordHistosProc(k,
                                      itrack_aos[i].E,
                                      itrack_aos[i].nint,
                                      itrack_aos[i].s,
                                      itrack_aos[i].lambda);
        } 
#endif    
      }
    }
  }
  for (int k = 0; k < kNumberPhysicsProcess; ++k) {
    if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
      printf("%s  Cuda Total time of %3d reps = %7.4f sec\n", GUPhysicsProcessName[k], fRepetitions, elapsedTotal[k]);
    }
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  free(data_h);
  free(itrack_aos);
  free(targetElements);

  cudaFree(data_d);
  cudaFree(itrack_d);
  cudaFree(targetElements_d);
  cudaFree(randomStates);

#ifdef VECPHYS_ROOT
  delete histogram;
#endif

}



} // end of vecphys namespace

