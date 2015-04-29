#include "GUBenchmarker.h"

#include "base/Stopwatch.h"
#include "backend/cuda/Backend.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

#include "GUCurand.h"

namespace vecphys {

inline namespace cuda {

__global__
void BenchmarkCudaKernel(Random_t* devStates,
			 //                         GUAliasTable* table,
                         GUAliasTableManager* table,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.,1001.,100,100,table);

  GUComptonKleinNishina model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

} // end namespace cuda

void GUBenchmarker::RunCuda()
{
  cudaDeviceReset();

  //set 1024 megabytes on the heap (global mememory) 
  //  cudaThreadSetLimit(cudaLimitMallocHeapSize, 1024*1024*1024);

#ifdef VECPHYS_ROOT
  GUHistogram* histogram = new GUHistogram("cuda.root", fMaxP);
#endif

  int *targetElements = new int [fNtracks];
  int *targetElements_d;

  //prepare table
  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  GUAliasTableManager* tableM_h = model->GetSampler()->GetAliasTableManager();
  GUAliasTableManager* tableM_d;
  cudaMalloc((void**)&tableM_d, tableM_h->SizeOfManager());

  tableM_h->Relocate(tableM_d);
  
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  //allocate memory for input/output tracks
  GUTrack *itrack_d;
  GUTrack *otrack_d;
  cudaMalloc((void**)&itrack_d, fNtracks*sizeof(GUTrack));
  cudaMalloc((void**)&otrack_d, fNtracks*sizeof(GUTrack));

  //set the default number of threads and thread blocks - should be setable
  int theNBlocks  =  26;
  int theNThreads = 192;

  //prepare random engines on the device
  Random_t* randomStates = 0;
  cudaMalloc(&randomStates, theNBlocks*theNThreads* sizeof(curandState));
  GUCurand_Init(randomStates, time(NULL), theNBlocks, theNThreads);

  Stopwatch timer;
  Precision elapsedCudaTotal= 0.0;

  Precision* incomingEn = new Precision[fNtracks];

  for (unsigned r = 0; r < fRepetitions; ++r) {

    PrepareTargetElements(targetElements, fNtracks);

    //H2D
    cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
    cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

    fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();
    cudaMemcpy(itrack_d, itrack_aos, fNtracks*sizeof(GUTrack), 
                 cudaMemcpyHostToDevice);

    for(int i = 0 ; i < fNtracks ; ++i) {
      incomingEn[i] = itrack_aos[i].E;
    }

    timer.Start();
    vecphys::cuda::BenchmarkCudaKernel<<<theNBlocks, theNThreads>>>(
				       randomStates,tableM_d,fNtracks,
    				       itrack_d, targetElements_d,otrack_d);
    cudaDeviceSynchronize();
    Precision elapsedCuda = timer.Stop();
    elapsedCudaTotal += elapsedCuda;

    if (fVerbosity > 0) {
      printf("CUDA   Task %d > %6.3f sec\n",r,elapsedCuda);
    }

    cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GUTrack), 
               cudaMemcpyDeviceToHost);
    cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GUTrack), 
               cudaMemcpyDeviceToHost);

#ifdef VECPHYS_ROOT
    histogram->RecordTime(elapsedCuda);

    for(int i = 0 ; i < fNtracks ; ++i) {
       histogram->RecordHistos( incomingEn[i],
                                itrack_aos[i].E,   
                                itrack_aos[i].pz/itrack_aos[i].E, 
                                otrack_aos[i].E,                  
                                otrack_aos[i].pz/otrack_aos[i].E);
    }
#endif    

  }

  printf("Cuda   Task Total time of %3d reps = %6.3f sec\n",fRepetitions ,elapsedCudaTotal);

  cudaFree(randomStates);
  cudaFree(itrack_d);
  cudaFree(otrack_d);
  cudaFree(tableM_d);
  cudaFree(targetElements_d);

  free(tableM_h);
  free(targetElements);
  free(otrack_aos);

  //  delete model;
  delete[] incomingEn; 
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

} // end of vecphys namespace

