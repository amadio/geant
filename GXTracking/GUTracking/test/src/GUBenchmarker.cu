#include "GUBenchmarker.h"

#include "base/Stopwatch.h"
#include "backend/cuda/Backend.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

#include "GUCurand.h"

namespace vecphys {
inline namespace cuda {

__global__
void BenchmarkCudaKernel(Random_t* devStates,
                         GUAliasTable* table,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,10,1.,1001.,100,100,table);
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
  cudaThreadSetLimit(cudaLimitMallocHeapSize, 1024*1024*1024);

  //prepare table
  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUAliasTable* table_h =  model->GetSampler()->GetAliasTable();

  GUAliasTable* table_d;
  cudaMalloc((void**)&table_d,table_h->SizeOfTable());
  table_h->Relocate(table_d);

  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }
  
  int *targetElements_d;
  cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
  cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

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

  for (unsigned r = 0; r < fRepetitions; ++r) {

    fTrackHandler->GenerateRandomTracks(fNtracks);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();
    cudaMemcpy(itrack_d, itrack_aos, fNtracks*sizeof(GUTrack), 
                 cudaMemcpyHostToDevice);

    timer.Start();
    vecphys::cuda::BenchmarkCudaKernel<<<theNBlocks, theNThreads>>>(
				       randomStates,table_d,fNtracks,
				       itrack_d, targetElements_d,otrack_d);
    cudaDeviceSynchronize();
    Precision elapsedCuda = timer.Stop();

    cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GUTrack), 
               cudaMemcpyDeviceToHost);

    if (fVerbosity > 0) {
      printf("CUDA   Task %d >: %6.3fs\n",r,elapsedCuda);
    }

  }

  if (fVerbosity > 1) {
    for(unsigned i = 0; i < 4 ; ++i) printf(" E[%d]= %f\n",i,otrack_aos[i].E);
  }

  cudaFree(randomStates);
  cudaFree(itrack_d);
  cudaFree(otrack_d);
  cudaFree(table_d);
  cudaFree(targetElements_d);

  free(table_h);
  free(targetElements);
  free(otrack_aos);
  //  delete model;
}

} // end of vecphys namespace

