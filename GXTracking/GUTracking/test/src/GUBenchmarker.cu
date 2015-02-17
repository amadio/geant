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
			 int nTrackSize, GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  vecphys::cuda::GUComptonKleinNishina model(devStates,tid);
  //construct model here

  while (tid < nTrackSize) {
    model.Interact(itrack[tid],targetElements[tid],&otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

} // end namespace cuda

void GUBenchmarker::RunCuda()
{
  GUTrackHandler *handler_in = new GUTrackHandler();
  handler_in->GenerateRandomTracks(fNtracks);

  GUTrack* track_aos = handler_in->GetAoSTracks();
  GUTrack* track_aos_out = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }
  
  int *targetElements_d;
  cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
  cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

  GUTrack *itrack_d;
  GUTrack *otrack_d;
  cudaMalloc((void**)&itrack_d, fNtracks*sizeof(GUTrack));
  cudaMalloc((void**)&otrack_d, fNtracks*sizeof(GUTrack));

  cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GUTrack), 
               cudaMemcpyHostToDevice);

  //set the default number of threads and thread blocks - should be setable
  int theNBlocks  =  26;
  int theNThreads = 192;

  //prepare random engines on the device
  Random_t* randomStates = 0;

  cudaMalloc(&randomStates, theNBlocks*theNThreads* sizeof(curandState));
  GUCurand_Init(randomStates, time(NULL), theNBlocks, theNThreads);

  Stopwatch timer;
  timer.Start();

  for (unsigned r = 0; r < fRepetitions; ++r) {
    vecphys::cuda::BenchmarkCudaKernel<<<theNBlocks, theNThreads>>>(
				       randomStates,fNtracks,
				       itrack_d, targetElements_d,otrack_d);
  }

  Precision elapsedCuda = timer.Stop();

  cudaMemcpy(otrack_d, track_aos_out, fNtracks*sizeof(GUTrack), 
               cudaMemcpyDeviceToHost);

  if (fVerbosity > 0) {
    printf("CUDA  : %6.3fs\n",elapsedCuda);
  }
  cudaDeviceSynchronize();
}

} // end of vecphys namespace

