#include "RngBenchmarker.h"
#include "RngBenchmarker_gpu.h"

#include <cuda.h>
#include <curand_kernel.h>

#include "base/Global.h"
#include "base/Stopwatch.h"
using vecgeom::Stopwatch;

namespace vecphys {

void RngBenchmarker::RunCuda()
{
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    printf("Waning: No Cuda Capable Device ...\n");
  }

  //cuda event timing
  cudaEvent_t start;
  cudaEvent_t stop;

  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  //set the default number of threads and thread blocks - should be setable
  int theNBlocks  =   26;
  int theNThreads =  192;

  //1. MRG32k3a:

  MRG32k3a<ScalarBackend> *mrg32k2a = new MRG32k3a<ScalarBackend>();
  MRG32k3a_t<ScalarBackend>* statesMRG32k3a_d = 0; 
  cudaMalloc((void**)&statesMRG32k3a_d, theNBlocks*theNThreads*sizeof(MRG32k3a_t<ScalarBackend>));
  mrg32k2a->Initialize(statesMRG32k3a_d, theNBlocks, theNThreads);

  //2. Threefry:
  Threefry<ScalarBackend> *threefry = new Threefry<ScalarBackend>();
  Threefry_t<ScalarBackend>* statesThreefry_d = 0; 
  cudaMalloc((void**)&statesThreefry_d, theNBlocks*theNThreads*sizeof(Threefry_t<ScalarBackend>));
  threefry->Initialize(statesThreefry_d, theNBlocks, theNThreads);

  //Philox:
  Philox<ScalarBackend> *philox = new Philox<ScalarBackend>();
  Philox_t<ScalarBackend>* statesPhilox_d = 0; 
  cudaMalloc((void**)&statesPhilox_d, theNBlocks*theNThreads*sizeof(Philox_t<ScalarBackend>));
  philox->Initialize(statesPhilox_d, theNBlocks, theNThreads);

  //4 curandStateMRG32k3a
  curandStateMRG32k3a* devStatesMRG32k3a = 0;
  cudaMalloc(&devStatesMRG32k3a,theNBlocks*theNThreads*sizeof(curandStateMRG32k3a));
  curand_setup_gpu(devStatesMRG32k3a, time(NULL), theNBlocks, theNThreads);

  //4 curandStatePhilox
  curandStatePhilox4_32_10_t* devStatesPhilox = 0;
  cudaMalloc(&devStatesPhilox,theNBlocks*theNThreads*sizeof(curandStateMRG32k3a));
  curand_setup_gpu(devStatesPhilox, time(NULL), theNBlocks, theNThreads);

  //return values for varification
  double *result_h;
  double *result_c;
  double *result_d;

  result_h = (double*) calloc(theNBlocks*theNThreads, sizeof(double));
  result_c = (double*) calloc(theNBlocks*theNThreads, sizeof(double));
  cudaMalloc((void**)&result_d,theNBlocks*theNThreads*sizeof(double));

  float elapsedEventTime[kNumberRng +2];
  float elapsedTotalTime[kNumberRng +2];
  double rngTotal[kNumberRng +2];
  double rngEvent[kNumberRng +2];

  for (int k = 0; k < kNumberRng + 2 ; ++k) {
    elapsedTotalTime[k] = 0.;
    rngTotal[k] = 0.;
  }

  for (unsigned r = 0; r < fRepetition; ++r) {

    for (int k = 0; k < kNumberRng + 2; ++k) {
      elapsedEventTime[k] = 0.0;
      rngEvent[k] = 0.0;

      cudaMemset(result_d,0,theNBlocks*theNThreads*sizeof(double));
  
      if(cudaEnabled) {
        cudaEventRecord (start,0);

        //call CUDA kernel
	if(k == 0) {
	  CudaMRG32k3a(statesMRG32k3a_d, result_d, fNSample, theNBlocks, theNThreads);
        }
	if(k == 1) {
     	  CudaThreefry(statesThreefry_d, result_d, fNSample, theNBlocks, theNThreads);
	}
	if(k == 2) {
	  CudaPhilox(statesPhilox_d, result_d, fNSample, theNBlocks, theNThreads);
	}

	if(k == 3) {
          CurandMRG32k3a(devStatesMRG32k3a,result_d,fNSample,theNBlocks,theNThreads);
	}

	if(k == 4) {
          CurandPhilox(devStatesPhilox,result_d,fNSample,theNBlocks,theNThreads);
	}

        cudaEventRecord (stop,0);
        cudaEventSynchronize (stop);
        cudaEventElapsedTime (&elapsedEventTime[k],start,stop);

        //copy the result for varification
        cudaMemcpy(result_h,result_d,theNBlocks*theNThreads*sizeof(double),
                cudaMemcpyDeviceToHost);
         
        for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) rngEvent[k] += result_h[k]; 
      }
      elapsedTotalTime[k] += elapsedEventTime[k]/1000.; //ms to second
      rngTotal[k] += rngEvent[k]; //ms to second
    }
  }

  for (int k = 0; k < kNumberRng + 2; ++k) {
    if(k < kNumberRng) {
      printf(" %s  CudaBackend   Total time = %7.4f sec CheckSum = %g\n", 
   	     RngName[k], elapsedTotalTime[k], rngTotal[k]);
    }
    if(k== kNumberRng) {
      printf("%s Nvidia   Total time = %7.4f sec CheckSum = %g\n", 
   	     "CurandMRG32k3a", elapsedTotalTime[k], rngTotal[k]);
    }
    if(k== kNumberRng+1) {
      printf("%s Nvidia   Total time = %7.4f sec CheckSum = %g\n", 
   	     "CurandPhilox  ", elapsedTotalTime[k], rngTotal[k]);
    }
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(statesMRG32k3a_d);
  cudaFree(statesThreefry_d);
  cudaFree(statesPhilox_d);
  cudaFree(devStatesMRG32k3a);
  cudaFree(devStatesPhilox);

  cudaFree(result_d);
  free(result_h);
  free(result_c);

  delete mrg32k2a;
  delete threefry;
  delete philox;
}

} // end of vecphys namespace

