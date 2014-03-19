#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "GXTrack.h"
#include "GPConstants.h"

//Magnetic Field
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include "stopwatch.h"

//Geometry
#include "GPVPhysicalVolume.h"

//FrameWork and Application
#include "GXVGeometry.h"
#include "GXUserGeometry.h"
#include "GXSimpleEcal.h"

#include <omp.h>

#include "GXRunManager.h"
#include "stepper_mic.h"

using namespace std;

// Variables may be decorated as a group using the push/pop method
   
#pragma offload_attribute(push, target(mic))
  GXTrack* track;
  GXFieldMap *bmap_h;
#pragma offload_attribute(pop)

// Or indvidually using either a __declspec or __attribute__

__attribute__ ((target(mic))) unsigned int ntracks = 32*4096;

int target_id = -1;


//#pragma offload_attribute(push,target(mic)) //{
//#include <iostream>
//#include <omp.h>
	 
__attribute__ ((target(mic))) void testThreadCount()
{
  int thread_count;
  
  #pragma omp parallel
  {
    #pragma omp single
    thread_count = omp_get_num_threads();
  }
  printf("--- Number of THREADS: %d\n",thread_count);
  //  cout << "Thread count: " << thread_count << endl;
}
	 
//#pragma offload_attribute(pop) //}


int main(int argc, char* argv[])
{
  //argument
  int stepperType = 1; //1 rk45, 2 rkf45, 3 nrk4
  if(argc >= 2) stepperType = atoi(argv[1]);
  if(stepperType < 1 || stepperType > 4) {
    std::cout << "Usage: stepper [1|2|3] ... " <<
      "1=RK4, 2=Felhberg, 3=Nystrom" << std::endl;
    return 0;
  }

  #pragma offload target(mic), if (argc > 1) 
  testThreadCount();

  // if offload compilation is enabled
  int num_devices = 0;
  int device_num = 0;

#ifdef __INTEL_OFFLOAD
  printf("--- Intel(R) Xeon Phi(TM) Devices ---\n");
  num_devices = _Offload_number_of_devices();
  printf("--- Number of Target Devices: %d\n",num_devices);

  device_num = _Offload_get_device_number();
  printf("--- Which Device number : %d\n",device_num);

#endif

  GXRunManager* runManager = GXRunManager::Instance();
  runManager->ConstructGeometry(new GXSimpleEcal());

  runManager->Initialization();

  //magnetic field on the device
  GXFieldMap** fieldMap = runManager->GetMagneticFieldMap(); 
  printf("Creating magnetic field map on the MIC device\n");
  
  bmap_h = (GXFieldMap *) _mm_malloc(nbinZ*nbinR*sizeof(GXFieldMap),128);

  for (int i = 0 ; i < nbinZ ; i++) {
    for (int j = 0 ; j < nbinR ; j++) {
      bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
      bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
    }
  }
  //free

  const int nthreads = 236;

  for(size_t ievt = 0 ; ievt < 10 ; ievt++){

    track = (GXTrack *) _mm_malloc(ntracks * sizeof(GXTrack),64);

    for(size_t i = 0 ; i < ntracks ; i++){
      track[i].x     = 300*(2.0*rand()/RAND_MAX-1.0);
      track[i].y     = 300*(2.0*rand()/RAND_MAX-1.0);
      track[i].z     = 300*(2.0*rand()/RAND_MAX-1.0);
      track[i].s     = 1.+100*(2.0*rand()/RAND_MAX-1.0);
      track[i].px    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track[i].py    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track[i].pz    = 1.+1000*(2.0*rand()/RAND_MAX-1.0);
      track[i].q     = -1.0;
    }

    StopWatch timer;
    timer.start();
    
    int k;
    if(stepperType == 1) {
    #pragma offload target(mic : target_id)		\
      in(bmap_h :length(nbinZ*nbinR))			\
      inout(track : length(ntracks)) 
      #pragma omp parallel for num_threads(nthreads)
      for (k = 0; k < ntracks ; ++k) {
        rk4_mic(bmap_h,&track[k],ntracks);
      }
    }
    else if(stepperType == 2) {
    #pragma offload target(mic : target_id) \
      in(bmap_h :length(nbinZ*nbinR))	    \
      inout(track : length(ntracks)) 
      #pragma omp parallel for num_threads(nthreads)
      for (k = 0; k < ntracks ; ++k) {
	rkf45_mic(bmap_h,&track[k],ntracks);
      }
    }
    else if(stepperType == 3) {
    #pragma offload target(mic : target_id) \
      in(bmap_h :length(nbinZ*nbinR)) \
      inout(track : length(ntracks)) 
      #pragma omp parallel for num_threads(nthreads)
      for (k = 0; k < ntracks ; ++k) {
	nrk4_mic(bmap_h,&track[k],ntracks);
      }
    }

    timer.stop();
    float elapsedTime_mic = timer.getTime();;

    timer.start();
    if(stepperType == 1)       rk4_cpu(bmap_h,track,ntracks);
    else if (stepperType == 2) rkf45_cpu(bmap_h,track,ntracks);
    else if (stepperType == 3) nrk4_cpu(bmap_h,track,ntracks);

    timer.stop();
    float elapsedTime_cpu = timer.getTime() - elapsedTime_mic;
    printf("Time Elapsed on MIC CPU : %6.3f %6.3f ms CPU/MIC = %6.3f\n",
	   elapsedTime_mic,elapsedTime_cpu,elapsedTime_cpu/elapsedTime_mic);

    // Deallocate
    _mm_free(track);

  }
  _mm_free(bmap_h);
}
