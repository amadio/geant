#include "GUBenchmarker.h"
#include "GUBenchmarker_gpu.h"
#include "GUPhysicsModelName.h"

#include "base/Global.h"
#include "base/Stopwatch.h"
using vecgeom::Stopwatch;

#include "GUAliasSampler.h"
//#include "GVComptonKleinNishina.h"
#include "ComptonKleinNishina.h"
#include "ConversionBetheHeitler.h"
#include "PhotoElectronSauterGavrila.h"
#include "IonisationMoller.h"
#include "BremSeltzerBerger.h"
#include "Physics2DVector.h"
#include "GUTrackHandler.h"

#include "SamplingMethod.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

#include "GUCurand.h"

namespace vecphys {

void GUBenchmarker::RunCuda()
{
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
    //    printf("CUDA Enabled with %d Device(s)\n",nDevice);
  }
  else {
    printf("Waning: No Cuda Capable Device ...");
  }

  //set 1024 megabytes on the heap (global mememory) 
  //cudaThreadSetLimit(cudaLimitMallocHeapSize, 1024*1024*1024);

  //cuda event timing
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

#ifdef VECPHYS_ROOT
  GUHistogram* histogram = new GUHistogram("cuda.root", fMaxP);
#endif

  int *targetElements = new int [fNtracks];
  int *targetElements_d;

  //prepare table - this step may be move to the physics list later
  GUAliasTableManager** tableM_h = 
    (GUAliasTableManager**) malloc(kNumberPhysicsModel*sizeof(GUAliasTableManager*)); 

  ComptonKleinNishina *KleinNishina;
  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kKleinNishina ) {      
    KleinNishina = new ComptonKleinNishina(0,-1);
    KleinNishina->Initialization();
    tableM_h[kKleinNishina] = KleinNishina->GetSampler()->GetAliasTableManager();
  }

  ComptonKleinNishina *HybridKleinNishina;
  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kHybridKleinNishina) {      
    HybridKleinNishina = new ComptonKleinNishina(0,-1);
    HybridKleinNishina->Initialization();
    tableM_h[kHybridKleinNishina] = HybridKleinNishina->GetSampler()->GetAliasTableManager();
  }

  ConversionBetheHeitler *BetheHeitler;
  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kBetheHeitler) {      
    BetheHeitler = new ConversionBetheHeitler(0,-1);
    BetheHeitler->Initialization();
    tableM_h[kBetheHeitler] = BetheHeitler->GetSampler()->GetAliasTableManager();
  }

  PhotoElectronSauterGavrila *SauterGavrila;
  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kSauterGavrila) {      
    SauterGavrila = new PhotoElectronSauterGavrila(0,-1);
    SauterGavrila->Initialization();
    tableM_h[kSauterGavrila] = SauterGavrila->GetSampler()->GetAliasTableManager();
  }

  IonisationMoller *MollerBhabha;
  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kMollerBhabha) {      
    MollerBhabha = new IonisationMoller(0,-1);
    MollerBhabha->Initialization();
    tableM_h[kMollerBhabha] = MollerBhabha->GetSampler()->GetAliasTableManager();
  }

  BremSeltzerBerger *SeltzerBerger;
  Physics2DVector* sbData;
  Physics2DVector* sbData_d;

  if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == kSeltzerBerger) {      
    SeltzerBerger = new BremSeltzerBerger(0,-1);
    tableM_h[kSeltzerBerger] = SeltzerBerger->GetSampler()->GetAliasTableManager();
  //SeltzerBerger data
    sbData = SeltzerBerger->GetSBData();
    cudaMalloc((void**)&sbData_d,maximumZ*sizeof(Physics2DVector));
    cudaMemcpy(sbData_d, sbData, maximumZ*sizeof(Physics2DVector),
	       cudaMemcpyHostToDevice);
  }

  GUAliasTableManager** tableM_d;
  cudaMalloc((void**)&tableM_d,kNumberPhysicsModel*sizeof(GUAliasTableManager*));

  GUAliasTableManager* temp_d[kNumberPhysicsModel];
  for(int i = 0 ; i < kNumberPhysicsModel ; ++i) {
    if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == i) {      
      cudaMalloc((void**)&temp_d[i],tableM_h[i]->SizeOfManager());
      tableM_h[i]->Relocate(temp_d[i]);
    }
  }

  cudaMemcpy(tableM_d,temp_d,sizeof(GUAliasTableManager*)*kNumberPhysicsModel,
             cudaMemcpyHostToDevice);

  GUTrack* itrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
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

  float elapsedTotal[kNumberPhysicsModel];
  float elapsedEventTime[kNumberPhysicsModel];

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) elapsedTotal[k] = 0.; 

  for (unsigned r = 0; r < fRepetitions; ++r) {

    //    PrepareTargetElements(targetElements, fNtracks);
    fMaterialHandler->PrepareTargetElements(targetElements, fNtracks, fMaterialMode);

    //H2D
    cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
    cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

    fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);

    GUTrack* track_aos = fTrackHandler->GetAoSTracks();
    fTrackHandler->SortAoSTracksByEnergy(track_aos,fNtracks);

    for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
      if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == k) {     

	fTrackHandler->CopyAoSTracks(track_aos,itrack_aos,fNtracks);
	cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GUTrack), 
                   cudaMemcpyHostToDevice);
	
	elapsedEventTime[k] = 0.0;
	
	if(cudaEnabled) {
	  cudaEventRecord (start,0);
	  //call CUDA kernels
	  CudaKernelFunc[k](theNBlocks, theNThreads, randomStates,
			    tableM_d,sbData_d,fNtracks, itrack_d, 
			    targetElements_d,otrack_d,fSampleType);
	  cudaEventRecord (stop,0);
	  cudaEventSynchronize (stop);
	  cudaEventElapsedTime (&elapsedEventTime[k],start,stop);
	}
	elapsedTotal[k] += elapsedEventTime[k]/1000.; //ms to second
	
	cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GUTrack), 
		   cudaMemcpyDeviceToHost);
	cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GUTrack), 
		   cudaMemcpyDeviceToHost);
	
#ifdef VECPHYS_ROOT
	histogram->RecordTime(k,elapsedEventTime[k]);
	for(int i = 0 ; i < fNtracks ; ++i) {
	  histogram->RecordHistos(k,track_aos[i].E,
				  itrack_aos[i].E,
				  itrack_aos[i].pz/itrack_aos[i].E,
				  otrack_aos[i].E,
				  otrack_aos[i].pz/otrack_aos[i].E);
	} 
#endif    
      }
    }
  }

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
    if(fEmModel == GUPhysicsModelIndex::kNullModel || fEmModel == k) {     
      printf("%s  CUDA   Total time of %3d reps = %7.4f sec\n",
	     GUPhysicsModelName[k], fRepetitions, elapsedTotal[k]);
    }
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(randomStates);
  cudaFree(itrack_d);
  cudaFree(otrack_d);
  cudaFree(tableM_d);
  cudaFree(targetElements_d);

  free(tableM_h);
  free(targetElements);
  free(itrack_aos);
  free(otrack_aos);

  //  delete model;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

} // end of vecphys namespace

