#include "GUBenchmarker.h"
#include "GUBenchmarker_gpu.h"
#include "GUPhysicsModelName.h"

#include "base/Stopwatch.h"
#include "backend/cuda/Backend.h"

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
  cudaDeviceReset();

  //set 1024 megabytes on the heap (global mememory) 
  //cudaThreadSetLimit(cudaLimitMallocHeapSize, 1024*1024*1024);

#ifdef VECPHYS_ROOT
  GUHistogram* histogram = new GUHistogram("cuda.root", fMaxP);
#endif

  int *targetElements = new int [fNtracks];
  int *targetElements_d;

  //prepare table - this step may be move to the physics list later
  ComptonKleinNishina *KleinNishina = new ComptonKleinNishina(0,-1);
  ConversionBetheHeitler *BetheHeitler = new ConversionBetheHeitler(0,-1);
  PhotoElectronSauterGavrila *SauterGavrila = new PhotoElectronSauterGavrila(0,-1);
  IonisationMoller *MollerBhabha = new IonisationMoller(0,-1);
  BremSeltzerBerger *SeltzerBerger = new BremSeltzerBerger(0,-1);

  GUAliasTableManager** tableM_h = 
    (GUAliasTableManager**) malloc(kNumberPhysicsModel*sizeof(GUAliasTableManager*)); 

  tableM_h[kKleinNishina]  = KleinNishina->GetSampler()->GetAliasTableManager();
  tableM_h[kBetheHeitler]  = BetheHeitler->GetSampler()->GetAliasTableManager();
  tableM_h[kSauterGavrila] = SauterGavrila->GetSampler()->GetAliasTableManager();
  tableM_h[kMollerBhabha]  = MollerBhabha->GetSampler()->GetAliasTableManager();
  tableM_h[kSeltzerBerger] = SeltzerBerger->GetSampler()->GetAliasTableManager();

  GUAliasTableManager** tableM_d;
  cudaMalloc((void**)&tableM_d,kNumberPhysicsModel*sizeof(GUAliasTableManager*));

  //  KleinNishina->GetSampler()->PrintTable();
  
  GUAliasTableManager* temp_d[kNumberPhysicsModel];
  for(int i = 0 ; i < kNumberPhysicsModel ; ++i) {
    cudaMalloc((void**)&temp_d[i],tableM_h[i]->SizeOfManager());
    tableM_h[i]->Relocate(temp_d[i]);
  }

  cudaMemcpy(tableM_d,temp_d,sizeof(GUAliasTableManager*)*kNumberPhysicsModel,
             cudaMemcpyHostToDevice);

  //SeltzerBerger data
  Physics2DVector* sbData = SeltzerBerger->GetSBData();
  Physics2DVector* sbData_d;
  cudaMalloc((void**)&sbData_d,maximumZ*sizeof(Physics2DVector));
  cudaMemcpy(sbData_d, sbData, maximumZ*sizeof(Physics2DVector),
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

  Precision elapsedTotal[kNumberPhysicsModel];
  Precision elapsedT[kNumberPhysicsModel];

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) elapsedTotal[k] = 0.; 

  for (unsigned r = 0; r < fRepetitions; ++r) {

    PrepareTargetElements(targetElements, fNtracks);

    //H2D
    cudaMalloc((void**)&targetElements_d, fNtracks*sizeof(int));
    cudaMemcpy(targetElements_d, targetElements, fNtracks*sizeof(int), 
               cudaMemcpyHostToDevice);

    fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);
    GUTrack* track_aos = fTrackHandler->GetAoSTracks();

    for(unsigned int k = 0 ; k < kNumberPhysicsModel ; ++k) {

      fTrackHandler->CopyAoSTracks(track_aos,itrack_aos);
      cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GUTrack), 
                   cudaMemcpyHostToDevice);

      elapsedT[k] = 0.0;
      elapsedT[k] = CudaKernelFunc[k](theNBlocks, theNThreads, randomStates,
                                      tableM_d,sbData_d,fNtracks, itrack_d, 
                                      targetElements_d,otrack_d,fSampleType);
      elapsedTotal[k] += elapsedT[k];

      cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GUTrack), 
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GUTrack), 
                 cudaMemcpyDeviceToHost);

#ifdef VECPHYS_ROOT
      histogram->RecordTime(k,elapsedT[k]);
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

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
    printf("%s  CUDA   Total time of %3d reps = %6.3f sec\n",
           GUPhysicsModelName[k], fRepetitions, elapsedTotal[k]);
  }

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

