/*
 * executable_test.C
 *
 */
// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

#include "TThread.h"
#include "Propagator.h"
#include "GeantScheduler.h"
#include "WorkloadManager.h"
#include "PhysicsProcess.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"
#include "TObject.h"
#include "UserApplication.h"
#include "ExN03Application.h"
#include "FactoryStore.h"
#include "Geant/Track.h"
#include "TROOT.h"
#include <iostream>
#include "GunGenerator.h"
#include "TTabPhysProcess.h"

void run(int nthreads=4,
         bool performance=true,
          const char *geomfile="ExN03.root",
          const char *xsec="xsec_FTFP_BERT.root",
          const char *fstate="fstate_FTFP_BERT.root",
         bool coprocessor = COPROCESSOR_REQUEST
        )
{
   // Those library used to need to be loaded explicitly and are now
   // automatically loaded by ROOT.
   // gSystem->Load("libPhysics");
   // gSystem->Load("libHist");
   // gSystem->Load("libThread");
   // gSystem->Load("libGeom");
   // gSystem->Load("libVMC");
   // gSystem->Load("../buildTGeo/lib/libGeant_v");
   // gSystem->Load("../buildTGeo/lib/libXsec");
   // gSystem->Load("../lib/libGeant_v");
   // gSystem->Load("../lib/libXsec");
   // gSystem->Load("../lib/libGeantExamples");
   // for vector physics - OFF now
   // gSystem->Load("../lib/libVphysproc");

//=============================================================================
// PERFORMANCE MODE SWITCH: no scoring, no memory cleanup thread, no monitoring
//=============================================================================
//   bool performance = true;

   int ntotal   = 50;  // Number of events to be transported
   int nbuffered  = 10;   // Number of buffered events (tunable [1,ntotal])
   TGeoManager::Import(geomfile);

   Propagator *prop = Propagator::Instance(ntotal, nbuffered, nthreads);
   // Monitor different features
   prop->SetNminThreshold(5*nthreads);
   if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
      CoprocessorBroker *gpuBroker = new CoprocessorBroker();
      gpuBroker->CudaSetup(32,128,1);
      prop->SetTaskBroker(gpuBroker);
#else
      std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
   }

   prop->SetMonitored(Propagator::kMonQueue,          true & (!performance));
   prop->SetMonitored(Propagator::kMonMemory,         false & (!performance));
   prop->SetMonitored(Propagator::kMonBasketsPerVol,  false & (!performance));
   prop->SetMonitored(Propagator::kMonVectors,        false & (!performance));
   prop->SetMonitored(Propagator::kMonConcurrency,    false & (!performance));
   prop->SetMonitored(Propagator::kMonTracksPerEvent, false & (!performance));
   bool graphics = (prop->GetMonFeatures()) ? true : false;
   prop->fUseMonitoring = graphics;
   prop->fNaverage = 500;   // Average number of tracks per event

   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
   prop->fPriorityThr = 0.05;

   // Initial vector size, this is no longer an important model parameter,
   // because is gets dynamically modified to accomodate the track flow
   prop->fNperBasket = 16;   // Initial vector size (tunable)

   // This is now the most important parameter for memory considerations
   prop->fMaxPerBasket = 256;   // Maximum vector size (tunable)

   prop->fEmin = 3.E-6; // [3 KeV] energy cut
   prop->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy

   // Create the tab. phys process.
   prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);

   // for vector physics -OFF now
   // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin, nthreads);

   prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
   prop->fLearnSteps = 0;
   if (performance) prop->fLearnSteps = 0;


   prop->fApplication = new ExN03Application();
   // Activate I/O
   prop->fFillTree = false;

// Activate debugging using -DBUG_HUNT=ON in your cmake build
   prop->fDebugEvt = 0;
   prop->fDebugTrk = 0;
   prop->fDebugStp = 0;
   prop->fDebugRep = 10;

// Activate standard scoring
   prop->fUseStdScoring = true;
   if (performance) prop->fUseStdScoring = false;
   // Monitor the application
   prop->fUseAppMonitoring = false;
   prop->PropagatorGeom(geomfile, nthreads, graphics);
   delete prop;
}

int main(int argc, char* argv[]) {
   if (argc < 6) {
      std::cerr << "Wrond usage - please add all possible values & paths to files";
      return 1;
   }
   run(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6]);
   return 0;
}
