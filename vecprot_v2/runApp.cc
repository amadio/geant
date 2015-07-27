#include "Rtypes.h"
#include "TGeoManager.h"

#include "Rtypes.h"
#include "TGeoManager.h"

#include "HepMC/HepMC.h"
#include "GunGenerator.h"
#include "base/messagelogger.h"
#include "HepMCGenerator.h"
#include "TTabPhysProcess.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/Particle.h"
#endif
#include "ExN03Application.h"

// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

int main() {
  Int_t nthreads = 15;
  const char *geomfile = "ExN03.root";
  const char *xsec = "xsec_FTFP_BERT.root";
  const char *fstate = "fstate_FTFP_BERT.root";
  bool coprocessor = false;

  //=============================================================================
  // PERFORMANCE MODE SWITCH: no scoring, no memory cleanup thread, no monitoring
  //=============================================================================
  Bool_t performance = true;
  double vt[3] = {-8, 0, 0};

  Int_t ntotal = 50;    // Number of events to be transported
  Int_t nbuffered = 10; // Number of buffered events (tunable [1,ntotal])
  TGeoManager::Import(geomfile);

  GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered, nthreads);
  prop->fVertex[0] = vt[0];
  prop->fVertex[1] = vt[1];
  prop->fVertex[2] = vt[2];
  // Monitor different features
  prop->SetNminThreshold(5 * nthreads);
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32, 128, 1);
    prop->SetTaskBroker(gpuBroker);
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  prop->SetMonitored(GeantPropagator::kMonQueue, true & (!performance));
  prop->SetMonitored(GeantPropagator::kMonMemory, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonBasketsPerVol, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonVectors, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonConcurrency, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonTracksPerEvent, false & (!performance));
  Bool_t graphics = (prop->GetMonFeatures()) ? true : false;
  prop->fUseMonitoring = graphics;
  prop->fNaverage = 500; // Average number of tracks per event

  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  prop->fPriorityThr = 0.05;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  prop->fNperBasket = 16; // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  prop->fMaxPerBasket = 256; // Maximum vector size (tunable)

  prop->fEmin = 3.E-6; // [3 KeV] energy cut
  prop->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy

// Create the tab. phys process.
#ifdef USE_VECGEOM_NAVIGATOR
  prop->LoadVecGeomGeometry();
#endif
  prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
  prop->fProcess->Initialize();

  // for vector physics -OFF now
  // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin, nthreads);

  prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, vt[0], vt[1], vt[2], 1, 0, 0);

  // Number of steps for learning phase (tunable [0, 1e6])
  // if set to 0 disable learning phase
  prop->fLearnSteps = 0;
  if (performance)
    prop->fLearnSteps = 0;

  prop->fApplication = new ExN03Application();

  // Activate debugging using -DBUG_HUNT=ON in your cmake build
  prop->fDebugEvt = 0;
  prop->fDebugTrk = 0;
  prop->fDebugStp = 0;
  prop->fDebugRep = 10;

  // Activate standard scoring
  prop->fUseStdScoring = true;
  if (performance)
    prop->fUseStdScoring = false;
  // Monitor the application
  prop->fUseAppMonitoring = false;
  prop->PropagatorGeom(geomfile, nthreads, graphics);
  delete prop;
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::messagelogger::I()->summary(std::cout, "a");
#endif
}
