#include "Rtypes.h"
#include "TGeoManager.h"

#include "Rtypes.h"
#include "TGeoManager.h"

#include "HepMC/HepMC.h"
#include "GunGenerator.h"
#include "base/messagelogger.h"
#include "HepMCGenerator.h"
#ifdef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#endif
#include "WorkloadManager.h"
#include "TTabPhysProcess.h"
#include "volumes/Particle.h"
#include "GeantPropagator.h"
#ifdef RESTORE_USE_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif
#include "ExN03Application.h"

// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

void loadvecgeomgeometry(GeantPropagator *prop);

int main() {
  Int_t nthreads = 15;
  const char *geomfile = "ExN03.root";
  const char *xsec = "xsec_FTFP_BERT.root";
  const char *fstate = "fstate_FTFP_BERT.root";
  bool coprocessor = false;

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
  Bool_t performance = true;
  double vt[3] = {-8, 0, 0};

  Int_t ntotal = 50;    // Number of events to be transported
  Int_t nbuffered = 10; // Number of buffered events (tunable [1,ntotal])
  TGeoManager::Import(geomfile);

  GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered);
  prop->fVertex[0] = vt[0];
  prop->fVertex[1] = vt[1];
  prop->fVertex[2] = vt[2];
  WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
  // Monitor different features
  wmgr->SetNminThreshold(5 * nthreads);
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32, 128, 1);
    wmgr->SetTaskBroker(gpuBroker);
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  wmgr->SetMonitored(WorkloadManager::kMonQueue, true & (!performance));
  wmgr->SetMonitored(WorkloadManager::kMonMemory, false & (!performance));
  wmgr->SetMonitored(WorkloadManager::kMonBasketsPerVol, false & (!performance));
  wmgr->SetMonitored(WorkloadManager::kMonVectors, false & (!performance));
  wmgr->SetMonitored(WorkloadManager::kMonConcurrency, false & (!performance));
  wmgr->SetMonitored(WorkloadManager::kMonTracksPerEvent, false & (!performance));
  Bool_t graphics = (wmgr->GetMonFeatures()) ? true : false;
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

#ifdef USE_VECGEOM_NAVIGATOR
  // Create particles
  Particle::CreateParticles();
#endif

  // Create the tab. phys process.
  loadvecgeomgeometry(prop);
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
