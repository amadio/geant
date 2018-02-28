// The following in ROOT v6 equivalent to gSystem->Load("../lib/libGeant_v");
// R__LOAD_LIBRARY(libGeant_v)

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

// Autoload the library early so that Propagator is defined when applicable.
namespace geant {
inline namespace cxx {
class TaskBroker;
class Propagator;
}
}

void run(int ncputhreads = 1, bool performance = true, const char *geomfile = "ExN03.root",
         const char *xsec = "xsec_FTFP_BERT.root", const char *fstate = "fstate_FTFP_BERT.root",
         bool coprocessor = COPROCESSOR_REQUEST)
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

  int nthreads     = ncputhreads;
  int ntotal       = 50; // Number of events to be transported
  int nbuffered    = 10; // Number of buffered events (tunable [1,ntotal])
  int npropagators = 1;

  geant::TaskBroker *broker = nullptr;
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32, 128, 1);
    broker = gpuBroker;

    nthreads += gpuBroker->GetNstream() + 1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  geant::GeantConfig *config = new geant::GeantConfig();
  config->fGeomFileName      = geomfile;

  // Monitor different features
  config->fNminThreshold = 5 * nthreads;
  config->SetMonitored(GeantConfig::kMonQueue, true & (!performance));
  config->SetMonitored(GeantConfig::kMonMemory, false & (!performance));
  config->SetMonitored(GeantConfig::kMonBasketsPerVol, false & (!performance));
  config->SetMonitored(GeantConfig::kMonVectors, false & (!performance));
  config->SetMonitored(GeantConfig::kMonConcurrency, false & (!performance));
  config->SetMonitored(GeantConfig::kMonTracksPerEvent, false & (!performance));
  bool graphics          = (config->GetMonFeatures()) ? true : false;
  config->fUseMonitoring = graphics;
  config->fNaverage      = 500; // Average number of tracks per event

  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = 0.05;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = 16; // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = 256; // Maximum vector size (tunable)

  // Minimum number of tracks in a basket that stay in the same volume and
  // therefore can be re-used with the same thread without re-basketizing.
  config->fNminReuse = 10000; // (default in propagator is 4)

  // Kill threshold - number of steps allowed before killing a track
  //                  (includes internal geometry steps)
  config->fNstepsKillThr = 100000;

  config->fEmin = 3.E-6; // [3 KeV] energy cut
  config->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy

  // Number of steps for learning phase (tunable [0, 1e6])
  // if set to 0 disable learning phase
  config->fLearnSteps                  = 0;
  if (performance) config->fLearnSteps = 0;

  // Activate I/O
  config->fFillTree = false;
  // Activate old version of single thread serialization/reading
  // config->fConcurrentWrite = false;

  // Activate debugging using -DBUG_HUNT=ON in your cmake build
  config->fDebugEvt = 0;
  config->fDebugTrk = 0;
  config->fDebugStp = 0;
  config->fDebugRep = 10;

  // Activate standard scoring
  config->fUseStdScoring                  = true;
  if (performance) config->fUseStdScoring = false;
  // Monitor the application
  config->fUseAppMonitoring = false;

  RunManager *runMgr = new RunManager(npropagators, nthreads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);

  // Create the tab. phys process.
  runMgr->SetPhysicsProcess(new geant::TTabPhysProcess("tab_phys", xsec, fstate));

  // for vector physics -OFF now
  // config->fVectorPhysicsProcess = new GVectorPhysicsProcess(config->fEmin, nthreads);

  runMgr->SetPrimaryGenerator(new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0));

  runMgr->SetUserApplication(new ExN03Application(runMgr));

  runMgr->RunSimulation();
  delete config;
}
