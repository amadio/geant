#ifdef USE_ROOT
#include "Rtypes.h"
#include "TGeoManager.h"
#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

#include <err.h>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "GeantRunManager.h"
#include "GunGenerator.h"
#include "TaskBroker.h"
#include "TTabPhysProcess.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "ExN03Application.h"
#include "ExN03DetectorConstruction.h"
#include "ExternalFramework.h"

using namespace Geant;

static int n_events = 50;
static int n_buffered = 16;
static int n_threads = 4;
static int n_track_max = 500;
static int n_learn_steps = 0;
static int n_reuse = 100000;
static int n_propagators = 1;
static bool monitor = false, score = false, debug = false, coprocessor = false;
static bool external_loop = false, usev3 = true, usenuma = false;

static struct option options[] = {{"events", required_argument, 0, 'e'},
                                  {"fstate", required_argument, 0, 'f'},
                                  {"geometry", required_argument, 0, 'g'},
                                  {"learn-steps", required_argument, 0, 'l'},
                                  {"max-tracks-per-basket", required_argument, 0, 'B'},
                                  {"monitor", no_argument, 0, 'm'},
                                  {"debug", no_argument, 0, 'd'},
                                  {"nbuffer", required_argument, 0, 'b'},
                                  {"score", no_argument, 0, 's'},
                                  {"threads", required_argument, 0, 't'},
                                  {"xsec", required_argument, 0, 'x'},
                                  {"coprocessor", required_argument, 0, 'r'},
                                  {"external_loop", required_argument, 0, 'i'},
                                  {"reuse", required_argument, 0, 'u'},
                                  {"propagators", required_argument, 0, 'p'},
                                  {"v3", required_argument, 0, 'v'},
                                  {"numa", required_argument, 0, 'n'},
                                  {0, 0, 0, 0}};

void help() {
  printf("\nUsage: cmsapp [OPTIONS] INPUT_FILE\n\n");

  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\t%s\n", options[i].val, options[i].name, options[i].has_arg ? options[i].name : "");
  }
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  std::cout << "Avoid ctest truncation of output: CTEST_FULL_OUTPUT" << std::endl;
  std::string exn03_geometry_filename("ExN03.root");
  std::string xsec_filename("xsec_FTFP_BERT.root");
  std::string fstate_filename("fstate_FTFP_BERT.root");

  if (argc == 1) {
    help();
    exit(0);
  }

  while (true) {
    int c, optidx = 0;

    c = getopt_long(argc, argv, "e:f:g:l:B:m:b:t:x:r:i:u:p:v:n:", options, &optidx);

    if (c == -1)
      break;

    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'e':
      n_events = (int)strtol(optarg, NULL, 10);

      if (n_events <= 0)
        errx(1, "number of events must be positive");
      break;

    case 'f':
      fstate_filename = optarg;
      break;

    case 'g':
      exn03_geometry_filename = optarg;
      break;

    case 'l':
      n_learn_steps = (int)strtol(optarg, NULL, 10);

      if (n_learn_steps <= 0)
        errx(1, "number of learning steps must be positive");

      break;

    case 'B':
      n_track_max = (int)strtol(optarg, NULL, 10);

      if (n_track_max < 1)
        errx(1, "max number of tracks per basket must be positive");

      break;

    case 'm':
      monitor = true;
      break;

    case 'b':
      n_buffered = (int)strtol(optarg, NULL, 10);

      if (n_buffered < 1)
        errx(1, "number of buffered events must be positive");
      break;

    case 't':
      n_threads = (int)strtol(optarg, NULL, 10);

      if (n_threads < 1)
        errx(1, "number of threads must be positive");
      break;

    case 's':
      score = true;
      break;

    case 'x':
      xsec_filename = optarg;
      break;

    case 'r':
      coprocessor = optarg;
      break;

    case 'i':
      external_loop = true;
      break;

    case 'u':
      n_reuse = (int)strtol(optarg, NULL, 10);
      break;

    case 'p':
      n_propagators = (int)strtol(optarg, NULL, 10);
      break;

    case 'v':
      usev3 = bool(strtol(optarg, NULL, 10));
      if (!usev3) n_buffered = 4;
      break;

    case 'n':
      usenuma = bool(strtol(optarg, NULL, 10));
      break;

    default:
      errx(1, "unknown option %c", c);
    }
  }
  bool performance = true;

  TaskBroker *broker = nullptr;
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32,128,1);
    broker = gpuBroker;
    nthreads += gpuBroker->GetNstream()+1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  GeantConfig* config=new GeantConfig();

  
//  TGeoManager::Import(exn03_geometry_filename.c_str());
  config->fGeomFileName = exn03_geometry_filename;
  config->fNtotal = n_events;
  config->fNbuff = n_buffered;
  config->fBmag = 1.; // 0.1 Tesla
  // V3 options
  config->fNstackLanes = 10;
  config->fNmaxBuffSpill = 128;  // New configuration parameter!!!
  config->fUseV3 = usev3;
  config->fUseNuma = usenuma;
  config->fUseMonitoring = monitor;
  config->fNminThreshold=5*n_threads;
  config->SetMonitored(GeantConfig::kMonQueue, false);
  config->SetMonitored(GeantConfig::kMonMemory, monitor);
  config->SetMonitored(GeantConfig::kMonBasketsPerVol, false);
  config->SetMonitored(GeantConfig::kMonVectors, false);
  config->SetMonitored(GeantConfig::kMonConcurrency, false);
  config->SetMonitored(GeantConfig::kMonTracksPerEvent, false);
  config->fNaverage = 500;   // Average number of tracks per event
  
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = 0.05;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = 16;   // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = n_track_max;   // Maximum vector size (tunable)
  config->fEmin = 3.E-6; // [3 KeV] energy cut
  config->fEmax = 0.3;  // [300MeV] used for now to select particle gun energy

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
  config->fLearnSteps = n_learn_steps;
  if (performance) config->fLearnSteps = 0;
   // Activate I/O
  config->fFillTree = false;
   // Activate debugging using -DBUG_HUNT=ON in your cmake build
  if (debug) {
    config->fUseDebug = true;
    config->fDebugTrk = 1;
  }
// Activate standard scoring   
  config->fUseStdScoring = true;
  if (performance) config->fUseStdScoring = false;
  // Monitor the application
  config->fUseAppMonitoring = false;

  // Set threshold for tracks to be reused in the same volume
  config->fNminReuse = n_reuse;

  // Create run manager
  GeantRunManager *runMgr = new GeantRunManager(n_propagators, n_threads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);
  // Create the tab. phys process.
  runMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str()));
  
// Create the tab. phys process.
#ifdef USE_VECGEOM_NAVIGATOR
//  runMgr->LoadVecGeomGeometry();
#endif

  // for vector physics -OFF now
  // runMgr->SetVectorPhysicsProcess(new GVectorPhysicsProcess(config->fEmin, nthreads));

  GunGenerator *generator = new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0);
  runMgr->SetPrimaryGenerator( generator );  
  runMgr->SetUserApplication ( new ExN03Application(runMgr) );
  runMgr->SetDetectorConstruction( new ExN03DetectorConstruction(exn03_geometry_filename.c_str(), runMgr) );
  if (external_loop) {
    userfw::Framework fw(n_propagators*n_threads, n_events, runMgr, generator);
    fw.Run();
  } else {
    runMgr->RunSimulation();
  }
  return 0;
}
#endif
