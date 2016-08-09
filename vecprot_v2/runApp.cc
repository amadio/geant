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

#include "GunGenerator.h"
#include "TaskBroker.h"
#include "TTabPhysProcess.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "ExN03Application.h"

#ifdef GEANT_TBB
#include "TaskMgrTBB.h"
#endif

static int n_events = 50;
static int n_buffered = 10;
static int n_threads = 4;
static int n_track_max = 500;
static int n_learn_steps = 0;
static int n_reuse = 100000;
static bool monitor = false, score = false, debug = false, coprocessor = false, tbbmode = false;

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
                                  {"tbbmode", required_argument, 0, 'i'},
                                  {"reuse", required_argument, 0, 'u'},
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

    c = getopt_long(argc, argv, "e:f:g:l:B:b:t:x:r:i:u:", options, &optidx);

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
      tbbmode = true;
      break;

    case 'u':
      n_reuse = (int)strtol(optarg, NULL, 10);
      break;

    default:
      errx(1, "unknown option %c", c);
    }
  }
  bool performance = true;
  TGeoManager::Import(exn03_geometry_filename.c_str());
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
  GeantPropagator *propagator = GeantPropagator::NewInstance(n_events, n_buffered,n_threads);
  WorkloadManager *wmanager = propagator->WorkloadManager();

  if (broker) propagator->SetTaskBroker(broker);
  wmanager->SetNminThreshold(5 * n_threads);
  propagator->fUseMonitoring = monitor;
#ifdef GEANT_TBB
  if (tbbmode)
    propagator->fTaskMgr = new TaskMgrTBB();
#endif

  // Monitor different features
  wmanager->SetNminThreshold(5*n_threads);
  wmanager->SetMonitored(GeantPropagator::kMonQueue, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonMemory, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonBasketsPerVol, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonVectors, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonConcurrency, monitor);
  wmanager->SetMonitored(GeantPropagator::kMonTracksPerEvent, monitor);
  propagator->fUseMonitoring = monitor;
  propagator->fNaverage = 500;   // Average number of tracks per event

  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  propagator->fPriorityThr = 0.05;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  propagator->fNperBasket = 16;   // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  propagator->fMaxPerBasket = n_track_max;   // Maximum vector size (tunable)
  propagator->fEmin = 3.E-6; // [3 KeV] energy cut
  propagator->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy

  // Create the tab. phys process.
  propagator->fProcess = new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str());
// Create the tab. phys process.
#ifdef USE_VECGEOM_NAVIGATOR
  propagator->LoadVecGeomGeometry();
#endif

  // for vector physics -OFF now
  // propagator->fVectorPhysicsProcess = new GVectorPhysicsProcess(propagator->fEmin, nthreads);
  propagator->fPrimaryGenerator = new GunGenerator(propagator->fNaverage, 11, propagator->fEmax, -8, 0, 0, 1, 0, 0);

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
  propagator->fLearnSteps = n_learn_steps;
  if (performance) propagator->fLearnSteps = 0;
  propagator->fApplication = new ExN03Application(propagator);
   // Activate I/O
  propagator->fFillTree = false;
  
  // Set threshold for tracks to be reused in the same volume
  propagator->fNminReuse = n_reuse;
  
   // Activate debugging using -DBUG_HUNT=ON in your cmake build
  if (debug) {
    propagator->fUseDebug = true;
    propagator->fDebugTrk = 1;
  }
// Activate standard scoring
  propagator->fUseStdScoring = true;
  if (performance) propagator->fUseStdScoring = false;
  // Monitor the application
  propagator->fUseAppMonitoring = false;

  propagator->PropagatorGeom(exn03_geometry_filename.c_str(), n_threads, monitor);
  return 0;
}
#endif
